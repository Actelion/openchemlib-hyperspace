package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.OnnxTensor;
import ai.onnxruntime.OrtException;
import ai.onnxruntime.OrtSession;
import java.nio.FloatBuffer;
import java.util.Map;
import java.util.Set;

public final class DeepSpaceV1Comparator implements EmbeddingComparator, AutoCloseable {
    private final DeepSpaceOnnxEnvironment runtime;
    private final OrtSession session;
    private final int targets;

    public DeepSpaceV1Comparator(DeepSpaceOnnxEnvironment runtime, DeepSpaceModelBundle bundle) {
        this.runtime = runtime;
        this.session = runtime.open(bundle.comparatorPath());
        this.targets = bundle.manifest().availableTargets.size();
        if (!session.getInputNames().equals(Set.of("query_embedding", "candidate_embedding"))
                || !session.getOutputNames().equals(Set.of("scores"))) {
            throw new DeepSpaceInferenceException("comparator ONNX names do not match the manifest contract");
        }
        OnnxContract.requireInput(session, "query_embedding", ai.onnxruntime.OnnxJavaType.FLOAT, 128);
        OnnxContract.requireInput(session, "candidate_embedding", ai.onnxruntime.OnnxJavaType.FLOAT, 128);
        OnnxContract.requireOutput(session, "scores", ai.onnxruntime.OnnxJavaType.FLOAT, targets);
    }

    @Override
    public synchronized float[][] compare(float[] queryEmbedding, float[][] candidateEmbeddings) {
        if (queryEmbedding.length != 128) throw new IllegalArgumentException("query embedding must be 128D");
        if (candidateEmbeddings.length == 0) return new float[0][targets];
        int b = candidateEmbeddings.length;
        float[] queries = new float[b * 128];
        float[] candidates = new float[b * 128];
        for (int i = 0; i < b; i++) {
            if (candidateEmbeddings[i].length != 128) {
                throw new IllegalArgumentException("candidate embedding " + i + " is not 128D");
            }
            System.arraycopy(queryEmbedding, 0, queries, i * 128, 128);
            System.arraycopy(candidateEmbeddings[i], 0, candidates, i * 128, 128);
        }
        try (OnnxTensor query = OnnxTensor.createTensor(runtime.environment(),
                     FloatBuffer.wrap(queries), new long[]{b, 128});
             OnnxTensor candidate = OnnxTensor.createTensor(runtime.environment(),
                     FloatBuffer.wrap(candidates), new long[]{b, 128});
             OrtSession.Result result = session.run(
                     Map.of("query_embedding", query, "candidate_embedding", candidate))) {
            float[][] scores = (float[][]) result.get("scores")
                    .orElseThrow(() -> new DeepSpaceInferenceException("comparator output missing"))
                    .getValue();
            if (scores.length != b || scores[0].length != targets) {
                throw new DeepSpaceInferenceException("comparator returned an incompatible shape");
            }
            return scores;
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("comparator inference failed", e);
        }
    }

    @Override public void close() {
        try { session.close(); } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot close comparator session", e);
        }
    }
}
