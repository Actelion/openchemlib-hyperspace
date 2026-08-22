package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.OnnxTensor;
import ai.onnxruntime.OrtException;
import ai.onnxruntime.OrtSession;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceFeatureSchema;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;
import java.nio.FloatBuffer;
import java.util.Map;
import java.util.Set;

public final class DeepSpaceV1Encoder implements ProductEmbeddingBatchEncoder, AutoCloseable {
    private final DeepSpaceOnnxEnvironment runtime;
    private final OrtSession session;

    public DeepSpaceV1Encoder(DeepSpaceOnnxEnvironment runtime, DeepSpaceModelBundle bundle) {
        this.runtime = runtime;
        this.session = runtime.open(bundle.encoderPath());
        validateNames(session.getInputNames(), Set.of("atom_x", "pair_x", "atom_mask"), "encoder inputs");
        validateNames(session.getOutputNames(), Set.of("embedding"), "encoder outputs");
        OnnxContract.requireInput(session, "atom_x", ai.onnxruntime.OnnxJavaType.FLOAT, 32, 56);
        OnnxContract.requireInput(session, "pair_x", ai.onnxruntime.OnnxJavaType.FLOAT, 32, 32, 36);
        OnnxContract.requireInput(session, "atom_mask", ai.onnxruntime.OnnxJavaType.BOOL, 32);
        OnnxContract.requireOutput(session, "embedding", ai.onnxruntime.OnnxJavaType.FLOAT, 128);
    }

    @Override
    public synchronized float[][] encode(DeepSpaceTensorBatch batch) {
        int b = batch.batchSize();
        int n = DeepSpaceFeatureSchema.MAX_ATOMS;
        boolean[][] mask = new boolean[b][n];
        for (int i = 0; i < b; i++) {
            System.arraycopy(batch.atomMask(), i * n, mask[i], 0, n);
        }
        try (OnnxTensor atoms = OnnxTensor.createTensor(runtime.environment(),
                     FloatBuffer.wrap(batch.atomFeatures()), new long[]{b, n, 56});
             OnnxTensor pairs = OnnxTensor.createTensor(runtime.environment(),
                     FloatBuffer.wrap(batch.pairFeatures()), new long[]{b, n, n, 36});
             OnnxTensor atomMask = OnnxTensor.createTensor(runtime.environment(), mask);
             OrtSession.Result result = session.run(
                     Map.of("atom_x", atoms, "pair_x", pairs, "atom_mask", atomMask))) {
            float[][] embeddings = (float[][]) result.get("embedding")
                    .orElseThrow(() -> new DeepSpaceInferenceException("encoder output missing"))
                    .getValue();
            if (embeddings.length != b || embeddings[0].length != 128) {
                throw new DeepSpaceInferenceException("encoder returned an incompatible shape");
            }
            return embeddings;
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("encoder inference failed", e);
        }
    }

    private static void validateNames(Set<String> observed, Set<String> expected, String label) {
        if (!observed.equals(expected)) {
            throw new DeepSpaceInferenceException(label + " mismatch: " + observed);
        }
    }

    @Override public void close() {
        try { session.close(); } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot close encoder session", e);
        }
    }
}
