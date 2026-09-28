package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.OnnxTensor;
import ai.onnxruntime.OrtException;
import ai.onnxruntime.OrtSession;
import java.nio.FloatBuffer;
import java.util.Map;
import java.util.Set;

public final class CompactSkelSpheresProjector implements AutoCloseable {
    private final DeepSpaceOnnxEnvironment runtime;
    private final OrtSession session;

    public CompactSkelSpheresProjector(
            DeepSpaceOnnxEnvironment runtime, CompactSkelSpheresModelBundle bundle) {
        this.runtime = runtime;
        this.session = runtime.open(bundle.projectionPath());
        if (!session.getInputNames().equals(Set.of("base_embedding"))
                || !session.getOutputNames().equals(Set.of("compact_embedding"))) {
            throw new DeepSpaceInferenceException("compact projection tensor names mismatch");
        }
        OnnxContract.requireInput(session, "base_embedding",
                ai.onnxruntime.OnnxJavaType.FLOAT, 128);
        OnnxContract.requireOutput(session, "compact_embedding",
                ai.onnxruntime.OnnxJavaType.FLOAT, 16);
    }

    public synchronized float[][] project(float[][] baseEmbeddings) {
        if (baseEmbeddings.length < 1) {
            throw new IllegalArgumentException("projection batch must not be empty");
        }
        float[] flattened = new float[baseEmbeddings.length * 128];
        for (int row = 0; row < baseEmbeddings.length; row++) {
            if (baseEmbeddings[row].length != 128) {
                throw new IllegalArgumentException("base embedding must be 128D");
            }
            System.arraycopy(baseEmbeddings[row], 0, flattened, row * 128, 128);
        }
        try (OnnxTensor input = OnnxTensor.createTensor(runtime.environment(),
                     FloatBuffer.wrap(flattened), new long[]{baseEmbeddings.length, 128});
             OrtSession.Result result = session.run(Map.of("base_embedding", input))) {
            float[][] projected = (float[][]) result.get("compact_embedding")
                    .orElseThrow(() -> new DeepSpaceInferenceException(
                            "compact projection output missing"))
                    .getValue();
            if (projected.length != baseEmbeddings.length
                    || projected[0].length != 16) {
                throw new DeepSpaceInferenceException(
                        "compact projection returned an incompatible shape");
            }
            return projected;
        } catch (OrtException error) {
            throw new DeepSpaceInferenceException("compact projection failed", error);
        }
    }

    public float[] project(float[] baseEmbedding) {
        return project(new float[][]{baseEmbedding})[0];
    }

    @Override
    public void close() {
        try {
            session.close();
        } catch (OrtException error) {
            throw new DeepSpaceInferenceException(
                    "cannot close compact projection session", error);
        }
    }
}
