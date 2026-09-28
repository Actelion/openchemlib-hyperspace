package com.idorsia.research.chem.hyperspace3d.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import java.util.Map;
import java.util.Objects;

@JsonIgnoreProperties(ignoreUnknown = true)
public final class CompactSkelSpheresManifest {
    public String artifactType;
    public int formatVersion;
    public String modelVersion;
    public String architecture;
    public String target;
    public int inputEmbeddingDim;
    public int hiddenDim;
    public int outputEmbeddingDim;
    public boolean l2Normalized;
    public int canonicalSeed;
    public int bestEpoch;
    public String recommendedStorageDtype;
    public Calibration calibration;
    public Map<String, Object> qualityGate;
    public String sourceFoundationCheckpointSha256;
    public String sourcePredictorCheckpointSha256;
    public String projectionCheckpointSha256;
    public String projectionSha256;
    public int onnxOpset;

    public void validate() {
        require("deepspace7-skelspheres16-onnx-bundle".equals(artifactType),
                "unsupported compact artifactType");
        require(formatVersion == 1, "unsupported compact formatVersion");
        require("mlp_128_128_16_l2".equals(architecture),
                "unsupported compact architecture");
        require("skelspheres_similarity".equals(target), "unsupported compact target");
        require(inputEmbeddingDim == 128 && hiddenDim == 128 && outputEmbeddingDim == 16,
                "incompatible compact dimensions");
        require(l2Normalized, "compact embeddings must be L2 normalized");
        require(canonicalSeed == 17, "compact model must be canonical seed 17");
        require(bestEpoch >= 0, "missing compact best epoch");
        require("float16".equals(recommendedStorageDtype),
                "unsupported compact storage recommendation");
        require(calibration != null
                        && "sigmoid_affine_cosine".equals(calibration.type)
                        && Double.isFinite(calibration.scale)
                        && calibration.scale > 0
                        && Double.isFinite(calibration.bias),
                "invalid compact score calibration");
        requireHash(sourceFoundationCheckpointSha256, "sourceFoundationCheckpointSha256");
        requireHash(sourcePredictorCheckpointSha256, "sourcePredictorCheckpointSha256");
        requireHash(projectionCheckpointSha256, "projectionCheckpointSha256");
        requireHash(projectionSha256, "projectionSha256");
        Objects.requireNonNull(modelVersion, "modelVersion");
    }

    public void validateSource(DeepSpaceModelManifest source) {
        validate();
        if (!sourceFoundationCheckpointSha256.equals(source.foundationCheckpointSha256)
                || !sourcePredictorCheckpointSha256.equals(source.predictorCheckpointSha256)) {
            throw new IllegalArgumentException(
                    "compact model was trained from a different primary model bundle");
        }
    }

    public double calibrate(double dotProduct) {
        return 1.0 / (1.0 + Math.exp(
                -(calibration.scale * dotProduct + calibration.bias)));
    }

    @JsonIgnoreProperties(ignoreUnknown = true)
    public static final class Calibration {
        public String type;
        public double scale;
        public double bias;
    }

    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    private static void requireHash(String value, String name) {
        require(value != null && value.matches("[0-9a-f]{64}"),
                "missing or invalid " + name);
    }
}
