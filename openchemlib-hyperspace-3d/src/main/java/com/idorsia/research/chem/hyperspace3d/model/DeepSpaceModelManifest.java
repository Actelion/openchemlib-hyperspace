package com.idorsia.research.chem.hyperspace3d.model;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import java.util.List;
import java.util.Objects;

@JsonIgnoreProperties(ignoreUnknown = true)
public final class DeepSpaceModelManifest {
    public String artifactType;
    public int formatVersion;
    public String modelVersion;
    public String architecture;
    public int maxAtoms;
    public int atomFeatureDim;
    public int pairFeatureDim;
    public int nodeLatentDim;
    public int pairLatentDim;
    public int embeddingDim;
    public int tokenHiddenDim;
    public int pairComparisonHiddenDim;
    public String graphFeatureContract;
    public boolean pharmacophoreAtomFeatures;
    public String representation;
    public List<String> availableTargets;
    public Double phesaPpWeight;
    public String encoderSha256;
    public String comparatorSha256;
    public String featureSchemaSha256;
    public String foundationCheckpointSha256;
    public String predictorCheckpointSha256;

    public void validate() {
        require("deepspace7-v1-onnx-bundle".equals(artifactType), "unsupported artifactType");
        require(formatVersion == 1, "unsupported formatVersion");
        require("compact_graph_v1".equals(architecture), "unsupported architecture");
        require("deepspace7-v3".equals(graphFeatureContract), "unsupported graphFeatureContract");
        require(maxAtoms == 32 && atomFeatureDim == 56 && pairFeatureDim == 36,
                "incompatible graph tensor dimensions");
        require(nodeLatentDim == 16 && pairLatentDim == 14 && embeddingDim == 128,
                "incompatible latent dimensions");
        require(tokenHiddenDim == 64 && pairComparisonHiddenDim == 128,
                "incompatible compact V1 hidden dimensions");
        require(pharmacophoreAtomFeatures, "V3 pharmacophore inputs must be enabled");
        require("node_pair".equals(representation), "unsupported representation");
        List<String> expected = List.of("ffp_similarity", "skelspheres_similarity",
                "flexophore_similarity", "phesa_total", "phesa_pharmacophore", "phesa_shape");
        require(expected.equals(availableTargets), "target ordering mismatch");
        require(phesaPpWeight != null && Double.isFinite(phesaPpWeight), "missing phesaPpWeight");
        requireHash(encoderSha256, "encoderSha256");
        requireHash(comparatorSha256, "comparatorSha256");
        requireHash(featureSchemaSha256, "featureSchemaSha256");
        Objects.requireNonNull(modelVersion, "modelVersion");
    }

    public int targetIndex(String target) {
        int index = availableTargets.indexOf(target);
        if (index < 0) throw new IllegalArgumentException("unknown target: " + target);
        return index;
    }

    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    private static void requireHash(String value, String name) {
        require(value != null && value.matches("[0-9a-f]{64}"), "missing or invalid " + name);
    }
}
