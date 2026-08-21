package com.idorsia.research.chem.hyperspace3d;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;
import java.util.List;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class DeepSpaceModelManifestTest {
    @Test void acceptsSelectedV1V3Contract() { valid().validate(); }

    @Test void rejectsFiftyDimensionalAndPoseBundles() {
        DeepSpaceModelManifest fifty = valid();
        fifty.atomFeatureDim = 50;
        assertThrows(IllegalArgumentException.class, fifty::validate);
        DeepSpaceModelManifest pose = valid();
        pose.architecture = "pose_conditioned_v2";
        assertThrows(IllegalArgumentException.class, pose::validate);
    }

    static DeepSpaceModelManifest valid() {
        DeepSpaceModelManifest m = new DeepSpaceModelManifest();
        m.artifactType = "deepspace7-v1-onnx-bundle";
        m.formatVersion = 1;
        m.modelVersion = "deepspace7-v1-v3graph";
        m.architecture = "compact_graph_v1";
        m.maxAtoms = 32; m.atomFeatureDim = 56; m.pairFeatureDim = 36;
        m.nodeLatentDim = 16; m.pairLatentDim = 14; m.embeddingDim = 128;
        m.tokenHiddenDim = 64; m.pairComparisonHiddenDim = 128;
        m.graphFeatureContract = "deepspace7-v3";
        m.pharmacophoreAtomFeatures = true; m.representation = "node_pair";
        m.availableTargets = List.of("ffp_similarity", "skelspheres_similarity",
                "flexophore_similarity", "phesa_total", "phesa_pharmacophore", "phesa_shape");
        m.phesaPpWeight = 0.5;
        m.encoderSha256 = "a".repeat(64); m.comparatorSha256 = "b".repeat(64);
        m.featureSchemaSha256 = "c".repeat(64);
        return m;
    }
}
