package com.idorsia.research.chem.hyperspace3d.index;

/** Training-checkpoint identity attached to cached molecule fingerprints. */
public record MoleculeFingerprintProvenance(String foundationCheckpointSha256,
        String predictorCheckpointSha256, String projectionCheckpointSha256) {}
