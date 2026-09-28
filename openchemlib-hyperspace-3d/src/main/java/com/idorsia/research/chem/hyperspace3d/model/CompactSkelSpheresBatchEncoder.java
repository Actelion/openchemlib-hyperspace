package com.idorsia.research.chem.hyperspace3d.model;

import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;

public final class CompactSkelSpheresBatchEncoder implements ProductEmbeddingBatchEncoder {
    private final ProductEmbeddingBatchEncoder baseEncoder;
    private final CompactSkelSpheresProjector projector;

    public CompactSkelSpheresBatchEncoder(
            ProductEmbeddingBatchEncoder baseEncoder,
            CompactSkelSpheresProjector projector) {
        this.baseEncoder = baseEncoder;
        this.projector = projector;
    }

    @Override
    public float[][] encode(DeepSpaceTensorBatch batch) {
        return projector.project(baseEncoder.encode(batch));
    }
}
