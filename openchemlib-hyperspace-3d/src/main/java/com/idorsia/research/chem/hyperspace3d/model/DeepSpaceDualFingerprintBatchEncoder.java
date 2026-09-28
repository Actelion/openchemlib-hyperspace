package com.idorsia.research.chem.hyperspace3d.model;

import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;
import java.util.Objects;

/** Runs the expensive base encoder once and projects its output to the compact space. */
public final class DeepSpaceDualFingerprintBatchEncoder
        implements MoleculeFingerprintBatchEncoder {
    private final ProductEmbeddingBatchEncoder baseEncoder;
    private final CompactSkelSpheresProjector compactProjector;

    public DeepSpaceDualFingerprintBatchEncoder(ProductEmbeddingBatchEncoder baseEncoder,
            CompactSkelSpheresProjector compactProjector) {
        this.baseEncoder = Objects.requireNonNull(baseEncoder);
        this.compactProjector = Objects.requireNonNull(compactProjector);
    }

    @Override public MoleculeFingerprintBatch encode(DeepSpaceTensorBatch batch) {
        float[][] base = baseEncoder.encode(batch);
        return new MoleculeFingerprintBatch(base, compactProjector.project(base));
    }
}
