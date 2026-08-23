package com.idorsia.research.chem.hyperspace3d.feature;

import java.util.Arrays;

public final class DeepSpaceTensorBatch {
    private final int batchSize;
    private final float[] atomFeatures;
    private final float[] pairFeatures;
    private final boolean[] atomMask;

    public DeepSpaceTensorBatch(int batchSize, float[] atomFeatures,
                                float[] pairFeatures, boolean[] atomMask) {
        int atoms = DeepSpaceFeatureSchema.MAX_ATOMS;
        if (batchSize < 1
                || atomFeatures.length != batchSize * atoms * DeepSpaceFeatureSchema.ATOM_FEATURE_DIM
                || pairFeatures.length != batchSize * atoms * atoms * DeepSpaceFeatureSchema.PAIR_FEATURE_DIM
                || atomMask.length != batchSize * atoms) {
            throw new IllegalArgumentException("tensor dimensions do not match the 50/36/32 contract");
        }
        this.batchSize = batchSize;
        this.atomFeatures = atomFeatures;
        this.pairFeatures = pairFeatures;
        this.atomMask = atomMask;
    }

    public int batchSize() { return batchSize; }
    public float[] atomFeatures() { return atomFeatures; }
    public float[] pairFeatures() { return pairFeatures; }
    public boolean[] atomMask() { return atomMask; }

    public DeepSpaceTensorBatch copy() {
        return new DeepSpaceTensorBatch(batchSize, Arrays.copyOf(atomFeatures, atomFeatures.length),
                Arrays.copyOf(pairFeatures, pairFeatures.length), Arrays.copyOf(atomMask, atomMask.length));
    }
}
