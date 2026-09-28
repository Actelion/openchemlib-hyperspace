package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;

/** Flat row-major fingerprint batch without per-molecule vector allocations. */
public record MoleculeVectorBatch(float[] values, int dimension,
        List<MoleculeVectorReference> references) {
    public MoleculeVectorBatch {
        if (dimension < 1 || values.length != Math.multiplyExact(references.size(), dimension)) {
            throw new IllegalArgumentException("inconsistent molecule vector batch");
        }
        references = List.copyOf(references);
    }

    public int size() { return references.size(); }

    public float[] vector(int row) {
        if (row < 0 || row >= size()) throw new IndexOutOfBoundsException(row);
        float[] result = new float[dimension];
        System.arraycopy(values, row * dimension, result, 0, dimension);
        return result;
    }
}
