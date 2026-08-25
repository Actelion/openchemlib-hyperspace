package com.idorsia.research.chem.hyperspace3d.index;

/** Vector columns shared by Java indexes and Deepspace7 supplier caches. */
public enum MoleculeFingerprintColumn {
    BASE_128(128), COMPACT_16(16);

    private final int dimension;

    MoleculeFingerprintColumn(int dimension) { this.dimension = dimension; }
    public int dimension() { return dimension; }
}
