package com.idorsia.research.chem.hyperspace3d.index;

/** Stable position of a molecule inside a flat fingerprint index. */
public record MoleculeVectorReference(int shardIndex, long localRow) {}
