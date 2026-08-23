package com.idorsia.research.chem.hyperspace3d.index;

/** Stable address of a product in a columnar index. */
public record ProductVectorReference(int shardIndex, long rowIndex, int reactionOrdinal,
                                     long tupleOffset) {}
