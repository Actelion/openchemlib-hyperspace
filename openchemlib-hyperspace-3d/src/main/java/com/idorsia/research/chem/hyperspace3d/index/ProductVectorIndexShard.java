package com.idorsia.research.chem.hyperspace3d.index;

/** The three independently checksummed columns of one version-2 product-vector shard. */
public record ProductVectorIndexShard(
        String vectorsPath, String rowsPath, String tuplesPath, long recordCount,
        String vectorsSha256, String rowsSha256, String tuplesSha256) {}
