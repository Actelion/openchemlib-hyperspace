package com.idorsia.research.chem.hyperspace3d.index;

import java.util.LinkedHashMap;
import java.util.Map;

/** Lightweight metadata for one atomically committed source-row shard. */
public final class MoleculeFingerprintIndexShard {
    public int shardIndex;
    public String directory;
    public long firstSourceRow;
    public long sourceRowCount;
    public long recordCount;
    public long rejectedCount;
    public Map<String, Long> rejections = new LinkedHashMap<>();
}
