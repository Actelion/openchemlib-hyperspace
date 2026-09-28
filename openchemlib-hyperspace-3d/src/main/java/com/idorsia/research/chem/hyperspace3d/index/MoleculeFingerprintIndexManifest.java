package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;
import java.util.Map;

/** Manifest for the flat, dual-vector molecule index. It intentionally has no hashes. */
public final class MoleculeFingerprintIndexManifest {
    public String artifactType = "hyperspace-molecule-fingerprint-index";
    public int formatVersion = 1;
    public String input;
    public String smilesColumn;
    public String idColumn;
    public String modelBundle;
    public String compactBundle;
    public String vectorDtype = "float16";
    public int baseDimension = 128;
    public int compactDimension = 16;
    public long sourceRowsPerShard;
    public long sourceRowCount;
    public long recordCount;
    public long rejectedCount;
    public List<MoleculeFingerprintIndexShard> shards;
    public Map<String, Object> runtime;
    public Map<String, Object> buildStatistics;

    public void validate() {
        require("hyperspace-molecule-fingerprint-index".equals(artifactType),
                "unsupported molecule-index artifact");
        require(formatVersion == 1, "unsupported molecule-index format version");
        require(text(input) && text(smilesColumn) && text(idColumn)
                && text(modelBundle) && text(compactBundle), "incomplete index description");
        require("float16".equals(vectorDtype) && baseDimension == 128
                && compactDimension == 16, "unsupported vector representation");
        require(sourceRowsPerShard > 0 && sourceRowsPerShard <= 8_000_000,
                "invalid source rows per shard");
        require(sourceRowCount >= 0 && recordCount >= 0 && rejectedCount >= 0
                && recordCount + rejectedCount == sourceRowCount,
                "inconsistent molecule-index counts");
        require(shards != null, "missing shard metadata");
        long sources = 0;
        long records = 0;
        long expectedFirst = 0;
        for (int i = 0; i < shards.size(); i++) {
            MoleculeFingerprintIndexShard shard = shards.get(i);
            require(shard != null && shard.shardIndex == i
                    && text(shard.directory) && shard.firstSourceRow == expectedFirst
                    && shard.sourceRowCount > 0 && shard.recordCount >= 0
                    && shard.rejectedCount >= 0
                    && shard.recordCount + shard.rejectedCount == shard.sourceRowCount,
                    "invalid molecule-index shard " + i);
            sources += shard.sourceRowCount;
            records += shard.recordCount;
            expectedFirst += shard.sourceRowCount;
        }
        require(sources == sourceRowCount && records == recordCount,
                "molecule-index shard totals do not match");
    }

    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }
}
