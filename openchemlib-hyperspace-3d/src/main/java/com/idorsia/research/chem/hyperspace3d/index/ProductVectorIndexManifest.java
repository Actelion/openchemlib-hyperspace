package com.idorsia.research.chem.hyperspace3d.index;

import java.util.HashSet;
import java.util.List;
import java.util.Map;

/** Columnar format-v2 manifest. Version 1 remains readable through ProductFingerprintIndexManifest. */
public final class ProductVectorIndexManifest {
    public String artifactType = "hyperspace-product-vector-index";
    public int formatVersion = 2;
    public String rawspaceIdentity;
    public String rawspaceSha256;
    public String downsampledRawspaceIdentity;
    public String downsampledRawspaceSha256;
    public String sourceModelBundleHash;
    public String projectionModelBundleHash;
    public String modelVersion;
    public String target = "skelspheres_similarity";
    public String samplingAlgorithm;
    public long samplingSeed;
    public Map<String, Object> reactionWeightConfiguration;
    public Map<String, Object> molecularFilters;
    public String vectorDtype = "float16";
    public int embeddingDimension = 16;
    public long recordCount;
    public List<String> reactionDictionary;
    public List<ProductVectorIndexShard> shards;
    public String buildConfigurationHash;
    public Map<String, Object> buildStatistics;

    public void validate() {
        require("hyperspace-product-vector-index".equals(artifactType), "unsupported vector-index artifact");
        require(formatVersion == 2, "unsupported vector-index format version");
        require("skelspheres_similarity".equals(target), "unsupported vector-index target");
        require("float16".equals(vectorDtype) && embeddingDimension == 16,
                "this reader requires 16D float16 vectors");
        require(text(rawspaceIdentity) && hash(rawspaceSha256)
                && text(downsampledRawspaceIdentity) && hash(downsampledRawspaceSha256)
                && hash(sourceModelBundleHash) && hash(projectionModelBundleHash)
                && text(modelVersion) && text(samplingAlgorithm), "incomplete index provenance");
        require(reactionWeightConfiguration != null && molecularFilters != null,
                "missing sampling or filter provenance");
        require(reactionDictionary != null && !reactionDictionary.isEmpty()
                && new HashSet<>(reactionDictionary).size() == reactionDictionary.size(),
                "invalid reaction dictionary");
        require(shards != null && recordCount >= 0, "invalid shard metadata");
        long sum = 0;
        for (ProductVectorIndexShard shard : shards) {
            require(shard != null && text(shard.vectorsPath()) && text(shard.rowsPath())
                    && text(shard.tuplesPath()) && shard.recordCount() > 0
                    && hash(shard.vectorsSha256()) && hash(shard.rowsSha256())
                    && hash(shard.tuplesSha256()), "invalid vector-index shard");
            sum += shard.recordCount();
        }
        require(sum == recordCount, "shard record count mismatch");
    }

    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static boolean hash(String value) { return value != null && value.matches("[0-9a-f]{64}"); }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }
}
