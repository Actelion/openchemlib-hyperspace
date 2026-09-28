package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;
import java.util.Map;

public final class ProductFingerprintIndexManifest {
    public String artifactType = "hyperspace3d-product-fingerprint-index";
    public int formatVersion = 1;
    public String rawspaceIdentity;
    public String rawspaceSha256;
    public String downsampledRawspaceIdentity;
    public String downsampledRawspaceSha256;
    public String modelBundleHash;
    public String samplingAlgorithm;
    public long samplingSeed;
    public Map<String, Object> reactionWeightConfiguration;
    public Map<String, Object> molecularFilters;
    public String fingerprintDtype = "float32";
    public int embeddingDimension = 128;
    public long recordCount;
    public List<ProductFingerprintIndexShard> shards;
    public String buildConfigurationHash;
    public Map<String, Object> buildStatistics;

    public void validate() {
        if (!"hyperspace3d-product-fingerprint-index".equals(artifactType) || formatVersion != 1
                || embeddingDimension != 128 || !"float32".equals(fingerprintDtype)
                || rawspaceIdentity == null || rawspaceSha256 == null
                || downsampledRawspaceIdentity == null || downsampledRawspaceSha256 == null
                || modelBundleHash == null || samplingAlgorithm == null
                || reactionWeightConfiguration == null || molecularFilters == null
                || shards == null || recordCount < 0) {
            throw new IllegalArgumentException("incomplete or incompatible product-index manifest");
        }
        if (shards.stream().anyMatch(shard -> shard == null || shard.path() == null
                || shard.path().isBlank() || shard.recordCount() <= 0
                || shard.sha256() == null || !shard.sha256().matches("[0-9a-f]{64}"))) {
            throw new IllegalArgumentException("invalid product-index shard metadata");
        }
        long shardRecords = shards.stream().mapToLong(ProductFingerprintIndexShard::recordCount).sum();
        if (shardRecords != recordCount) throw new IllegalArgumentException("shard record count mismatch");
    }
}
