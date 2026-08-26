package com.idorsia.research.chem.hyperspace3d.mining;

import java.util.List;
import java.util.Map;

/** Provenance contract for an immutable query-centric exact-PheSA cache. */
public final class PheSAMiningManifest {
    public String artifactType = "hyperspace-phesa-query-pair-cache";
    public int formatVersion = 1;
    public int embeddingDim = 128;
    public String descriptorDtype = "float16";
    public double phesaPpWeight = 0.5;
    public int maxConformers = 8;
    public int miningRound;
    public long samplingSeed = 17;
    public String fingerprintModelHash;
    public String sourceIndexHash;
    public String candidateUniverseHash;
    public String queryRegistryHash;
    public String oclVersion;
    public long moleculeCount;
    public long queryCount;
    public long pairCount;
    public List<Double> rankBoundaries;
    public List<Integer> rankQuotas;
    public int randomQuota;
    public Map<String, Integer> chemistryStratumQuotas;
    public Map<String, Long> selectionMask;
    public Map<String, FileEntry> files;

    public void validate() {
        require("hyperspace-phesa-query-pair-cache".equals(artifactType), "unsupported artifactType");
        require(formatVersion == 1 && embeddingDim == 128
                && "float16".equals(descriptorDtype), "unsupported cache representation");
        require(maxConformers == 8 && Double.isFinite(phesaPpWeight)
                && phesaPpWeight >= 0 && phesaPpWeight <= 1, "invalid PheSA contract");
        require(miningRound >= 0 && moleculeCount >= 0 && queryCount >= 0 && pairCount >= 0,
                "invalid cache counts");
        require(hash(fingerprintModelHash) && hash(sourceIndexHash) && hash(candidateUniverseHash)
                && hash(queryRegistryHash),
                "cache provenance hashes are incomplete");
        require(oclVersion != null && !oclVersion.isBlank(), "OCL version is required");
        require(rankBoundaries != null && rankBoundaries.size() == 7
                && rankQuotas != null && rankQuotas.size() == 6 && randomQuota >= 0
                && chemistryStratumQuotas != null,
                "rank selection contract is incomplete");
        require(selectionMask != null && selectionMask.equals(MiningSelection.manifestNames()),
                "selection-mask contract differs from format version 1");
        require(files != null && files.keySet().containsAll(
                List.of("molecules", "descriptors", "query_index")), "file manifest is incomplete");
        for (FileEntry entry : files.values()) {
            require(entry != null && entry.path != null && !entry.path.isBlank()
                    && hash(entry.sha256) && entry.sizeBytes >= 0, "invalid file entry");
        }
    }

    private static boolean hash(String value) {
        return value != null && value.matches("[0-9a-f]{64}");
    }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    public static final class FileEntry {
        public String path;
        public String sha256;
        public long sizeBytes;
        public FileEntry() {}
        public FileEntry(String path, String sha256, long sizeBytes) {
            this.path = path; this.sha256 = sha256; this.sizeBytes = sizeBytes;
        }
    }
}
