package com.idorsia.research.chem.hyperspace3d.workflow;

import java.util.*;

public final class FingerprintWorkflowManifest {
    public int formatVersion = 1;
    public String workflowId;
    public String configHash;
    public String inputSha256;
    public String modelBundleHash;
    public String compactBundleHash;
    public String environmentSetupSha256;
    public Map<String, String> applicationFiles = new TreeMap<>();
    public long sourceRowCount;
    public List<Partition> partitions = new ArrayList<>();

    public record Partition(int index, String input, long firstSourceRow, long sourceRowCount, String sha256) {}
    public record Completion(String workflowId, int partition, Map<String, String> files) {}
}
