package com.idorsia.research.chem.hyperspace3d.index.build;

import java.util.LinkedHashMap;
import java.util.Map;

/** Mutable run metrics; snapshots are safe to serialize into manifests. */
public final class ProductFingerprintIndexBuildMetrics {
    public long proposals;
    public long accepted;
    public long assemblyNanos;
    public long featurizationNanos;
    public long tensorPackingNanos;
    public long encodingNanos;
    public long writingNanos;
    public long hashingNanos;
    public Map<String, Long> rejections = new LinkedHashMap<>();
    public Map<String, Long> acceptedByReaction = new LinkedHashMap<>();

    public void rejected(String reason) { rejections.merge(reason, 1L, Long::sum); }
    public void accepted(String reaction) {
        accepted++;
        acceptedByReaction.merge(reaction, 1L, Long::sum);
    }

    public Map<String, Object> manifestView(long elapsedNanos) {
        Map<String, Object> map = new LinkedHashMap<>();
        map.put("proposals", proposals);
        map.put("accepted", accepted);
        map.put("rejections", new LinkedHashMap<>(rejections));
        map.put("acceptedByReaction", new LinkedHashMap<>(acceptedByReaction));
        map.put("elapsedSeconds", elapsedNanos / 1_000_000_000d);
        map.put("assemblySeconds", assemblyNanos / 1_000_000_000d);
        map.put("featurizationSeconds", featurizationNanos / 1_000_000_000d);
        map.put("tensorPackingSeconds", tensorPackingNanos / 1_000_000_000d);
        map.put("encodingSeconds", encodingNanos / 1_000_000_000d);
        map.put("writingSeconds", writingNanos / 1_000_000_000d);
        map.put("hashingSeconds", hashingNanos / 1_000_000_000d);
        return map;
    }
}
