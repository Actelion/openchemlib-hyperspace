package com.idorsia.research.chem.hyperspace3d.index.build;

import java.util.LinkedHashMap;
import java.util.Map;

/** Timings and counters for diagnosing the source-preparation/GPU balance. */
public final class MoleculeFingerprintIndexBuildMetrics {
    public long sourceRows;
    public long accepted;
    public long rejected;
    public long parsingAndFeaturizationNanos;
    public long tensorPackingNanos;
    public long encodingNanos;
    public long writingNanos;
    public long encoderWaitNanos;
    public Map<String, Long> rejections = new LinkedHashMap<>();
    public transient long lastProgressNanos;

    public void rejected(String reason) {
        rejected++;
        rejections.merge(reason, 1L, Long::sum);
    }

    public Map<String, Object> manifestView(long elapsedNanos) {
        Map<String, Object> result = new LinkedHashMap<>();
        double elapsed = elapsedNanos / 1_000_000_000d;
        result.put("elapsedSeconds", elapsed);
        result.put("sourceRowsPerSecond", sourceRows / Math.max(0.001, elapsed));
        result.put("acceptedPerSecond", accepted / Math.max(0.001, elapsed));
        result.put("sourceRows", sourceRows);
        result.put("accepted", accepted);
        result.put("rejected", rejected);
        result.put("rejections", new LinkedHashMap<>(rejections));
        result.put("parsingAndFeaturizationSeconds",
                parsingAndFeaturizationNanos / 1_000_000_000d);
        result.put("tensorPackingSeconds", tensorPackingNanos / 1_000_000_000d);
        result.put("encodingSeconds", encodingNanos / 1_000_000_000d);
        result.put("writingSeconds", writingNanos / 1_000_000_000d);
        result.put("encoderWaitForPreparedBatchSeconds",
                encoderWaitNanos / 1_000_000_000d);
        return result;
    }
}
