package com.idorsia.research.chem.hyperspace3d.screening;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class ReactionEliteCollector {
    private final int capacity;
    private final Map<String, GlobalEliteCollector> collectors = new LinkedHashMap<>();

    public ReactionEliteCollector(int capacity) {
        if (capacity < 1) throw new IllegalArgumentException("capacity must be positive");
        this.capacity = capacity;
    }

    public void offer(ScoredProductFingerprint hit) {
        collectors.computeIfAbsent(hit.tuple().reactionId(), ignored ->
                new GlobalEliteCollector(capacity)).offer(hit);
    }

    public Map<String, List<ScoredProductFingerprint>> results() {
        Map<String, List<ScoredProductFingerprint>> result = new LinkedHashMap<>();
        collectors.forEach((reaction, collector) -> result.put(reaction, collector.results()));
        return result;
    }
}
