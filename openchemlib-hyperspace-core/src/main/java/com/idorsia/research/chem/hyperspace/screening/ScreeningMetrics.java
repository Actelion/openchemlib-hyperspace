package com.idorsia.research.chem.hyperspace.screening;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

public final class ScreeningMetrics {

    private final AtomicLong sampled = new AtomicLong();
    private final AtomicLong duplicates = new AtomicLong();
    private final AtomicLong microOptimized = new AtomicLong();
    private final AtomicLong submitted = new AtomicLong();
    private final AtomicLong hits = new AtomicLong();
    private final Map<String, AtomicLong> reactionHits = new ConcurrentHashMap<>();

    public void incrementSampled() {
        sampled.incrementAndGet();
    }

    public void incrementDuplicates() {
        duplicates.incrementAndGet();
    }

    public void incrementMicroOptimized() {
        microOptimized.incrementAndGet();
    }

    public void incrementSubmitted() {
        submitted.incrementAndGet();
    }

    public void recordHit(String reactionId) {
        hits.incrementAndGet();
        reactionHits.computeIfAbsent(reactionId, key -> new AtomicLong()).incrementAndGet();
    }

    public long getSampled() {
        return sampled.get();
    }

    public long getDuplicates() {
        return duplicates.get();
    }

    public long getMicroOptimized() {
        return microOptimized.get();
    }

    public long getSubmitted() {
        return submitted.get();
    }

    public long getHits() {
        return hits.get();
    }

    public Map<String, AtomicLong> getReactionHits() {
        return reactionHits;
    }
}
