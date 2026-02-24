package com.idorsia.research.chem.hyperspace.screening;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.concurrent.atomic.AtomicLong;

public final class ScreeningMetrics {

    private final AtomicLong sampled = new AtomicLong();
    private final AtomicLong duplicateSeeds = new AtomicLong();
    private final AtomicLong microOptimized = new AtomicLong();
    private final AtomicLong submitted = new AtomicLong();
    private final AtomicLong hits = new AtomicLong();
    private final Map<String, AtomicLong> reactionHits = new ConcurrentHashMap<>();
    private final LongAdder candidateComparisons = new LongAdder();
    private final LongAdder microOptComparisons = new LongAdder();
    private final LongAdder fullOptComparisons = new LongAdder();
    private final Map<String, ReactionSamplingStats> reactionSampling = new ConcurrentHashMap<>();

    public void incrementSampled() {
        sampled.incrementAndGet();
    }

    public void incrementDuplicateSeeds() {
        duplicateSeeds.incrementAndGet();
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

    void recordSamplingAttempt(String reactionId) {
        reactionSampling.computeIfAbsent(reactionId, key -> new ReactionSamplingStats())
                .attempts.increment();
    }

    void recordSamplingSuccess(String reactionId) {
        reactionSampling.computeIfAbsent(reactionId, key -> new ReactionSamplingStats())
                .successes.increment();
    }

    void incrementCandidateComparisons() {
        candidateComparisons.increment();
    }

    void incrementMicroOptComparisons() {
        microOptComparisons.increment();
    }

    void incrementFullOptComparisons() {
        fullOptComparisons.increment();
    }

    LongAdder candidateComparisonCounter() {
        return candidateComparisons;
    }

    LongAdder microOptComparisonCounter() {
        return microOptComparisons;
    }

    LongAdder fullOptComparisonCounter() {
        return fullOptComparisons;
    }

    public long getSampled() {
        return sampled.get();
    }

    public long getDuplicateSeeds() {
        return duplicateSeeds.get();
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

    public Map<String, ReactionSamplingStats> getReactionSampling() {
        return reactionSampling;
    }

    public long getCandidateComparisons() {
        return candidateComparisons.sum();
    }

    public long getMicroOptComparisons() {
        return microOptComparisons.sum();
    }

    public long getFullOptComparisons() {
        return fullOptComparisons.sum();
    }

    public long getTotalComparisons() {
        return getCandidateComparisons() + getMicroOptComparisons() + getFullOptComparisons();
    }

    public static final class ReactionSamplingStats {
        private final LongAdder attempts = new LongAdder();
        private final LongAdder successes = new LongAdder();

        public long getAttempts() {
            return attempts.sum();
        }

        public long getSuccesses() {
            return successes.sum();
        }
    }
}
