package com.idorsia.research.chem.hyperspace.screening;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.DoubleAccumulator;
import java.util.concurrent.atomic.DoubleAdder;
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
    private final ScoreStats sampledScoreStats = new ScoreStats();
    private final ScoreStats preMicroScoreStats = new ScoreStats();
    private final ScoreStats postMicroScoreStats = new ScoreStats();
    private final ScoreStats preFullScoreStats = new ScoreStats();
    private final ScoreStats postFullAllScoreStats = new ScoreStats();
    private final ScoreStats postFullReportedScoreStats = new ScoreStats();
    private final ScoreStats microGainStats = new ScoreStats();
    private final ScoreStats fullGainStats = new ScoreStats();

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

    public void recordSampledScore(double score) {
        sampledScoreStats.add(score);
    }

    public void recordPreMicroScore(double score) {
        preMicroScoreStats.add(score);
    }

    public void recordPostMicroScore(double score) {
        postMicroScoreStats.add(score);
    }

    public void recordPreFullScore(double score) {
        preFullScoreStats.add(score);
    }

    public void recordPostFullAllScore(double score) {
        postFullAllScoreStats.add(score);
    }

    public void recordPostFullReportedScore(double score) {
        postFullReportedScoreStats.add(score);
    }

    public void recordMicroGain(double gain) {
        microGainStats.add(gain);
    }

    public void recordFullGain(double gain) {
        fullGainStats.add(gain);
    }

    public ScoreStatsSnapshot getSampledScoreStats() {
        return sampledScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getPreMicroScoreStats() {
        return preMicroScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getPostMicroScoreStats() {
        return postMicroScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getPreFullScoreStats() {
        return preFullScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getPostFullAllScoreStats() {
        return postFullAllScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getPostFullReportedScoreStats() {
        return postFullReportedScoreStats.snapshot();
    }

    public ScoreStatsSnapshot getMicroGainStats() {
        return microGainStats.snapshot();
    }

    public ScoreStatsSnapshot getFullGainStats() {
        return fullGainStats.snapshot();
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

    private static final class ScoreStats {
        private final LongAdder count = new LongAdder();
        private final DoubleAdder sum = new DoubleAdder();
        private final DoubleAdder sumSquares = new DoubleAdder();
        private final DoubleAccumulator min = new DoubleAccumulator(Double::min, Double.POSITIVE_INFINITY);
        private final DoubleAccumulator max = new DoubleAccumulator(Double::max, Double.NEGATIVE_INFINITY);

        private void add(double value) {
            if (!Double.isFinite(value)) {
                return;
            }
            count.increment();
            sum.add(value);
            sumSquares.add(value * value);
            min.accumulate(value);
            max.accumulate(value);
        }

        private ScoreStatsSnapshot snapshot() {
            long n = count.sum();
            if (n == 0L) {
                return ScoreStatsSnapshot.empty();
            }
            double sumValue = sum.sum();
            double mean = sumValue / n;
            double variance = Math.max(0.0, (sumSquares.sum() / n) - (mean * mean));
            return new ScoreStatsSnapshot(n, mean, Math.sqrt(variance), min.get(), max.get());
        }
    }

    public static final class ScoreStatsSnapshot {
        private static final ScoreStatsSnapshot EMPTY =
                new ScoreStatsSnapshot(0L, Double.NaN, Double.NaN, Double.NaN, Double.NaN);

        private final long count;
        private final double mean;
        private final double stdev;
        private final double min;
        private final double max;

        private ScoreStatsSnapshot(long count, double mean, double stdev, double min, double max) {
            this.count = count;
            this.mean = mean;
            this.stdev = stdev;
            this.min = min;
            this.max = max;
        }

        private static ScoreStatsSnapshot empty() {
            return EMPTY;
        }

        public long getCount() {
            return count;
        }

        public double getMean() {
            return mean;
        }

        public double getStdev() {
            return stdev;
        }

        public double getMin() {
            return min;
        }

        public double getMax() {
            return max;
        }
    }
}
