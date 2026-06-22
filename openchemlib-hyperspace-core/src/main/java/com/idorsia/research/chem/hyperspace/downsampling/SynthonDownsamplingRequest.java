package com.idorsia.research.chem.hyperspace.downsampling;

import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Immutable descriptor for how a SynthonDownsampler should behave.
 */
public class SynthonDownsamplingRequest implements Serializable {

    private final int maxCenters;
    private final double sizeCapScale;
    private final double sizeCapOffset;
    private final double minSimilarity;
    private final long randomSeed;
    private final boolean enforceConnectorEquivalence;
    private final boolean includeClusterMembers;
    private final transient ProgressReporter progressReporter;
    private final int progressReportIntervalSeconds;
    private final Map<String, Object> attributes;

    private SynthonDownsamplingRequest(Builder builder) {
        this.maxCenters = builder.maxCenters;
        this.sizeCapScale = builder.sizeCapScale;
        this.sizeCapOffset = builder.sizeCapOffset;
        this.minSimilarity = builder.minSimilarity;
        this.randomSeed = builder.randomSeed;
        this.enforceConnectorEquivalence = builder.enforceConnectorEquivalence;
        this.includeClusterMembers = builder.includeClusterMembers;
        this.progressReporter = builder.progressReporter;
        this.progressReportIntervalSeconds = builder.progressReportIntervalSeconds;
        this.attributes = Collections.unmodifiableMap(new HashMap<>(builder.attributes));
    }

    public int getMaxCenters() {
        return maxCenters;
    }

    public double getSizeCapScale() {
        return sizeCapScale;
    }

    public double getSizeCapOffset() {
        return sizeCapOffset;
    }

    public int getEffectiveMaxCenters(int totalSynthons) {
        int effective = maxCenters;
        int sizeBased = computeSizeBasedCap(totalSynthons);
        if (sizeBased > 0 && (effective <= 0 || sizeBased < effective)) {
            effective = sizeBased;
        }
        return effective;
    }

    public double getMinSimilarity() {
        return minSimilarity;
    }

    public long getRandomSeed() {
        return randomSeed;
    }

    public boolean isEnforceConnectorEquivalence() {
        return enforceConnectorEquivalence;
    }

    public boolean isIncludeClusterMembers() {
        return includeClusterMembers;
    }

    public Map<String, Object> getAttributes() {
        return attributes;
    }

    public ProgressReporter getProgressReporter() {
        return progressReporter == null ? ProgressReporter.NO_OP : progressReporter;
    }

    public int getProgressReportIntervalSeconds() {
        return progressReportIntervalSeconds;
    }

    public static Builder builder() {
        return new Builder();
    }

    public static class Builder {
        private int maxCenters = 0;
        private double sizeCapScale = 0.0;
        private double sizeCapOffset = 0.0;
        private double minSimilarity = 0.0;
        private long randomSeed = 1L;
        private boolean enforceConnectorEquivalence = true;
        private boolean includeClusterMembers = false;
        private ProgressReporter progressReporter = ProgressReporter.NO_OP;
        private int progressReportIntervalSeconds = 60;
        private final Map<String, Object> attributes = new HashMap<>();

        public Builder withMaxCenters(int maxCenters) {
            this.maxCenters = maxCenters;
            return this;
        }

        public Builder withSizeCapScale(double scale) {
            this.sizeCapScale = scale;
            return this;
        }

        public Builder withSizeCapOffset(double offset) {
            this.sizeCapOffset = offset;
            return this;
        }

        public Builder withMinSimilarity(double minSimilarity) {
            this.minSimilarity = minSimilarity;
            return this;
        }

        public Builder withRandomSeed(long randomSeed) {
            this.randomSeed = randomSeed;
            return this;
        }

        public Builder enforceConnectorEquivalence(boolean enforce) {
            this.enforceConnectorEquivalence = enforce;
            return this;
        }

        public Builder includeClusterMembers(boolean includeClusterMembers) {
            this.includeClusterMembers = includeClusterMembers;
            return this;
        }

        public Builder withProgressReporter(ProgressReporter progressReporter) {
            this.progressReporter = progressReporter == null ? ProgressReporter.NO_OP : progressReporter;
            return this;
        }

        public Builder withProgressReportIntervalSeconds(int seconds) {
            if (seconds < 0) {
                throw new IllegalArgumentException("Progress report interval must be >= 0");
            }
            this.progressReportIntervalSeconds = seconds;
            return this;
        }

        public Builder putAttribute(String key, Object value) {
            if (key != null && value != null) {
                this.attributes.put(key, value);
            }
            return this;
        }

        public SynthonDownsamplingRequest build() {
            return new SynthonDownsamplingRequest(this);
        }
    }

    public interface ProgressReporter {
        ProgressReporter NO_OP = event -> { };

        void onProgress(ProgressEvent event);
    }

    public static final class ProgressEvent {
        private final String reactionId;
        private final int fragmentIndex;
        private final int processedSynthons;
        private final int totalSynthons;
        private final int centerCount;
        private final int effectiveMaxCenters;
        private final long elapsedMillis;
        private final int skippedDescriptors;
        private final boolean completed;

        public ProgressEvent(String reactionId,
                             int fragmentIndex,
                             int processedSynthons,
                             int totalSynthons,
                             int centerCount,
                             int effectiveMaxCenters,
                             long elapsedMillis,
                             int skippedDescriptors,
                             boolean completed) {
            this.reactionId = reactionId;
            this.fragmentIndex = fragmentIndex;
            this.processedSynthons = processedSynthons;
            this.totalSynthons = totalSynthons;
            this.centerCount = centerCount;
            this.effectiveMaxCenters = effectiveMaxCenters;
            this.elapsedMillis = elapsedMillis;
            this.skippedDescriptors = skippedDescriptors;
            this.completed = completed;
        }

        public String getReactionId() {
            return reactionId;
        }

        public int getFragmentIndex() {
            return fragmentIndex;
        }

        public int getProcessedSynthons() {
            return processedSynthons;
        }

        public int getTotalSynthons() {
            return totalSynthons;
        }

        public int getCenterCount() {
            return centerCount;
        }

        public int getEffectiveMaxCenters() {
            return effectiveMaxCenters;
        }

        public long getElapsedMillis() {
            return elapsedMillis;
        }

        public int getSkippedDescriptors() {
            return skippedDescriptors;
        }

        public boolean isCompleted() {
            return completed;
        }
    }

    private int computeSizeBasedCap(int totalSynthons) {
        if (totalSynthons <= 0) {
            return 0;
        }
        if (sizeCapScale <= 0.0 && sizeCapOffset <= 0.0) {
            return 0;
        }
        double value = sizeCapScale * Math.sqrt(totalSynthons) + sizeCapOffset;
        int cap = (int) Math.ceil(value);
        cap = Math.max(1, cap);
        return Math.min(cap, totalSynthons);
    }
}
