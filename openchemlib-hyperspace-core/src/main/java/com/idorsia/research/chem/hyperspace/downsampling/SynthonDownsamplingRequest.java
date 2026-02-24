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
    private final Map<String, Object> attributes;

    private SynthonDownsamplingRequest(Builder builder) {
        this.maxCenters = builder.maxCenters;
        this.sizeCapScale = builder.sizeCapScale;
        this.sizeCapOffset = builder.sizeCapOffset;
        this.minSimilarity = builder.minSimilarity;
        this.randomSeed = builder.randomSeed;
        this.enforceConnectorEquivalence = builder.enforceConnectorEquivalence;
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

    public Map<String, Object> getAttributes() {
        return attributes;
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
