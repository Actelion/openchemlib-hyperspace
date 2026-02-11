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
    private final double minSimilarity;
    private final long randomSeed;
    private final boolean enforceConnectorEquivalence;
    private final Map<String, Object> attributes;

    private SynthonDownsamplingRequest(Builder builder) {
        this.maxCenters = builder.maxCenters;
        this.minSimilarity = builder.minSimilarity;
        this.randomSeed = builder.randomSeed;
        this.enforceConnectorEquivalence = builder.enforceConnectorEquivalence;
        this.attributes = Collections.unmodifiableMap(new HashMap<>(builder.attributes));
    }

    public int getMaxCenters() {
        return maxCenters;
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
        private double minSimilarity = 0.0;
        private long randomSeed = 1L;
        private boolean enforceConnectorEquivalence = true;
        private final Map<String, Object> attributes = new HashMap<>();

        public Builder withMaxCenters(int maxCenters) {
            this.maxCenters = maxCenters;
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
}
