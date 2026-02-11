package com.idorsia.research.chem.hyperspace.localopt;

import java.io.Serializable;

public class LocalOptimizationRequest implements Serializable {

    private final int beamSize;
    private final int neighborPoolSize;
    private final int sampledNeighbors;
    private final int perPositionCap;
    private final int maxRounds;
    private final int patience;
    private final double minPhesaSimilarity;
    private final double improvementTolerance;
    private final long randomSeed;
    private final LocalOptimizationLogLevel logLevel;

    private LocalOptimizationRequest(Builder builder) {
        this.beamSize = builder.beamSize;
        this.neighborPoolSize = builder.neighborPoolSize;
        this.sampledNeighbors = builder.sampledNeighbors;
        this.perPositionCap = builder.perPositionCap;
        this.maxRounds = builder.maxRounds;
        this.patience = builder.patience;
        this.minPhesaSimilarity = builder.minPhesaSimilarity;
        this.improvementTolerance = builder.improvementTolerance;
        this.randomSeed = builder.randomSeed;
        this.logLevel = builder.logLevel;
    }

    public int getBeamSize() {
        return beamSize;
    }

    public int getNeighborPoolSize() {
        return neighborPoolSize;
    }

    public int getSampledNeighbors() {
        return sampledNeighbors;
    }

    public int getPerPositionCap() {
        return perPositionCap;
    }

    public int getMaxRounds() {
        return maxRounds;
    }

    public int getPatience() {
        return patience;
    }

    public double getMinPhesaSimilarity() {
        return minPhesaSimilarity;
    }

    public double getImprovementTolerance() {
        return improvementTolerance;
    }

    public long getRandomSeed() {
        return randomSeed;
    }

    public LocalOptimizationLogLevel getLogLevel() {
        return logLevel;
    }

    public static Builder builder() {
        return new Builder();
    }

    public static final class Builder {
        private int beamSize = 10;
        private int neighborPoolSize = 100;
        private int sampledNeighbors = 10;
        private int perPositionCap = 2;
        private int maxRounds = 6;
        private int patience = 3;
        private double minPhesaSimilarity = 0.6;
        private double improvementTolerance = 1e-4;
        private long randomSeed = 13L;
        private LocalOptimizationLogLevel logLevel = LocalOptimizationLogLevel.NONE;

        public Builder beamSize(int beamSize) {
            this.beamSize = beamSize;
            return this;
        }

        public Builder neighborPoolSize(int neighborPoolSize) {
            this.neighborPoolSize = neighborPoolSize;
            return this;
        }

        public Builder sampledNeighbors(int sampledNeighbors) {
            this.sampledNeighbors = sampledNeighbors;
            return this;
        }

        public Builder perPositionCap(int perPositionCap) {
            this.perPositionCap = perPositionCap;
            return this;
        }

        public Builder maxRounds(int maxRounds) {
            this.maxRounds = maxRounds;
            return this;
        }

        public Builder patience(int patience) {
            this.patience = patience;
            return this;
        }

        public Builder minPhesaSimilarity(double minPhesaSimilarity) {
            this.minPhesaSimilarity = minPhesaSimilarity;
            return this;
        }

        public Builder improvementTolerance(double improvementTolerance) {
            this.improvementTolerance = improvementTolerance;
            return this;
        }

        public Builder randomSeed(long randomSeed) {
            this.randomSeed = randomSeed;
            return this;
        }

        public Builder logLevel(LocalOptimizationLogLevel logLevel) {
            this.logLevel = logLevel;
            return this;
        }

        public LocalOptimizationRequest build() {
            return new LocalOptimizationRequest(this);
        }
    }
}
