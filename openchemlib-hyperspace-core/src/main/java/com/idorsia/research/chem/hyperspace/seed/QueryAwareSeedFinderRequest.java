package com.idorsia.research.chem.hyperspace.seed;

import java.io.Serializable;
import java.nio.file.Path;

/**
 * Encapsulates all configuration parameters for the query-aware seed finder.
 */
public class QueryAwareSeedFinderRequest implements Serializable {

    private final long attemptsPerReaction;
    private final int minAtoms;
    private final int maxAtoms;
    private final int maxRotatableBonds;
    private final double minPhesaSimilarity;
    private final long randomSeed;
    private final int maxHitsPerFile;
    private final int numThreads;
    private final Path outputPath;

    private QueryAwareSeedFinderRequest(Builder builder) {
        this.attemptsPerReaction = builder.attemptsPerReaction;
        this.minAtoms = builder.minAtoms;
        this.maxAtoms = builder.maxAtoms;
        this.maxRotatableBonds = builder.maxRotatableBonds;
        this.minPhesaSimilarity = builder.minPhesaSimilarity;
        this.randomSeed = builder.randomSeed;
        this.maxHitsPerFile = builder.maxHitsPerFile;
        this.numThreads = builder.numThreads;
        this.outputPath = builder.outputPath;
    }

    public long getAttemptsPerReaction() {
        return attemptsPerReaction;
    }

    public int getMinAtoms() {
        return minAtoms;
    }

    public int getMaxAtoms() {
        return maxAtoms;
    }

    public int getMaxRotatableBonds() {
        return maxRotatableBonds;
    }

    public double getMinPhesaSimilarity() {
        return minPhesaSimilarity;
    }

    public long getRandomSeed() {
        return randomSeed;
    }

    public int getMaxHitsPerFile() {
        return maxHitsPerFile;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public Path getOutputPath() {
        return outputPath;
    }

    public static Builder builder(Path outputPath) {
        return new Builder(outputPath);
    }

    public static final class Builder {
        private long attemptsPerReaction = 1000;
        private int minAtoms = 1;
        private int maxAtoms = 0;
        private int maxRotatableBonds = Integer.MAX_VALUE;
        private double minPhesaSimilarity = 0.0;
        private long randomSeed = 13L;
        private int maxHitsPerFile = 0;
        private int numThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        private final Path outputPath;

        private Builder(Path outputPath) {
            this.outputPath = outputPath;
        }

        public Builder attemptsPerReaction(long attemptsPerReaction) {
            this.attemptsPerReaction = attemptsPerReaction;
            return this;
        }

        public Builder minAtoms(int minAtoms) {
            this.minAtoms = minAtoms;
            return this;
        }

        public Builder maxAtoms(int maxAtoms) {
            this.maxAtoms = maxAtoms;
            return this;
        }

        public Builder maxRotatableBonds(int maxRotatableBonds) {
            this.maxRotatableBonds = maxRotatableBonds;
            return this;
        }

        public Builder minPhesaSimilarity(double minPhesaSimilarity) {
            this.minPhesaSimilarity = minPhesaSimilarity;
            return this;
        }

        public Builder randomSeed(long randomSeed) {
            this.randomSeed = randomSeed;
            return this;
        }

        public Builder maxHitsPerFile(int maxHitsPerFile) {
            this.maxHitsPerFile = maxHitsPerFile;
            return this;
        }

        public Builder numThreads(int numThreads) {
            this.numThreads = numThreads;
            return this;
        }

        public QueryAwareSeedFinderRequest build() {
            return new QueryAwareSeedFinderRequest(this);
        }
    }
}
