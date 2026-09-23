package com.idorsia.research.chem.hyperspace.stats;

public final class RawSynthonSpaceStatsOptions {
    private final int examplesPerReaction;
    private final int productSamplesPerReaction;
    private final long seed;
    private final int threads;
    private final int maxReactions;
    private final ProgressListener progressListener;

    private RawSynthonSpaceStatsOptions(Builder builder) {
        this.examplesPerReaction = Math.max(0, builder.examplesPerReaction);
        this.productSamplesPerReaction = Math.max(0, builder.productSamplesPerReaction);
        this.seed = builder.seed;
        this.threads = Math.max(1, builder.threads);
        this.maxReactions = Math.max(0, builder.maxReactions);
        this.progressListener = builder.progressListener == null ? ProgressListener.NOOP : builder.progressListener;
    }

    public int getExamplesPerReaction() {
        return examplesPerReaction;
    }

    public int getProductSamplesPerReaction() {
        return productSamplesPerReaction;
    }

    public long getSeed() {
        return seed;
    }

    public int getThreads() {
        return threads;
    }

    public int getMaxReactions() {
        return maxReactions;
    }

    public ProgressListener getProgressListener() {
        return progressListener;
    }

    public static Builder builder() {
        return new Builder();
    }

    public interface ProgressListener {
        ProgressListener NOOP = (completed, total, reactionId) -> { };

        void onReactionCompleted(int completed, int total, String reactionId);
    }

    public static final class Builder {
        private int examplesPerReaction = 10;
        private int productSamplesPerReaction = 100;
        private long seed = 13L;
        private int threads = Runtime.getRuntime().availableProcessors();
        private int maxReactions = 0;
        private ProgressListener progressListener;

        private Builder() {
        }

        public Builder examplesPerReaction(int examplesPerReaction) {
            this.examplesPerReaction = examplesPerReaction;
            return this;
        }

        public Builder productSamplesPerReaction(int productSamplesPerReaction) {
            this.productSamplesPerReaction = productSamplesPerReaction;
            return this;
        }

        public Builder seed(long seed) {
            this.seed = seed;
            return this;
        }

        public Builder threads(int threads) {
            this.threads = threads;
            return this;
        }

        public Builder maxReactions(int maxReactions) {
            this.maxReactions = maxReactions;
            return this;
        }

        public Builder progressListener(ProgressListener progressListener) {
            this.progressListener = progressListener;
            return this;
        }

        public RawSynthonSpaceStatsOptions build() {
            return new RawSynthonSpaceStatsOptions(this);
        }
    }
}
