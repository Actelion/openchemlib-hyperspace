package com.idorsia.research.chem.hyperspace.seed;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class QueryAwareSeedFinderResult implements Serializable {

    private final List<ReactionStats> reactionStats;

    public QueryAwareSeedFinderResult(List<ReactionStats> reactionStats) {
        this.reactionStats = Collections.unmodifiableList(new ArrayList<>(reactionStats));
    }

    public List<ReactionStats> getReactionStats() {
        return reactionStats;
    }

    public long getTotalAttempts() {
        return reactionStats.stream().mapToLong(ReactionStats::getAttemptCount).sum();
    }

    public long getTotalAssemblies() {
        return reactionStats.stream().mapToLong(ReactionStats::getAssembledCount).sum();
    }

    public long getTotalThresholdPasses() {
        return reactionStats.stream().mapToLong(ReactionStats::getThresholdPassCount).sum();
    }

    public long getTotalPhesaEvaluations() {
        return reactionStats.stream().mapToLong(ReactionStats::getPhesaEvaluationCount).sum();
    }

    public long getTotalHits() {
        return reactionStats.stream().mapToLong(ReactionStats::getHitCount).sum();
    }

    public static final class ReactionStats implements Serializable {
        private final String reactionId;
        private long attemptCount;
        private long assembledCount;
        private long thresholdPassCount;
        private long phesaEvaluationCount;
        private long hitCount;

        public ReactionStats(String reactionId) {
            this.reactionId = reactionId;
        }

        void incrementAttempt() {
            attemptCount++;
        }

        void incrementAssemblies() {
            assembledCount++;
        }

        void incrementThresholdPass() {
            thresholdPassCount++;
        }

        void incrementPhesaEvaluation() {
            phesaEvaluationCount++;
        }

        void incrementHit() {
            hitCount++;
        }

        public String getReactionId() {
            return reactionId;
        }

        public long getAttemptCount() {
            return attemptCount;
        }

        public long getAssembledCount() {
            return assembledCount;
        }

        public long getThresholdPassCount() {
            return thresholdPassCount;
        }

        public long getPhesaEvaluationCount() {
            return phesaEvaluationCount;
        }

        public long getHitCount() {
            return hitCount;
        }
    }
}
