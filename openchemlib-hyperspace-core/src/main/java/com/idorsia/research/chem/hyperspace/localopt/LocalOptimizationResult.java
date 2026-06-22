package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class LocalOptimizationResult implements Serializable {

    private final String reactionId;
    private final List<BeamEntry> beamEntries;
    private final List<String> seedFragments;
    private final double bestObservedScore;
    private final OptimizationStats stats;

    public LocalOptimizationResult(String reactionId,
                                   List<BeamEntry> beamEntries,
                                   List<String> seedFragments) {
        this(reactionId, beamEntries, seedFragments, Double.NaN);
    }

    public LocalOptimizationResult(String reactionId,
                                   List<BeamEntry> beamEntries,
                                   List<String> seedFragments,
                                   double bestObservedScore) {
        this(reactionId, beamEntries, seedFragments, bestObservedScore, OptimizationStats.unavailable());
    }

    public LocalOptimizationResult(String reactionId,
                                   List<BeamEntry> beamEntries,
                                   List<String> seedFragments,
                                   double bestObservedScore,
                                   OptimizationStats stats) {
        this.reactionId = reactionId;
        this.beamEntries = Collections.unmodifiableList(new ArrayList<>(beamEntries));
        this.seedFragments = Collections.unmodifiableList(new ArrayList<>(seedFragments));
        this.bestObservedScore = bestObservedScore;
        this.stats = stats == null ? OptimizationStats.unavailable() : stats;
    }

    public String getReactionId() {
        return reactionId;
    }

    public List<BeamEntry> getBeamEntries() {
        return beamEntries;
    }

    public List<String> getSeedFragments() {
        return seedFragments;
    }

    /**
     * Best score seen during optimization, regardless of reporting filters.
     */
    public double getBestObservedScore() {
        return bestObservedScore;
    }

    public OptimizationStats getStats() {
        return stats;
    }

    public enum StopReason {
        UNAVAILABLE,
        SEED_SCORE_FAILED,
        EMPTY_BEAM,
        PATIENCE,
        MAX_ROUNDS
    }

    public static final class OptimizationStats implements Serializable {
        private static final OptimizationStats UNAVAILABLE = new OptimizationStats(
                Double.NaN,
                Double.NaN,
                0,
                0,
                0,
                StopReason.UNAVAILABLE);

        private final double initialScore;
        private final double bestScore;
        private final int roundsAttempted;
        private final int roundsCompleted;
        private final int scoredCandidates;
        private final StopReason stopReason;

        public OptimizationStats(double initialScore,
                                 double bestScore,
                                 int roundsAttempted,
                                 int roundsCompleted,
                                 int scoredCandidates,
                                 StopReason stopReason) {
            this.initialScore = initialScore;
            this.bestScore = bestScore;
            this.roundsAttempted = roundsAttempted;
            this.roundsCompleted = roundsCompleted;
            this.scoredCandidates = scoredCandidates;
            this.stopReason = stopReason == null ? StopReason.UNAVAILABLE : stopReason;
        }

        public static OptimizationStats unavailable() {
            return UNAVAILABLE;
        }

        public double getInitialScore() {
            return initialScore;
        }

        public double getBestScore() {
            return bestScore;
        }

        public int getRoundsAttempted() {
            return roundsAttempted;
        }

        public int getRoundsCompleted() {
            return roundsCompleted;
        }

        public int getScoredCandidates() {
            return scoredCandidates;
        }

        public StopReason getStopReason() {
            return stopReason;
        }
    }

    public static final class BeamEntry implements Serializable {
        private final List<SynthonSpace.FragId> fragments;
        private final double phesaSimilarity;
        private final int atomCount;
        private final int rotatableBonds;
        private final String idcode;
        private final int originatingRound;

        public BeamEntry(List<SynthonSpace.FragId> fragments,
                         double phesaSimilarity,
                         int atomCount,
                         int rotatableBonds,
                         String idcode,
                         int originatingRound) {
            this.fragments = Collections.unmodifiableList(new ArrayList<>(fragments));
            this.phesaSimilarity = phesaSimilarity;
            this.atomCount = atomCount;
            this.rotatableBonds = rotatableBonds;
            this.idcode = idcode;
            this.originatingRound = originatingRound;
        }

        public List<SynthonSpace.FragId> getFragments() {
            return fragments;
        }

        public double getPhesaSimilarity() {
            return phesaSimilarity;
        }

        /**
         * Generic accessor for the primary optimization score (PheSA similarity by default).
         */
        public double getScore() {
            return phesaSimilarity;
        }

        public int getAtomCount() {
            return atomCount;
        }

        public int getRotatableBonds() {
            return rotatableBonds;
        }

        public String getIdcode() {
            return idcode;
        }

        public int getOriginatingRound() {
            return originatingRound;
        }
    }
}
