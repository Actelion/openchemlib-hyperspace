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

    public LocalOptimizationResult(String reactionId,
                                   List<BeamEntry> beamEntries,
                                   List<String> seedFragments) {
        this.reactionId = reactionId;
        this.beamEntries = Collections.unmodifiableList(new ArrayList<>(beamEntries));
        this.seedFragments = Collections.unmodifiableList(new ArrayList<>(seedFragments));
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
