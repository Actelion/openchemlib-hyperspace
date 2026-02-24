package com.idorsia.research.chem.hyperspace.screening;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Intermediate representation that flows between sampling, micro optimization, and full optimization.
 */
public final class ScreeningCandidate {
    private final String reactionId;
    private final List<String> fragmentIds;
    private final String assembledIdcode;
    private final double similarity;
    private final int atomCount;
    private final int rotatableBonds;

    public ScreeningCandidate(String reactionId,
                              List<String> fragmentIds,
                              String assembledIdcode,
                              double similarity,
                              int atomCount,
                              int rotatableBonds) {
        this.reactionId = reactionId;
        this.fragmentIds = Collections.unmodifiableList(new ArrayList<>(fragmentIds));
        this.assembledIdcode = assembledIdcode;
        this.similarity = similarity;
        this.atomCount = atomCount;
        this.rotatableBonds = rotatableBonds;
    }

    public String getReactionId() {
        return reactionId;
    }

    public List<String> getFragmentIds() {
        return fragmentIds;
    }

    public String getAssembledIdcode() {
        return assembledIdcode;
    }

    public double getSimilarity() {
        return similarity;
    }

    public int getAtomCount() {
        return atomCount;
    }

    public int getRotatableBonds() {
        return rotatableBonds;
    }
}
