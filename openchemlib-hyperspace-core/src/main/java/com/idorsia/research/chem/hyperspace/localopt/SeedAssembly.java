package com.idorsia.research.chem.hyperspace.localopt;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SeedAssembly implements Serializable {
    private final String reactionId;
    private final List<String> fragmentIds;
    private final double initialScore;

    public SeedAssembly(String reactionId, List<String> fragmentIds, double initialScore) {
        this.reactionId = reactionId;
        this.fragmentIds = Collections.unmodifiableList(new ArrayList<>(fragmentIds));
        this.initialScore = initialScore;
    }

    public String getReactionId() {
        return reactionId;
    }

    public List<String> getFragmentIds() {
        return fragmentIds;
    }

    public double getInitialScore() {
        return initialScore;
    }
}
