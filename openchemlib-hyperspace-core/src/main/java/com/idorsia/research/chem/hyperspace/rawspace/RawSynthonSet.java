package com.idorsia.research.chem.hyperspace.rawspace;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * Immutable view of a single synthon set inside a {@link RawSynthonSpace}.
 */
public final class RawSynthonSet implements Serializable {

    private final String reactionId;
    private final int fragmentIndex;
    private final List<RawSynthon> synthons;

    public RawSynthonSet(String reactionId, int fragmentIndex, List<RawSynthon> synthons) {
        this.reactionId = Objects.requireNonNull(reactionId, "reactionId");
        this.fragmentIndex = fragmentIndex;
        this.synthons = Collections.unmodifiableList(Objects.requireNonNull(synthons, "synthons"));
    }

    public String getReactionId() {
        return reactionId;
    }

    public int getFragmentIndex() {
        return fragmentIndex;
    }

    public List<RawSynthon> getSynthons() {
        return synthons;
    }

    public SynthonSpace.FragType toFragType() {
        return new SynthonSpace.FragType(reactionId, fragmentIndex);
    }
}
