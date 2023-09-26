package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

public class LoadedSynthonSpace {

    private SynthonSpace space;
    private String name;

    public LoadedSynthonSpace(SynthonSpace space, String name) {
        this.space = space;
        this.name = name;
    }

    public SynthonSpace getSpace() {
        return space;
    }

    public void setSpace(SynthonSpace space) {
        this.space = space;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
