package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

public class LoadedSynthonSpace {

    private SynthonSpace space;
    private String name;
    private int threads;

    public LoadedSynthonSpace(SynthonSpace space, String name) {
        this(space, name, Runtime.getRuntime().availableProcessors());
    }

    public LoadedSynthonSpace(SynthonSpace space, String name, int threads) {
        if (threads < 1) throw new IllegalArgumentException("threads must be positive");
        this.space = java.util.Objects.requireNonNull(space);
        this.name = name;
        this.threads = threads;
    }

    public int getThreads() { return threads; }

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
