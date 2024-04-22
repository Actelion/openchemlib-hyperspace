package com.idorsia.research.chem.hyperspace.service;

public class SubstructureSearchTaskConfiguration {

    int maxSplits = 3;
    int numThreads = -1;
    int maxCombinatorialHits = 64;

    public SubstructureSearchTaskConfiguration() {}

    public SubstructureSearchTaskConfiguration(int maxSplits, int numThreads, int maxCombinatorialHits) {
        this.maxSplits = maxSplits;
        this.numThreads = numThreads;
        this.maxCombinatorialHits = maxCombinatorialHits;
    }

    public int getMaxSplits() {
        return maxSplits;
    }

    public void setMaxSplits(int maxSplits) {
        this.maxSplits = maxSplits;
    }

    public int getNumThreads() {
        return numThreads;
    }

    public void setNumThreads(int numThreads) {
        this.numThreads = numThreads;
    }

    public int getMaxCombinatorialHits() {
        return maxCombinatorialHits;
    }

    public void setMaxCombinatorialHits(int maxCombinatorialHits) {
        this.maxCombinatorialHits = maxCombinatorialHits;
    }
}
