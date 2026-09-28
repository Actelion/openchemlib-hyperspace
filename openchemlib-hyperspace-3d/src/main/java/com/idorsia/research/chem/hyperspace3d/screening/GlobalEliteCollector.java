package com.idorsia.research.chem.hyperspace3d.screening;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

public final class GlobalEliteCollector {
    private final int capacity;
    private final PriorityQueue<ScoredProductFingerprint> heap =
            new PriorityQueue<>(Comparator.comparingDouble(ScoredProductFingerprint::objective));

    public GlobalEliteCollector(int capacity) {
        if (capacity < 1) throw new IllegalArgumentException("capacity must be positive");
        this.capacity = capacity;
    }

    public void offer(ScoredProductFingerprint hit) {
        if (heap.size() < capacity) heap.add(hit);
        else if (hit.objective() > heap.element().objective()) {
            heap.remove();
            heap.add(hit);
        }
    }

    public List<ScoredProductFingerprint> results() {
        List<ScoredProductFingerprint> result = new ArrayList<>(heap);
        result.sort(Comparator.comparingDouble(ScoredProductFingerprint::objective).reversed());
        return result;
    }
}
