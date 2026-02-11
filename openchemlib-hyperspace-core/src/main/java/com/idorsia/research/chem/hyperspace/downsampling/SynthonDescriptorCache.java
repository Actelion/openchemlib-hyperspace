package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

/**
 * Thread-safe cache that stores descriptors by FragId idcode.
 */
public class SynthonDescriptorCache<T> implements Serializable {

    private final Map<String, T> cache = new ConcurrentHashMap<>();
    private final Function<SynthonSpace.FragId, T> computer;

    public SynthonDescriptorCache(Function<SynthonSpace.FragId, T> computer) {
        this.computer = Objects.requireNonNull(computer, "computer");
    }

    public T getOrCompute(SynthonSpace.FragId fragId) {
        return cache.computeIfAbsent(fragId.idcode, key -> computer.apply(fragId));
    }

    public void clear() {
        cache.clear();
    }

    public int size() {
        return cache.size();
    }
}
