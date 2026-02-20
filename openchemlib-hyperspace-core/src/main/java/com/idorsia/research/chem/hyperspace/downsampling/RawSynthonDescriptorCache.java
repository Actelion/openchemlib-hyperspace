package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;

import java.io.Serializable;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;

/**
 * Thread-safe descriptor cache keyed by raw synthon identifiers.
 */
public final class RawSynthonDescriptorCache<T> implements Serializable {

    private final Map<String, T> cache = new ConcurrentHashMap<>();
    private final Function<RawSynthon, T> computer;

    public RawSynthonDescriptorCache(Function<RawSynthon, T> computer) {
        this.computer = computer;
    }

    public T getOrCompute(RawSynthon synthon) {
        String key = synthon.getReactionId() + ":" + synthon.getFragmentId();
        return cache.computeIfAbsent(key, k -> computer.apply(synthon));
    }

    public void clear() {
        cache.clear();
    }
}
