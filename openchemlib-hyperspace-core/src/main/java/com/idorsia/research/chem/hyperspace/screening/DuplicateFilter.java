package com.idorsia.research.chem.hyperspace.screening;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

/**
 * Tracks assembled IDCodes to avoid redundant work.
 */
public final class DuplicateFilter {

    private final ConcurrentMap<String, Boolean> cache = new ConcurrentHashMap<>();
    private final int maxEntries;

    public DuplicateFilter(int maxEntries) {
        this.maxEntries = Math.max(1, maxEntries);
    }

    public boolean markIfDuplicate(String idcode) {
        if (idcode == null || idcode.isBlank()) {
            return false;
        }
        boolean duplicate = cache.putIfAbsent(idcode, Boolean.TRUE) != null;
        evictIfNeeded();
        return duplicate;
    }

    private void evictIfNeeded() {
        if (cache.size() <= maxEntries) {
            return;
        }
        // simple coarse eviction: clear half the cache
        int target = maxEntries / 2;
        int removed = 0;
        for (String key : cache.keySet()) {
            if (cache.remove(key) != null) {
                removed++;
            }
            if (removed >= (cache.size() - target)) {
                break;
            }
        }
    }
}
