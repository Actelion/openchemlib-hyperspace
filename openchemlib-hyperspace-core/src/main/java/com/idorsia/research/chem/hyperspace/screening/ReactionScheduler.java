package com.idorsia.research.chem.hyperspace.screening;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Simple weighted scheduler that picks reactions proportional to their synthon counts.
 */
public final class ReactionScheduler {

    private final List<Entry> entries;
    private final double totalWeight;

    public ReactionScheduler(Map<String, Integer> reactionSizes) {
        if (reactionSizes.isEmpty()) {
            throw new IllegalArgumentException("No reactions available for scheduling");
        }
        List<Entry> list = new ArrayList<>();
        double sum = 0d;
        for (Map.Entry<String, Integer> entry : reactionSizes.entrySet()) {
            double weight = Math.max(1d, Math.sqrt(Math.max(1, entry.getValue())));
            sum += weight;
            list.add(new Entry(entry.getKey(), weight));
        }
        this.entries = Collections.unmodifiableList(list);
        this.totalWeight = sum;
    }

    public String pick(Random random) {
        double target = random.nextDouble() * totalWeight;
        double cumulative = 0d;
        for (Entry entry : entries) {
            cumulative += entry.weight;
            if (target <= cumulative) {
                return entry.reactionId;
            }
        }
        return entries.get(entries.size() - 1).reactionId;
    }

    private static final class Entry {
        private final String reactionId;
        private final double weight;

        private Entry(String reactionId, double weight) {
            this.reactionId = reactionId;
            this.weight = weight;
        }
    }
}
