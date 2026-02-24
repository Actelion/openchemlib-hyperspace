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

    public ReactionScheduler(Map<String, List<Integer>> reactionSizes,
                             double minWeight,
                             double exponent) {
        if (reactionSizes.isEmpty()) {
            throw new IllegalArgumentException("No reactions available for scheduling");
        }
        if (exponent <= 0) {
            throw new IllegalArgumentException("Exponent must be positive");
        }
        List<Entry> list = new ArrayList<>();
        double maxLog = Double.NEGATIVE_INFINITY;
        List<ReactionWeight> weights = new ArrayList<>();
        for (Map.Entry<String, List<Integer>> entry : reactionSizes.entrySet()) {
            double logProduct = 0d;
            for (int size : entry.getValue()) {
                if (size > 0) {
                    logProduct += Math.log(size);
                }
            }
            double logWeight = exponent * logProduct;
            weights.add(new ReactionWeight(entry.getKey(), logWeight));
            maxLog = Math.max(maxLog, logWeight);
        }
        double sum = 0d;
        for (ReactionWeight weight : weights) {
            double scaled = Math.exp(weight.logWeight - maxLog);
            if (Double.isNaN(scaled) || Double.isInfinite(scaled)) {
                scaled = 0d;
            }
            double finalWeight = Math.max(minWeight, scaled);
            sum += finalWeight;
            list.add(new Entry(weight.reactionId, finalWeight));
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

    private static final class ReactionWeight {
        private final String reactionId;
        private final double logWeight;

        private ReactionWeight(String reactionId, double logWeight) {
            this.reactionId = reactionId;
            this.logWeight = logWeight;
        }
    }
}
