package com.idorsia.research.chem.hyperspace.screening;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Weighted scheduler for selecting reactions during continuous screening.
 */
public final class ReactionScheduler {

    private final List<Entry> entries;
    private final Map<String, Double> weightsByReaction;
    private final double totalWeight;

    public ReactionScheduler(Map<String, List<Integer>> reactionSizes,
                             double minWeight,
                             double exponent) {
        this(reactionSizes, Weighting.exponent(minWeight, exponent));
    }

    public ReactionScheduler(Map<String, List<Integer>> reactionSizes,
                             Weighting weighting) {
        if (reactionSizes.isEmpty()) {
            throw new IllegalArgumentException("No reactions available for scheduling");
        }
        Weighting effectiveWeighting = weighting == null ? Weighting.exponent(0.01, 1.0) : weighting;
        List<Entry> list = new ArrayList<>();
        Map<String, Double> weightMap = new LinkedHashMap<>();
        double sum = effectiveWeighting.apply(reactionSizes, list, weightMap);
        this.entries = Collections.unmodifiableList(list);
        this.weightsByReaction = Collections.unmodifiableMap(weightMap);
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

    public Map<String, Double> getWeightsByReaction() {
        return weightsByReaction;
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

    public static final class Weighting {
        private final Mode mode;
        private final double minWeight;
        private final double exponent;
        private final List<Bucket> buckets;

        private Weighting(Mode mode, double minWeight, double exponent, List<Bucket> buckets) {
            this.mode = mode;
            this.minWeight = minWeight;
            this.exponent = exponent;
            this.buckets = buckets == null ? Collections.emptyList() : Collections.unmodifiableList(new ArrayList<>(buckets));
            validate();
        }

        public static Weighting exponent(double minWeight, double exponent) {
            return new Weighting(Mode.EXPONENT, minWeight, exponent, Collections.emptyList());
        }

        public static Weighting bucketedProduct(List<Bucket> buckets) {
            return new Weighting(Mode.BUCKETED_PRODUCT, 0.0, 1.0, buckets);
        }

        private void validate() {
            if (mode == null) {
                throw new IllegalArgumentException("Reaction weighting mode is required");
            }
            if (mode == Mode.EXPONENT) {
                if (exponent <= 0.0 || !Double.isFinite(exponent)) {
                    throw new IllegalArgumentException("Exponent must be positive and finite");
                }
                if (minWeight < 0.0 || !Double.isFinite(minWeight)) {
                    throw new IllegalArgumentException("Minimum weight must be >= 0 and finite");
                }
                return;
            }
            if (buckets.isEmpty()) {
                throw new IllegalArgumentException("Bucketed reaction weighting requires at least one bucket");
            }
            double previous = 0.0;
            boolean sawFallback = false;
            for (int idx = 0; idx < buckets.size(); idx++) {
                Bucket bucket = buckets.get(idx);
                if (bucket == null) {
                    throw new IllegalArgumentException("Reaction weighting bucket must not be null");
                }
                if (bucket.weight < 0.0 || !Double.isFinite(bucket.weight)) {
                    throw new IllegalArgumentException("Reaction weighting bucket weight must be >= 0 and finite");
                }
                if (bucket.maxProductExclusive == null) {
                    if (idx != buckets.size() - 1) {
                        throw new IllegalArgumentException("Open-ended reaction weighting bucket must be last");
                    }
                    sawFallback = true;
                    continue;
                }
                if (!Double.isFinite(bucket.maxProductExclusive) || bucket.maxProductExclusive <= previous) {
                    throw new IllegalArgumentException("Reaction weighting bucket thresholds must be finite and increasing");
                }
                previous = bucket.maxProductExclusive;
            }
            if (!sawFallback) {
                throw new IllegalArgumentException("Bucketed reaction weighting requires a final open-ended bucket");
            }
        }

        private double apply(Map<String, List<Integer>> reactionSizes,
                             List<Entry> entries,
                             Map<String, Double> weightsByReaction) {
            if (mode == Mode.BUCKETED_PRODUCT) {
                return applyBucketed(reactionSizes, entries, weightsByReaction);
            }
            return applyExponent(reactionSizes, entries, weightsByReaction);
        }

        private double applyExponent(Map<String, List<Integer>> reactionSizes,
                                     List<Entry> entries,
                                     Map<String, Double> weightsByReaction) {
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
                entries.add(new Entry(weight.reactionId, finalWeight));
                weightsByReaction.put(weight.reactionId, finalWeight);
            }
            return sum;
        }

        private double applyBucketed(Map<String, List<Integer>> reactionSizes,
                                     List<Entry> entries,
                                     Map<String, Double> weightsByReaction) {
            double sum = 0d;
            for (Map.Entry<String, List<Integer>> entry : reactionSizes.entrySet()) {
                double product = product(entry.getValue());
                double weight = bucketWeight(product);
                sum += weight;
                entries.add(new Entry(entry.getKey(), weight));
                weightsByReaction.put(entry.getKey(), weight);
            }
            if (sum <= 0.0) {
                throw new IllegalArgumentException("Total reaction weight must be positive");
            }
            return sum;
        }

        private double bucketWeight(double product) {
            for (Bucket bucket : buckets) {
                if (bucket.maxProductExclusive == null || product < bucket.maxProductExclusive) {
                    return bucket.weight;
                }
            }
            return buckets.get(buckets.size() - 1).weight;
        }

        private double product(List<Integer> sizes) {
            double product = 1.0;
            for (int size : sizes) {
                if (size <= 0) {
                    return 0.0;
                }
                product *= size;
                if (Double.isInfinite(product)) {
                    return Double.POSITIVE_INFINITY;
                }
            }
            return product;
        }
    }

    public static final class Bucket {
        private final Double maxProductExclusive;
        private final double weight;

        public Bucket(Double maxProductExclusive, double weight) {
            this.maxProductExclusive = maxProductExclusive;
            this.weight = weight;
        }

        public Double getMaxProductExclusive() {
            return maxProductExclusive;
        }

        public double getWeight() {
            return weight;
        }
    }

    public enum Mode {
        EXPONENT,
        BUCKETED_PRODUCT
    }
}
