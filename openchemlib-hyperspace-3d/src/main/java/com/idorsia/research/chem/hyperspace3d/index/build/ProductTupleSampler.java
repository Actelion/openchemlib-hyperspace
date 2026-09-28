package com.idorsia.research.chem.hyperspace3d.index.build;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.screening.ReactionScheduler;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Objects;

/**
 * Deterministic, weighted sampling of downsampled product tuples. Each reaction
 * uses a seeded affine permutation of its flattened Cartesian product, making
 * tuple generation exactly without replacement with O(reaction-count) state.
 */
public final class ProductTupleSampler {
    public static final String ALGORITHM = "weighted-downsampled-affine-without-replacement-v1";

    private final Map<String, ReactionState> states;
    private final List<String> reactionIds;
    private final Map<String, Double> weights;
    private final StatefulRandom random;

    public ProductTupleSampler(RawSynthonSpace full, RawSynthonSpace downsampled,
            long seed, ReactionScheduler.Weighting weighting) {
        Objects.requireNonNull(full, "full");
        Objects.requireNonNull(downsampled, "downsampled");
        Map<String, ReactionState> mutable = new LinkedHashMap<>();
        List<String> ids = new ArrayList<>(downsampled.getReactions().keySet());
        ids.sort(Comparator.naturalOrder());
        Map<String, List<Integer>> sizes = new LinkedHashMap<>();
        for (String reactionId : ids) {
            ReactionState state = buildState(full, downsampled, reactionId, seed);
            if (state != null) {
                mutable.put(reactionId, state);
                sizes.put(reactionId, state.poolSizes());
            }
        }
        if (mutable.isEmpty()) throw new IllegalArgumentException("no eligible downsampled reactions");
        ReactionScheduler scheduler = new ReactionScheduler(sizes, weighting);
        this.states = mutable;
        this.reactionIds = List.copyOf(mutable.keySet());
        this.weights = scheduler.getWeightsByReaction();
        this.random = new StatefulRandom(seed ^ 0x6a09e667f3bcc909L);
    }

    public List<String> reactionIds() { return reactionIds; }
    public long cardinality(String reactionId) { return required(reactionId).cardinality; }
    public boolean exhausted(String reactionId) { return required(reactionId).exhausted(); }

    public Sample next(String reactionId) {
        ReactionState state = required(reactionId);
        if (state.exhausted()) return null;
        return state.next();
    }

    public Sample nextWeighted() {
        double total = 0d;
        for (String id : reactionIds) if (!states.get(id).exhausted()) total += weights.get(id);
        if (!(total > 0d)) return null;
        double target = random.nextDouble() * total;
        double cumulative = 0d;
        for (String id : reactionIds) {
            ReactionState state = states.get(id);
            if (state.exhausted()) continue;
            cumulative += weights.get(id);
            if (target <= cumulative) return state.next();
        }
        for (int i = reactionIds.size() - 1; i >= 0; i--) {
            ReactionState state = states.get(reactionIds.get(i));
            if (!state.exhausted()) return state.next();
        }
        return null;
    }

    public Snapshot snapshot() {
        Snapshot snapshot = new Snapshot();
        snapshot.randomState = random.state;
        snapshot.reactions = new LinkedHashMap<>();
        states.forEach((id, state) -> snapshot.reactions.put(id,
                new Cursor(state.cursor, state.current)));
        return snapshot;
    }

    public void restore(Snapshot snapshot) {
        if (snapshot == null || snapshot.reactions == null
                || !snapshot.reactions.keySet().equals(states.keySet())) {
            throw new IllegalArgumentException("sampler checkpoint does not match reactions");
        }
        random.state = snapshot.randomState;
        snapshot.reactions.forEach((id, saved) -> {
            ReactionState state = states.get(id);
            if (saved.cursor < 0 || saved.cursor > state.cardinality
                    || saved.current < 0 || saved.current >= state.cardinality) {
                throw new IllegalArgumentException("invalid sampler cursor for " + id);
            }
            state.cursor = saved.cursor;
            state.current = saved.current;
        });
    }

    private ReactionState required(String id) {
        ReactionState state = states.get(id);
        if (state == null) throw new IllegalArgumentException("unknown reaction: " + id);
        return state;
    }

    private static ReactionState buildState(RawSynthonSpace full, RawSynthonSpace downsampled,
            String reactionId, long seed) {
        RawSynthonSpace.ReactionData fullData = full.getReactions().get(reactionId);
        RawSynthonSpace.ReactionData reducedData = downsampled.getReactions().get(reactionId);
        if (fullData == null) throw new IllegalArgumentException(
                "downsampled reaction missing from full rawspace: " + reactionId);
        if (!fullData.getRawFragmentSets().keySet().equals(
                reducedData.getRawFragmentSets().keySet())) {
            throw new IllegalArgumentException("downsampled reaction positions differ from full rawspace: "
                    + reactionId);
        }
        List<Integer> positions = new ArrayList<>(reducedData.getRawFragmentSets().keySet());
        positions.sort(Comparator.naturalOrder());
        if (positions.size() < 2) return null;
        if (positions.size() > 16) throw new IllegalArgumentException(
                "product-index shards support at most 16 synthon positions: " + reactionId);
        List<List<Choice>> pools = new ArrayList<>();
        long cardinality = 1;
        for (int position : positions) {
            List<RawSynthon> fullSet = fullData.getRawFragmentSets().get(position);
            List<RawSynthon> reducedSet = reducedData.getRawFragmentSets().get(position);
            if (fullSet == null || reducedSet == null || reducedSet.isEmpty()) {
                throw new IllegalArgumentException("missing synthon set " + reactionId + "/" + position);
            }
            Map<String, Integer> ordinalById = new LinkedHashMap<>();
            for (int i = 0; i < fullSet.size(); i++) {
                RawSynthon fullSynthon = fullSet.get(i);
                if (!reactionId.equals(fullSynthon.getReactionId())
                        || position != fullSynthon.getFragmentIndex()) {
                    throw new IllegalArgumentException("invalid full-rawspace synthon location: "
                            + fullSynthon.getFragmentId());
                }
                if (ordinalById.put(fullSet.get(i).getFragmentId(), i) != null) {
                    throw new IllegalArgumentException("duplicate full-rawspace synthon ID: "
                            + reactionId + "/" + fullSet.get(i).getFragmentId());
                }
            }
            List<Choice> choices = new ArrayList<>(reducedSet.size());
            Set<String> reducedIds = new HashSet<>();
            for (RawSynthon reduced : reducedSet) {
                if (!reactionId.equals(reduced.getReactionId())
                        || position != reduced.getFragmentIndex()) {
                    throw new IllegalArgumentException("invalid downsampled synthon location: "
                            + reduced.getFragmentId());
                }
                if (!reducedIds.add(reduced.getFragmentId())) {
                    throw new IllegalArgumentException("duplicate downsampled synthon ID: "
                            + reactionId + "/" + reduced.getFragmentId());
                }
                Integer ordinal = ordinalById.get(reduced.getFragmentId());
                if (ordinal == null) throw new IllegalArgumentException(
                        "downsampled synthon missing from full rawspace: " + reduced.getFragmentId());
                RawSynthon authoritative = fullSet.get(ordinal);
                if (!authoritative.getIdcode().equals(reduced.getIdcode())
                        || !authoritative.getConnectors().equals(reduced.getConnectors())) {
                    throw new IllegalArgumentException("downsampled synthon differs from full rawspace: "
                            + reactionId + "/" + reduced.getFragmentId());
                }
                choices.add(new Choice(authoritative, ordinal));
            }
            pools.add(List.copyOf(choices));
            try {
                cardinality = Math.multiplyExact(cardinality, choices.size());
            } catch (ArithmeticException overflow) {
                throw new IllegalArgumentException("downsampled reaction cardinality exceeds int64: " + reactionId);
            }
        }
        long key = mix64(seed ^ mixString(reactionId));
        long offset = floorMod(key, cardinality);
        long step = cardinality == 1 ? 0 : 1 + floorMod(mix64(key), cardinality - 1);
        while (cardinality > 1 && gcd(step, cardinality) != 1) {
            step = step + 1 == cardinality ? 1 : step + 1;
        }
        return new ReactionState(reactionId, List.copyOf(positions), List.copyOf(pools),
                cardinality, offset, step);
    }

    private static long addMod(long value, long increment, long modulus) {
        if (modulus == 1) return 0;
        return value >= modulus - increment ? value - (modulus - increment) : value + increment;
    }

    private static long floorMod(long value, long modulus) {
        return (value & Long.MAX_VALUE) % modulus;
    }

    private static long gcd(long a, long b) {
        while (b != 0) { long next = a % b; a = b; b = next; }
        return a;
    }

    private static long mixString(String value) {
        long hash = 0xcbf29ce484222325L;
        for (int i = 0; i < value.length(); i++) {
            hash ^= value.charAt(i);
            hash *= 0x100000001b3L;
        }
        return hash;
    }

    private static long mix64(long value) {
        value = (value ^ (value >>> 30)) * 0xbf58476d1ce4e5b9L;
        value = (value ^ (value >>> 27)) * 0x94d049bb133111ebL;
        return value ^ (value >>> 31);
    }

    public record Sample(ProductTuple tuple, List<RawSynthon> synthons) {
        public Sample { synthons = List.copyOf(synthons); }
    }

    public static final class Snapshot {
        public long randomState;
        public Map<String, Cursor> reactions;
    }

    public static final class Cursor {
        public long cursor;
        public long current;
        public Cursor() {}
        public Cursor(long cursor, long current) { this.cursor = cursor; this.current = current; }
    }

    private record Choice(RawSynthon synthon, int fullOrdinal) {}

    private static final class ReactionState {
        private final String reactionId;
        private final List<Integer> positions;
        private final List<List<Choice>> pools;
        private final long cardinality;
        private final long step;
        private long cursor;
        private long current;

        private ReactionState(String reactionId, List<Integer> positions, List<List<Choice>> pools,
                long cardinality, long current, long step) {
            this.reactionId = reactionId;
            this.positions = positions;
            this.pools = pools;
            this.cardinality = cardinality;
            this.current = current;
            this.step = step;
        }

        private List<Integer> poolSizes() {
            List<Integer> result = new ArrayList<>(pools.size());
            pools.forEach(pool -> result.add(pool.size()));
            return Collections.unmodifiableList(result);
        }

        private boolean exhausted() { return cursor >= cardinality; }

        private Sample next() {
            long flattened = current;
            current = addMod(current, step, cardinality);
            cursor++;
            List<String> ids = new ArrayList<>(pools.size());
            List<Integer> ordinals = new ArrayList<>(pools.size());
            List<RawSynthon> synthons = new ArrayList<>(pools.size());
            for (int i = pools.size() - 1; i >= 0; i--) {
                List<Choice> pool = pools.get(i);
                int choiceIndex = (int) (flattened % pool.size());
                flattened /= pool.size();
                Choice choice = pool.get(choiceIndex);
                ids.add(0, choice.synthon.getFragmentId());
                ordinals.add(0, choice.fullOrdinal);
                synthons.add(0, choice.synthon);
            }
            return new Sample(new ProductTuple(reactionId, ids, ordinals), synthons);
        }
    }

    private static final class StatefulRandom {
        private long state;
        private StatefulRandom(long seed) { this.state = seed; }
        private double nextDouble() {
            state += 0x9e3779b97f4a7c15L;
            return (mix64(state) >>> 11) * 0x1.0p-53;
        }
    }
}
