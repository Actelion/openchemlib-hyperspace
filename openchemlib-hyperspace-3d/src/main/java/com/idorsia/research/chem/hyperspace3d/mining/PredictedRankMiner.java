package com.idorsia.research.chem.hyperspace3d.mining;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/** Deterministic independent total/shape/PP rank-bucket mining. */
public final class PredictedRankMiner {
    public static final double[] DEFAULT_BOUNDARIES = {0.0, 0.0001, 0.001, 0.01, 0.05, 0.20, 1.0};
    public static final int[] DEFAULT_QUOTAS = {128, 192, 192, 128, 96, 64};
    public static final int DEFAULT_RANDOM_QUOTA = 384;

    private final double[] boundaries;
    private final int[] quotas;
    private final int randomQuota;
    private final long seed;

    public PredictedRankMiner(double[] boundaries, int[] quotas, int randomQuota, long seed) {
        if (boundaries == null || quotas == null || boundaries.length != quotas.length + 1
                || quotas.length != MiningSelection.BUCKETS || boundaries[0] != 0.0
                || boundaries[boundaries.length - 1] != 1.0 || randomQuota < 0) {
            throw new IllegalArgumentException("invalid rank-bucket configuration");
        }
        for (int i = 0; i < quotas.length; i++) {
            if (!(boundaries[i] >= 0.0 && boundaries[i] < boundaries[i + 1]) || quotas[i] < 0) {
                throw new IllegalArgumentException("invalid rank-bucket configuration");
            }
        }
        this.boundaries = boundaries.clone();
        this.quotas = quotas.clone();
        this.randomQuota = randomQuota;
        this.seed = seed;
    }

    public static PredictedRankMiner defaults(long seed) {
        return new PredictedRankMiner(DEFAULT_BOUNDARIES, DEFAULT_QUOTAS,
                DEFAULT_RANDOM_QUOTA, seed);
    }

    /** Scores are row-major [candidate][total, shape, pharmacophore]. */
    public List<Selection> select(float[][] scores) {
        return select(scores, Map.of());
    }

    /** Merges independent retrieval-channel masks with rank buckets. */
    public List<Selection> select(float[][] scores, Map<Integer, Long> baseMasks) {
        if (scores == null || scores.length == 0) throw new IllegalArgumentException("scores are empty");
        int count = scores.length;
        for (float[] row : scores) {
            if (row == null || row.length != 3 || !Float.isFinite(row[0])
                    || !Float.isFinite(row[1]) || !Float.isFinite(row[2])) {
                throw new IllegalArgumentException("scores must be finite Nx3 values");
            }
        }
        int[][] ranks = new int[3][];
        Map<Integer, MutableSelection> selected = new HashMap<>();
        if (baseMasks != null) {
            for (var entry : baseMasks.entrySet()) {
                if (entry.getKey() < 0 || entry.getKey() >= count || entry.getValue() == null
                        || entry.getValue() == 0L) {
                    throw new IllegalArgumentException("invalid base selection mask");
                }
                selected.computeIfAbsent(entry.getKey(), MutableSelection::new).mask |= entry.getValue();
            }
        }
        int[] offsets = {MiningSelection.TOTAL_OFFSET, MiningSelection.SHAPE_OFFSET,
                MiningSelection.PHARMACOPHORE_OFFSET};
        for (int channel = 0; channel < 3; channel++) {
            ranks[channel] = ranks(scores, channel);
            for (int bucket = 0; bucket < quotas.length; bucket++) {
                int first = (int) Math.floor(boundaries[bucket] * count);
                int stop = bucket == quotas.length - 1 ? count
                        : (int) Math.floor(boundaries[bucket + 1] * count);
                for (int candidate : reservoir(ranks[channel], first, stop, quotas[bucket],
                        mix(seed, channel, bucket))) {
                    selected.computeIfAbsent(candidate, MutableSelection::new).mask |=
                            MiningSelection.rankBucket(offsets[channel], bucket);
                }
            }
        }
        int[] identity = new int[count];
        for (int index = 0; index < count; index++) identity[index] = index;
        for (int candidate : reservoir(identity, 0, count, randomQuota,
                mix(seed, 17, 31))) {
            selected.computeIfAbsent(candidate, MutableSelection::new).mask |= MiningSelection.RANDOM;
        }
        List<Selection> result = new ArrayList<>(selected.size());
        for (MutableSelection value : selected.values()) {
            int candidate = value.candidate;
            result.add(new Selection(candidate, scores[candidate].clone(),
                    ranks[0][candidate] + 1L, ranks[1][candidate] + 1L,
                    ranks[2][candidate] + 1L, value.mask));
        }
        result.sort(Comparator.comparingInt(Selection::candidateOrdinal));
        return result;
    }

    /** Returns zero-based rank indexed by candidate ordinal. */
    static int[] ranks(float[][] scores, int channel) {
        int count = scores.length;
        long[] ordered = new long[count];
        for (int index = 0; index < count; index++) {
            int bits = Float.floatToRawIntBits(scores[index][channel]);
            int sortable = bits ^ ((bits >> 31) & 0x7fff_ffff);
            ordered[index] = ((long) sortable << 32) | (0xffff_ffffL - index);
        }
        Arrays.sort(ordered);
        int[] result = new int[count];
        for (int rank = 0; rank < count; rank++) {
            int candidate = (int) (0xffff_ffffL - (ordered[count - rank - 1] & 0xffff_ffffL));
            result[candidate] = rank;
        }
        return result;
    }

    private static int[] reservoir(int[] rankByCandidateOrIdentity, int first, int stop,
            int quota, long seed) {
        if (quota == 0 || stop <= first) return new int[0];
        boolean ranks = !isIdentity(rankByCandidateOrIdentity);
        PriorityQueue<Choice> heap = new PriorityQueue<>((left, right) -> {
            int priority = compareUnsigned(right.priority, left.priority);
            return priority != 0 ? priority : Integer.compare(right.candidate, left.candidate);
        });
        for (int candidate = 0; candidate < rankByCandidateOrIdentity.length; candidate++) {
            int rank = ranks ? rankByCandidateOrIdentity[candidate] : candidate;
            if (rank < first || rank >= stop) continue;
            Choice choice = new Choice(candidate, mix(seed, candidate, rank));
            if (heap.size() < quota) heap.add(choice);
            else if (compareUnsigned(choice.priority, heap.peek().priority) < 0
                    || (choice.priority == heap.peek().priority
                    && candidate < heap.peek().candidate)) {
                heap.poll(); heap.add(choice);
            }
        }
        return heap.stream().mapToInt(Choice::candidate).sorted().toArray();
    }

    private static boolean isIdentity(int[] values) {
        for (int i = 0; i < values.length; i++) if (values[i] != i) return false;
        return true;
    }

    private static int compareUnsigned(long left, long right) {
        return Long.compare(left ^ Long.MIN_VALUE, right ^ Long.MIN_VALUE);
    }

    static long mix(long seed, long first, long second) {
        long value = seed ^ (first * 0x9E3779B97F4A7C15L)
                ^ (second * 0xD1B54A32D192ED03L);
        value ^= value >>> 30; value *= 0xBF58476D1CE4E5B9L;
        value ^= value >>> 27; value *= 0x94D049BB133111EBL;
        return value ^ (value >>> 31);
    }

    private record Choice(int candidate, long priority) {}
    private static final class MutableSelection {
        private final int candidate;
        private long mask;
        private MutableSelection(int candidate) { this.candidate = candidate; }
    }

    public record Selection(int candidateOrdinal, float[] scores, long rankTotal,
            long rankShape, long rankPharmacophore, long selectionMask) {}
}
