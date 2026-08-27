package com.idorsia.research.chem.hyperspace3d.mining;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

/**
 * Fixed-model four-stratum enrichment: 400 channel-balanced absolute-top,
 * 400 from the remaining top 1%, 400 from 1-5%, and 400 random candidates.
 */
public final class FixedBaselineEnrichmentMiner {
    public static final String POLICY_ID = "fixed_baseline_1600_v1";
    public static final int ABSOLUTE_TOP_STOP = 400;
    public static final double TOP_ONE_PERCENT = 0.01;
    public static final double TOP_FIVE_PERCENT = 0.05;
    public static final int[] CHANNEL_QUOTAS = {200, 100, 100};
    public static final int RANDOM_QUOTA = 400;
    public static final int TARGET_PAIRS = 1600;

    private final long seed;

    public FixedBaselineEnrichmentMiner(long seed) {
        this.seed = seed;
    }

    /** Scores are row-major [candidate][total, shape, pharmacophore]. */
    public List<PredictedRankMiner.Selection> select(float[][] scores) {
        validate(scores);
        int count = scores.length;
        if (count < 100_000) {
            throw new IllegalArgumentException(
                    "fixed-baseline enrichment requires at least 100,000 candidates");
        }
        int[][] ranks = new int[3][];
        for (int channel = 0; channel < ranks.length; channel++) {
            ranks[channel] = PredictedRankMiner.ranks(scores, channel);
        }
        Set<Integer> selected = new HashSet<>(TARGET_PAIRS * 2);
        Map<Integer, Long> masks = new HashMap<>(TARGET_PAIRS * 2);
        int topOneStop = Math.max(ABSOLUTE_TOP_STOP + 1,
                (int) Math.ceil(TOP_ONE_PERCENT * count));
        int topFiveStop = Math.max(topOneStop + 1,
                (int) Math.ceil(TOP_FIVE_PERCENT * count));
        selectRankStratum(ranks, 0, ABSOLUTE_TOP_STOP, 0, selected, masks);
        selectRankStratum(ranks, ABSOLUTE_TOP_STOP, topOneStop, 1, selected, masks);
        selectRankStratum(ranks, topOneStop, topFiveStop, 2, selected, masks);
        for (int candidate : reservoir(null, 0, count, RANDOM_QUOTA,
                PredictedRankMiner.mix(seed, 17, 31), selected)) {
            selected.add(candidate);
            masks.put(candidate, MiningSelection.RANDOM);
        }
        if (selected.size() != TARGET_PAIRS) {
            throw new IllegalStateException("fixed-baseline policy did not select 1,600 unique pairs");
        }
        List<PredictedRankMiner.Selection> result = new ArrayList<>(TARGET_PAIRS);
        for (int candidate : selected) {
            result.add(new PredictedRankMiner.Selection(candidate, scores[candidate].clone(),
                    ranks[0][candidate] + 1L, ranks[1][candidate] + 1L,
                    ranks[2][candidate] + 1L, masks.get(candidate)));
        }
        result.sort(Comparator.comparingInt(PredictedRankMiner.Selection::candidateOrdinal));
        return result;
    }

    private void selectRankStratum(int[][] ranks, int first, int stop, int bucket,
            Set<Integer> selected, Map<Integer, Long> masks) {
        int[] offsets = {MiningSelection.TOTAL_OFFSET, MiningSelection.SHAPE_OFFSET,
                MiningSelection.PHARMACOPHORE_OFFSET};
        for (int channel = 0; channel < ranks.length; channel++) {
            int before = selected.size();
            for (int candidate : reservoir(ranks[channel], first, stop,
                    CHANNEL_QUOTAS[channel], PredictedRankMiner.mix(seed, channel, bucket), selected)) {
                selected.add(candidate);
                masks.put(candidate, MiningSelection.rankBucket(offsets[channel], bucket));
            }
            if (selected.size() - before != CHANNEL_QUOTAS[channel]) {
                throw new IllegalStateException("rank stratum cannot satisfy its channel quota");
            }
        }
    }

    private static int[] reservoir(int[] rankByCandidate, int first, int stop, int quota,
            long seed, Set<Integer> excluded) {
        int count = rankByCandidate == null ? stop : rankByCandidate.length;
        PriorityQueue<Choice> heap = new PriorityQueue<>((left, right) -> {
            int priority = Long.compareUnsigned(right.priority, left.priority);
            return priority != 0 ? priority : Integer.compare(right.candidate, left.candidate);
        });
        for (int candidate = 0; candidate < count; candidate++) {
            if (excluded.contains(candidate)) continue;
            int rank = rankByCandidate == null ? candidate : rankByCandidate[candidate];
            if (rank < first || rank >= stop) continue;
            Choice choice = new Choice(candidate, PredictedRankMiner.mix(seed, candidate, rank));
            if (heap.size() < quota) heap.add(choice);
            else if (Long.compareUnsigned(choice.priority, heap.peek().priority) < 0
                    || (choice.priority == heap.peek().priority
                    && candidate < heap.peek().candidate)) {
                heap.poll(); heap.add(choice);
            }
        }
        if (heap.size() != quota) {
            throw new IllegalStateException("selection range cannot satisfy quota " + quota);
        }
        return heap.stream().mapToInt(Choice::candidate).sorted().toArray();
    }

    private static void validate(float[][] scores) {
        if (scores == null || scores.length == 0) {
            throw new IllegalArgumentException("scores are empty");
        }
        for (float[] row : scores) {
            if (row == null || row.length != 3 || !Float.isFinite(row[0])
                    || !Float.isFinite(row[1]) || !Float.isFinite(row[2])) {
                throw new IllegalArgumentException("scores must be finite Nx3 values");
            }
        }
    }

    private record Choice(int candidate, long priority) {}
}
