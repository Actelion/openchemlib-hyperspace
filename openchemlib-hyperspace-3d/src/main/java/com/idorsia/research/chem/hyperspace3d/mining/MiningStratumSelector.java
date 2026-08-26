package com.idorsia.research.chem.hyperspace3d.mining;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Deterministic, bounded chemistry-stratum reservoirs over a mining universe. */
public final class MiningStratumSelector {
    public static final int DEFAULT_FLEXIBLE_QUOTA = 64;
    public static final int DEFAULT_HIGH_SP3_QUOTA = 64;
    public static final int DEFAULT_STEREOCHEMICAL_QUOTA = 64;
    public static final int DEFAULT_SIZE_EXTREME_QUOTA_PER_TAIL = 64;

    private MiningStratumSelector() {}

    public static Map<Integer, Long> defaults(List<MiningCandidateUniverse.Candidate> candidates,
            long seed) {
        Map<Integer, Long> result = new HashMap<>();
        add(result, sample(candidates, DEFAULT_FLEXIBLE_QUOTA, seed ^ 0x11L,
                value -> value.rotatableBonds() >= 6), MiningSelection.FLEXIBLE);
        add(result, sample(candidates, DEFAULT_HIGH_SP3_QUOTA, seed ^ 0x22L,
                value -> value.sp3Fraction() >= 0.5), MiningSelection.HIGH_SP3);
        add(result, sample(candidates, DEFAULT_STEREOCHEMICAL_QUOTA, seed ^ 0x33L,
                value -> value.stereoCenters() >= 1), MiningSelection.STEREOCHEMICAL);
        add(result, sample(candidates, DEFAULT_SIZE_EXTREME_QUOTA_PER_TAIL, seed ^ 0x44L,
                value -> value.heavyAtoms() <= 12), MiningSelection.SIZE_EXTREME);
        add(result, sample(candidates, DEFAULT_SIZE_EXTREME_QUOTA_PER_TAIL, seed ^ 0x55L,
                value -> value.heavyAtoms() >= 28), MiningSelection.SIZE_EXTREME);
        return result;
    }

    private static List<Integer> sample(List<MiningCandidateUniverse.Candidate> candidates,
            int quota, long seed, Criterion criterion) {
        List<Choice> choices = new ArrayList<>();
        for (var candidate : candidates) {
            if (criterion.test(candidate)) {
                choices.add(new Choice(candidate.id(), PredictedRankMiner.mix(seed, candidate.id(), 0)));
            }
        }
        choices.sort((left, right) -> {
            int order = Long.compareUnsigned(left.priority(), right.priority());
            return order != 0 ? order : Integer.compare(left.candidate(), right.candidate());
        });
        return choices.subList(0, Math.min(quota, choices.size())).stream()
                .map(Choice::candidate).sorted(Comparator.naturalOrder()).toList();
    }

    private static void add(Map<Integer, Long> result, List<Integer> candidates, long mask) {
        for (int candidate : candidates) result.merge(candidate, mask, (left, right) -> left | right);
    }

    private interface Criterion { boolean test(MiningCandidateUniverse.Candidate candidate); }
    private record Choice(int candidate, long priority) {}
}
