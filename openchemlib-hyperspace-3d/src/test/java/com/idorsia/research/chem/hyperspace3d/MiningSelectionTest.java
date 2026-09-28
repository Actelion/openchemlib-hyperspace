package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.mining.CandidateReferenceSampler;
import com.idorsia.research.chem.hyperspace3d.mining.MiningCandidateUniverse;
import com.idorsia.research.chem.hyperspace3d.mining.MiningSelection;
import com.idorsia.research.chem.hyperspace3d.mining.MiningStratumSelector;
import com.idorsia.research.chem.hyperspace3d.mining.PredictedRankMiner;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import org.junit.jupiter.api.Test;

class MiningSelectionTest {
    @Test
    void bottomHashReferenceSamplingIsBoundedStableAndUnique() {
        var first = new CandidateReferenceSampler(100, 17);
        var second = new CandidateReferenceSampler(100, 17);
        for (int row = 0; row < 10_000; row++) first.offer(row % 3, row);
        for (int row = 9_999; row >= 0; row--) second.offer(row % 3, row);
        assertEquals(first.results(), second.results());
        assertEquals(100, first.results().size());
        assertEquals(100, new HashSet<>(first.results().stream()
                .map(CandidateReferenceSampler.SampledReference::reference).toList()).size());
    }

    @Test
    void ranksAllFiniteFloatsNumericallyAndMergesChemistryChannels() {
        float[][] scores = {
                {-2.0f, 0.1f, 0.1f}, {0.0f, 0.2f, 0.2f},
                {-1.0f, 0.3f, 0.3f}, {2.0f, 0.4f, 0.4f}
        };
        var miner = new PredictedRankMiner(
                new double[]{0, .25, .5, .75, .875, .9375, 1},
                new int[]{1, 1, 1, 1, 1, 1}, 0, 17);
        var selections = miner.select(scores, java.util.Map.of(0, MiningSelection.FLEXIBLE));
        var ranks = selections.stream().collect(java.util.stream.Collectors.toMap(
                PredictedRankMiner.Selection::candidateOrdinal,
                PredictedRankMiner.Selection::rankTotal));
        assertEquals(1L, ranks.get(3));
        assertEquals(2L, ranks.get(1));
        assertEquals(3L, ranks.get(2));
        assertEquals(4L, ranks.get(0));
        assertTrue(selections.stream().filter(value -> value.candidateOrdinal() == 0)
                .findFirst().orElseThrow().selectionMask() != 0L);
    }

    @Test
    void chemistryStrataAreDeterministicAndCanOverlap() {
        List<MiningCandidateUniverse.Candidate> values = new ArrayList<>();
        for (int index = 0; index < 100; index++) {
            values.add(new MiningCandidateUniverse.Candidate(index,
                    new MoleculeVectorReference(0, index), index, "m" + index, "CC",
                    index < 50 ? 10 : 30, index % 8, index / 100.0, index % 3));
        }
        var first = MiningStratumSelector.defaults(values, 17);
        assertEquals(first, MiningStratumSelector.defaults(values, 17));
        assertTrue(first.values().stream().anyMatch(mask ->
                Long.bitCount(mask & (MiningSelection.FLEXIBLE | MiningSelection.HIGH_SP3
                        | MiningSelection.STEREOCHEMICAL | MiningSelection.SIZE_EXTREME)) > 1));
    }
}
