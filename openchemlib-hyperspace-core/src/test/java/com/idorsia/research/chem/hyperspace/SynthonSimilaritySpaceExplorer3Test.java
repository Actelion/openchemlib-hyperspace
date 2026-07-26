package com.idorsia.research.chem.hyperspace;

import com.actelion.research.calc.combinatorics.CombinationGenerator;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

class SynthonSimilaritySpaceExplorer3Test {

    @Test
    void cappedCutCombinationsMatchExistingGeneratorBelowCap() {
        List<int[]> expected = CombinationGenerator.getAllOutOf(6, 3);
        List<int[]> actual = SynthonSimilaritySpaceExplorer3.generateCutCombinations(6, 3, 100);

        assertEquals(expected.size(), actual.size());
        for(int i = 0; i < expected.size(); i++) {
            assertArrayEquals(expected.get(i), actual.get(i));
        }
    }

    @Test
    void cappedCutCombinationsStopAtConfiguredLimit() {
        List<int[]> actual = SynthonSimilaritySpaceExplorer3.generateCutCombinations(30, 3, 5);

        assertEquals(5, actual.size());
        assertArrayEquals(new int[]{0, 1, 2}, actual.get(0));
        assertArrayEquals(new int[]{0, 1, 6}, actual.get(4));
    }

    @Test
    void oldSimilaritySearchConfigConstructorKeepsUnlimitedCutCombinations() {
        SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3 config =
                new SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3(4, 3, 3, 0.8, 4000);

        assertEquals(Integer.MAX_VALUE, config.max_cut_combinations_per_split_level);
    }

    @Test
    void newSimilaritySearchConfigConstructorStoresCutCombinationLimit() {
        SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3 config =
                new SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3(4, 3, 3, 0.8, 4000, 10_000);

        assertEquals(10_000, config.max_cut_combinations_per_split_level);
    }
}
