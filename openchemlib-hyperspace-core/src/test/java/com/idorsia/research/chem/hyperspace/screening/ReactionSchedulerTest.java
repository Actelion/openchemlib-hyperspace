package com.idorsia.research.chem.hyperspace.screening;

import org.junit.jupiter.api.Test;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class ReactionSchedulerTest {

    @Test
    void bucketedProductWeightsUseExclusiveUpperBounds() {
        Map<String, List<Integer>> sizes = new LinkedHashMap<>();
        sizes.put("small", List.of(20, 99));      // 1980
        sizes.put("medium-low", List.of(20, 100)); // 2000
        sizes.put("medium-high", List.of(100, 499)); // 49900
        sizes.put("large", List.of(100, 500));    // 50000

        ReactionScheduler scheduler = new ReactionScheduler(sizes, ReactionScheduler.Weighting.bucketedProduct(List.of(
                new ReactionScheduler.Bucket(2000.0, 0.01),
                new ReactionScheduler.Bucket(50000.0, 0.1),
                new ReactionScheduler.Bucket(null, 1.0)
        )));

        Map<String, Double> weights = scheduler.getWeightsByReaction();
        assertEquals(0.01, weights.get("small"), 1e-12);
        assertEquals(0.1, weights.get("medium-low"), 1e-12);
        assertEquals(0.1, weights.get("medium-high"), 1e-12);
        assertEquals(1.0, weights.get("large"), 1e-12);
    }

    @Test
    void exponentSchedulerConstructorStillWorks() {
        Map<String, List<Integer>> sizes = new LinkedHashMap<>();
        sizes.put("a", List.of(10, 10));
        sizes.put("b", List.of(100, 100));

        ReactionScheduler scheduler = new ReactionScheduler(sizes, 0.01, 1.0);

        assertEquals(0.01, scheduler.getWeightsByReaction().get("a"), 1e-12);
        assertEquals(1.0, scheduler.getWeightsByReaction().get("b"), 1e-12);
    }

    @Test
    void bucketedProductRequiresOpenEndedFinalBucket() {
        assertThrows(IllegalArgumentException.class, () -> ReactionScheduler.Weighting.bucketedProduct(List.of(
                new ReactionScheduler.Bucket(2000.0, 0.01),
                new ReactionScheduler.Bucket(50000.0, 0.1)
        )));
    }
}
