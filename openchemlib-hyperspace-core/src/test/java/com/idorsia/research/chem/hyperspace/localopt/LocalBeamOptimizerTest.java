package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

class LocalBeamOptimizerTest {

    @Test
    void recordsStatsWithoutUnconditionalBestOutput() {
        RawSynthonSpace rawSpace = RawSynthonSpace.builder("toy")
                .addRawFragments("rxn", 0, List.of(raw("rxn", 0, "a0"), raw("rxn", 0, "b0")))
                .addRawFragments("rxn", 1, List.of(raw("rxn", 1, "a1"), raw("rxn", 1, "b1")))
                .build();
        LocalBeamOptimizer optimizer = new LocalBeamOptimizer(
                new SynthonSetAccessor(rawSpace),
                scoreByFragments(),
                sampleBetterNeighbor());
        LocalOptimizationRequest request = LocalOptimizationRequest.builder()
                .beamSize(4)
                .neighborPoolSize(2)
                .sampledNeighbors(1)
                .perPositionCap(4)
                .maxRounds(3)
                .patience(1)
                .minPhesaSimilarity(0.0)
                .reportAllCandidates(true)
                .randomSeed(1L)
                .logLevel(LocalOptimizationLogLevel.NONE)
                .build();

        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        PrintStream previousOut = System.out;
        LocalOptimizationResult result;
        try {
            System.setOut(new PrintStream(stdout, true, StandardCharsets.UTF_8));
            result = optimizer.optimize(new SeedAssembly("rxn", List.of("a0", "a1"), 0.5), request);
        } finally {
            System.setOut(previousOut);
        }

        String output = stdout.toString(StandardCharsets.UTF_8);
        assertFalse(output.contains("Best:"));
        LocalOptimizationResult.OptimizationStats stats = result.getStats();
        assertEquals(0.5, stats.getInitialScore(), 1e-12);
        assertEquals(0.8, stats.getBestScore(), 1e-12);
        assertEquals(2, stats.getRoundsAttempted());
        assertEquals(2, stats.getRoundsCompleted());
        assertEquals(4, stats.getScoredCandidates());
        assertEquals(LocalOptimizationResult.StopReason.PATIENCE, stats.getStopReason());
    }


    @Test
    void seedScoreFailureIsGatedByLogLevel() {
        RawSynthonSpace rawSpace = RawSynthonSpace.builder("toy")
                .addRawFragments("rxn", 0, List.of(raw("rxn", 0, "a0")))
                .build();
        LocalBeamOptimizer optimizer = new LocalBeamOptimizer(
                new SynthonSetAccessor(rawSpace),
                (reactionId, fragments, originatingRound, logger) -> null,
                sampleBetterNeighbor());

        String quietOutput = captureStdout(() -> optimizer.optimize(
                new SeedAssembly("rxn", List.of("a0"), 0.42),
                LocalOptimizationRequest.builder()
                        .maxRounds(1)
                        .logLevel(LocalOptimizationLogLevel.NONE)
                        .build()));
        assertEquals("", quietOutput);

        String summaryOutput = captureStdout(() -> optimizer.optimize(
                new SeedAssembly("rxn", List.of("a0"), 0.42),
                LocalOptimizationRequest.builder()
                        .maxRounds(1)
                        .logLevel(LocalOptimizationLogLevel.SUMMARY)
                        .build()));
        assertFalse(summaryOutput.contains("scoreCandidate returned null"));
        org.junit.jupiter.api.Assertions.assertTrue(summaryOutput.contains("stop=SEED_SCORE_FAILED"));
        org.junit.jupiter.api.Assertions.assertTrue(summaryOutput.contains("initial=0.4200"));
    }

    private static RawSynthon raw(String reactionId, int fragmentIndex, String fragmentId) {
        return new RawSynthon(reactionId, fragmentIndex, fragmentId, fragmentId + "_idcode", new BitSet());
    }

    private static AssemblyScorer scoreByFragments() {
        return (reactionId, fragments, originatingRound, logger) -> {
            double score = 0.5;
            for (SynthonSpace.FragId fragment : fragments) {
                if ("b0".equals(fragment.fragment_id)) {
                    score += 0.1;
                }
                if ("b1".equals(fragment.fragment_id)) {
                    score += 0.2;
                }
            }
            List<String> fragmentIds = new ArrayList<>();
            for (SynthonSpace.FragId fragment : fragments) {
                fragmentIds.add(fragment.fragment_id);
            }
            logger.logCandidate(fragmentIds, score, "idcode", originatingRound);
            return new LocalOptimizationResult.BeamEntry(fragments, score, 10, 1, "idcode", originatingRound);
        };
    }

    private static NeighborSampler sampleBetterNeighbor() {
        return (String reactionId,
                int fragIdx,
                SynthonSpace.FragId center,
                Map<Integer, List<SynthonSpace.FragId>> synthonSets,
                LocalOptimizationRequest request,
                Random rng) -> {
            for (SynthonSpace.FragId fragment : synthonSets.get(fragIdx)) {
                if (fragment.fragment_id.startsWith("b")) {
                    return List.of(fragment);
                }
            }
            return List.of();
        };
    }

    private static String captureStdout(Runnable runnable) {
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        PrintStream previousOut = System.out;
        try {
            System.setOut(new PrintStream(stdout, true, StandardCharsets.UTF_8));
            runnable.run();
        } finally {
            System.setOut(previousOut);
        }
        return stdout.toString(StandardCharsets.UTF_8);
    }
}
