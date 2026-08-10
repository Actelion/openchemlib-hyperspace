package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class LocalOptimizationResultWriterTest {

    @TempDir
    Path tempDir;

    @Test
    void writesOptionalSourceSpaceColumn() throws Exception {
        Path output = tempDir.resolve("hits.tsv");
        LocalOptimizationResult result = new LocalOptimizationResult(
                "rxnA",
                List.of(new LocalOptimizationResult.BeamEntry(
                        List.of(new SynthonSpace.FragId("rxnA", 0, "idcode1", "frag1", new BitSet(), new BitSet(), new BitSet())),
                        0.72,
                        18,
                        4,
                        "assembled",
                        2)),
                List.of("seed1"));

        try (LocalOptimizationResultWriter writer = new LocalOptimizationResultWriter(output,
                reactionId -> "Vendor Space")) {
            writer.write(result);
        }

        List<String> lines = Files.readAllLines(output);
        assertEquals("rxnId	sourceSpace	fragIds	Structure [idcode]	atoms	rotatableBonds	phesaSimilarity	attemptIndex	seedFragIds", lines.get(0));
        assertEquals("rxnA	Vendor Space	frag1	assembled	18	4	0.72	2	seed1", lines.get(1));
    }
}
