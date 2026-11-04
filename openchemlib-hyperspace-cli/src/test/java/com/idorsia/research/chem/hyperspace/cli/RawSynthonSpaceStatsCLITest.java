package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import org.junit.jupiter.api.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

class RawSynthonSpaceStatsCLITest {
    @Test
    void cliWritesStatsBundle() throws Exception {
        RawSynthonSpace rawSpace = toyRawSpace();
        Path rawIn = Files.createTempFile("toy-stats", ".rawspace.gz");
        Path outDir = Files.createTempDirectory("toy-stats-out");
        RawSynthonSpaceIO.write(rawSpace, rawIn);

        RawSynthonSpaceStatsCLI.main(new String[]{
                "--rawIn", rawIn.toString(),
                "--outDir", outDir.toString(),
                "--examplesPerReaction", "1",
                "--productSamplesPerReaction", "1",
                "--threads", "1"
        });

        assertTrue(Files.exists(outDir.resolve("summary.md")));
        assertTrue(Files.readString(outDir.resolve("reaction_stats.tsv")).contains("rxnA"));
        assertTrue(Files.readString(outDir.resolve("source_stats.tsv")).contains("VendorA"));
    }

    private static RawSynthonSpace toyRawSpace() throws Exception {
        return RawSynthonSpace.builder("toy_stats")
                .addRawFragments("rxnA", 0, List.of(
                        raw("rxnA", 0, "a0", "C[U]"),
                        raw("rxnA", 0, "a1", "CC[U]")))
                .addRawFragments("rxnA", 1, List.of(
                        raw("rxnA", 1, "b0", "N[U]")))
                .addReactionMetadata("rxnA", "source.spaceName", "VendorA")
                .build();
    }

    private static RawSynthon raw(String reactionId, int setIndex, String fragmentId, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule(reactionId, setIndex, fragmentId, molecule);
    }
}
