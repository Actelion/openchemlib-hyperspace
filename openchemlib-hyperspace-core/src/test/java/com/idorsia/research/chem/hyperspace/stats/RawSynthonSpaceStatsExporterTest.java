package com.idorsia.research.chem.hyperspace.stats;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.junit.jupiter.api.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class RawSynthonSpaceStatsExporterTest {
    @Test
    void computesExactReactionAndSourceStatsAndExamples() throws Exception {
        RawSynthonSpace rawSpace = toyRawSpace();
        RawSynthonSpaceStatsReport report = RawSynthonSpaceStatsExporter.analyze(rawSpace,
                RawSynthonSpaceStatsOptions.builder()
                        .examplesPerReaction(2)
                        .productSamplesPerReaction(3)
                        .threads(2)
                        .seed(7L)
                        .build());

        assertEquals(1, report.getSpaceSummary().reactionCount());
        assertEquals(1, report.getSpaceSummary().twoSetReactionCount());
        assertEquals("2", report.getSpaceSummary().totalProductCount().toString());
        assertEquals(3, report.getSpaceSummary().totalSynthons());
        assertEquals(1, report.getSourceStats().size());
        assertEquals("VendorA", report.getSourceStats().get(0).sourceSpace());
        assertEquals("2", report.getReactionStats().get(0).productCount().toString());
        assertEquals("0:2,1:1", report.getReactionStats().get(0).setSizes());
        assertEquals(3, report.getReactionStats().get(0).productSampleCount());
        assertFalse(report.getExampleProducts().isEmpty());
        assertTrue(report.getExampleProducts().size() <= 2);
    }

    @Test
    void writesReportBundle() throws Exception {
        RawSynthonSpaceStatsReport report = RawSynthonSpaceStatsExporter.analyze(toyRawSpace(),
                RawSynthonSpaceStatsOptions.builder()
                        .examplesPerReaction(1)
                        .productSamplesPerReaction(1)
                        .threads(1)
                        .build());
        Path out = Files.createTempDirectory("rawspace-stats");
        report.writeToDirectory(out);

        assertTrue(Files.exists(out.resolve("summary.md")));
        assertTrue(Files.exists(out.resolve("space_summary.tsv")));
        assertTrue(Files.exists(out.resolve("reaction_stats.tsv")));
        assertTrue(Files.exists(out.resolve("source_stats.tsv")));
        assertTrue(Files.exists(out.resolve("synthon_set_stats.tsv")));
        assertTrue(Files.exists(out.resolve("example_products.tsv")));
        assertTrue(Files.readString(out.resolve("example_products.tsv")).contains("Structure [idcode]"));
    }

    static RawSynthonSpace toyRawSpace() throws Exception {
        return RawSynthonSpace.builder("toy_stats")
                .addRawFragments("rxnA", 0, List.of(
                        raw("rxnA", 0, "a0", "C[U]"),
                        raw("rxnA", 0, "a1", "CC[U]")))
                .addRawFragments("rxnA", 1, List.of(
                        raw("rxnA", 1, "b0", "N[U]")))
                .addReactionMetadata("rxnA", "source.spaceName", "VendorA")
                .addReactionMetadata("rxnA", "source.spacePath", "/data/vendor.rawspace.gz")
                .addReactionMetadata("rxnA", "source.originalReactionId", "orig_rxnA")
                .build();
    }

    private static RawSynthon raw(String reactionId, int setIndex, String fragmentId, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule(reactionId, setIndex, fragmentId, molecule);
    }
}
