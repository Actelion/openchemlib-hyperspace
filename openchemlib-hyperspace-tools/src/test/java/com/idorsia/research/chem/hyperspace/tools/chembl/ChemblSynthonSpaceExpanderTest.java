package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ChemblSynthonSpaceExpanderTest {

    @TempDir
    Path tempDir;

    @Test
    void expandsObservedRawSpaceWithSmallSimilarSynthons() throws Exception {
        Path seedInput = tempDir.resolve("seed.tsv");
        Files.write(seedInput, List.of(
                "chembl_id\tidcode",
                "SEED_1\t" + idcode("CCOC"),
                "SEED_2\t" + idcode("CCCOC"),
                "SEED_3\t" + idcode("CCNCC")
        ));

        ObservedSynthonSpaceMiner.Result seed = new ObservedSynthonSpaceMiner().mine(
                ObservedSynthonSpaceMiner.Options.builder()
                        .input(seedInput)
                        .spaceName("seed")
                        .topReactions(5)
                        .minSetSize(1)
                        .minSourceMolecules(1)
                        .sampleProducts(0)
                        .build());

        Path catalogInput = tempDir.resolve("catalog.tsv");
        Files.write(catalogInput, List.of(
                "chembl_id\tidcode",
                "CAT_1\t" + idcode("CCOC"),
                "CAT_2\t" + idcode("CCCCOC"),
                "CAT_3\t" + idcode("CCCOCC"),
                "CAT_4\t" + idcode("CCNCCC"),
                "CAT_5\t" + idcode("CCSCC")
        ));

        ChemblSynthonSpaceExpander.Result expanded = new ChemblSynthonSpaceExpander().expand(seed.rawSpace(),
                ChemblSynthonSpaceExpander.Options.builder()
                        .chemblInput(catalogInput)
                        .catalogMaxMolecules(0)
                        .maxAddedPerSet(5)
                        .maxCandidateHeavyAtoms(12)
                        .minFragFpSimilarity(0.0)
                        .progressInterval(0)
                        .build());

        assertEquals("chembl-observed-synthon-space-expanded-v1",
                expanded.rawSpace().getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
        assertTrue(expanded.catalogStats().catalogFragments() > 0);
        assertTrue(expanded.reactions().stream().anyMatch(summary -> summary.addedSynthons() > 0));
        assertTrue(expanded.rawSpace().getReactions().values().stream()
                .flatMap(reaction -> reaction.getFragmentAttributes().values().stream())
                .anyMatch(attributes -> attributes.containsKey("expansion.score")));
    }


    @Test
    void capOnlyExpansionMinesSingleSplitsAndLeavesTwoConnectorSetsUnchanged() throws Exception {
        Path seedInput = tempDir.resolve("seed_caps.tsv");
        Files.write(seedInput, List.of(
                "chembl_id\tidcode",
                "SEED_1\t" + idcode("CCOC"),
                "SEED_2\t" + idcode("CCCOC"),
                "SEED_3\t" + idcode("CCNCC")
        ));

        ObservedSynthonSpaceMiner.Result seed = new ObservedSynthonSpaceMiner().mine(
                ObservedSynthonSpaceMiner.Options.builder()
                        .input(seedInput)
                        .spaceName("seed_caps")
                        .topReactions(5)
                        .minSetSize(1)
                        .minSourceMolecules(1)
                        .sampleProducts(0)
                        .build());

        Path catalogInput = tempDir.resolve("catalog_caps.tsv");
        Files.write(catalogInput, List.of(
                "chembl_id\tidcode",
                "CAT_1\t" + idcode("CCOC"),
                "CAT_2\t" + idcode("CCCCOC"),
                "CAT_3\t" + idcode("CCCOCC"),
                "CAT_4\t" + idcode("CCNCCC"),
                "CAT_5\t" + idcode("CCSCC")
        ));

        Map<String, Map<Integer, Integer>> beforeCounts = countsByReactionAndSet(seed.rawSpace());
        ChemblSynthonSpaceExpander.Result expanded = new ChemblSynthonSpaceExpander().expand(seed.rawSpace(),
                ChemblSynthonSpaceExpander.Options.builder()
                        .chemblInput(catalogInput)
                        .catalogMaxMolecules(0)
                        .catalogSplitFilter(ObservedSynthonSpaceMiner.SplitFilter.ONE)
                        .targetConnectorCount(ChemblSynthonSpaceExpander.TargetConnectorCount.ONE)
                        .maxAddedPerSet(5)
                        .maxCandidateHeavyAtoms(12)
                        .minFragFpSimilarity(0.0)
                        .progressInterval(0)
                        .build());

        assertEquals(0, expanded.catalogStats().twoConnectorFragments());
        assertTrue(expanded.reactions().stream().anyMatch(summary -> summary.addedSynthons() > 0));
        expanded.rawSpace().getReactions().forEach((reactionId, reaction) ->
                reaction.getRawFragmentSets().forEach((setIdx, synthons) -> {
                    if (connectorCount(synthons) == 2) {
                        assertEquals(beforeCounts.get(reactionId).get(setIdx), synthons.size());
                    }
                }));
    }

    @Test
    void scoreWeightsAreNormalizedAndSizeBiasCanFavorSmallerSynthons() {
        ChemblSynthonSpaceExpander.Options options = ChemblSynthonSpaceExpander.Options.builder()
                .chemblInput(tempDir.resolve("unused.tsv"))
                .similarityWeight(0.50)
                .smallnessWeight(0.45)
                .occurrenceWeight(0.05)
                .sizeBiasPower(1.5)
                .build();
        ChemblSynthonSpaceExpander.Options scaled = ChemblSynthonSpaceExpander.Options.builder()
                .chemblInput(tempDir.resolve("unused.tsv"))
                .similarityWeight(5.0)
                .smallnessWeight(4.5)
                .occurrenceWeight(0.5)
                .sizeBiasPower(1.5)
                .build();

        double score = ChemblSynthonSpaceExpander.score(options, 0.8, Math.pow(0.8, options.sizeBiasPower()), 0.2);
        double scaledScore = ChemblSynthonSpaceExpander.score(scaled, 0.8, Math.pow(0.8, scaled.sizeBiasPower()), 0.2);
        double smaller = ChemblSynthonSpaceExpander.score(options, 0.8, Math.pow(0.8, options.sizeBiasPower()), 0.2);
        double medium = ChemblSynthonSpaceExpander.score(options, 0.8, Math.pow(0.3, options.sizeBiasPower()), 0.2);

        assertEquals(score, scaledScore, 1.0e-12);
        assertTrue(smaller > medium);
    }

    private static Map<String, Map<Integer, Integer>> countsByReactionAndSet(RawSynthonSpace rawSpace) {
        Map<String, Map<Integer, Integer>> counts = new HashMap<>();
        rawSpace.getReactions().forEach((reactionId, reaction) -> {
            Map<Integer, Integer> setCounts = new HashMap<>();
            reaction.getRawFragmentSets().forEach((setIdx, synthons) -> setCounts.put(setIdx, synthons.size()));
            counts.put(reactionId, setCounts);
        });
        return counts;
    }

    private static int connectorCount(List<RawSynthon> synthons) {
        if (synthons.isEmpty()) {
            return 0;
        }
        return synthons.get(0).getConnectors().cardinality();
    }


    private static String idcode(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule.getIDCode();
    }
}
