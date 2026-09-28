package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import org.json.JSONObject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class LatentAssemblyMinerTest {
    @TempDir
    Path tempDir;

    @Test
    void minesTwoOverlappingViewsWithExactCanonicalAtomMappings() throws Exception {
        LatentAssemblyMiner.Result result = mine(
                List.of(new Input("AB_1", "CCCOc1ccc(NCCC)cc1")),
                LatentAssemblyMiner.ArmMode.TWO);

        assertEquals(1, result.decompositions().size());
        LatentAssemblyDecomposition decomposition = result.decompositions().getFirst();
        assertEquals(2, decomposition.arms().size());
        assertEquals(2, decomposition.assemblyViews().size());
        assertEquals(2, decomposition.cutBonds().size());
        assertEquals(2, decomposition.scaffoldView().openScaffoldAttachmentViewAtoms().size());
        assertTrue(decomposition.productNonHydrogenAtomCount() <= 32);

        assertPartitionAndViewMappings(decomposition);
        for (int arm = 0; arm < decomposition.arms().size(); arm++) {
            Set<Integer> expected = new HashSet<>(decomposition.scaffoldProductAtoms());
            expected.addAll(decomposition.arms().get(arm).productAtoms());
            assertEquals(expected,
                    new HashSet<>(decomposition.assemblyViews().get(arm).productAtomByViewAtom()));
            assertEquals(1, decomposition.assemblyViews().get(arm)
                    .openScaffoldAttachmentViewAtoms().size());
        }
    }

    @Test
    void minesThreeArmStarWhenRequested() throws Exception {
        LatentAssemblyMiner.Result result = mine(
                List.of(new Input("ABC_1", "CCCc1cc(CCC)cc(CCC)c1")),
                LatentAssemblyMiner.ArmMode.THREE);

        assertEquals(1, result.decompositions().size());
        LatentAssemblyDecomposition decomposition = result.decompositions().getFirst();
        assertEquals(3, decomposition.arms().size());
        assertEquals(3, decomposition.cutBonds().size());
        assertEquals(3, decomposition.scaffoldView().openScaffoldAttachmentViewAtoms().size());
        decomposition.assemblyViews().forEach(view ->
                assertEquals(2, view.openScaffoldAttachmentViewAtoms().size()));
        assertPartitionAndViewMappings(decomposition);
    }

    @Test
    void noRingBondIsSplitAcrossComponents() throws Exception {
        LatentAssemblyMiner.Result result = mine(
                List.of(new Input("RINGS_1", "CCCc1ccc(-c2ccc(CCC)cc2)cc1")),
                LatentAssemblyMiner.ArmMode.TWO);
        LatentAssemblyDecomposition decomposition = result.decompositions().getFirst();
        StereoMolecule product = parseIdcode(decomposition.productIdcode());

        int[] componentByAtom = new int[product.getAtoms()];
        decomposition.scaffoldProductAtoms().forEach(atom -> componentByAtom[atom] = 1);
        for (int i = 0; i < decomposition.arms().size(); i++) {
            int component = i + 2;
            decomposition.arms().get(i).productAtoms().forEach(atom -> componentByAtom[atom] = component);
        }
        for (int bond = 0; bond < product.getBonds(); bond++) {
            if (product.isRingBond(bond)) {
                assertEquals(componentByAtom[product.getBondAtom(0, bond)],
                        componentByAtom[product.getBondAtom(1, bond)]);
            }
        }
    }

    @Test
    void rejectsProductsAboveHardThirtyTwoNonHydrogenAtomLimit() throws Exception {
        String thirtyThreeCarbonChain = "C".repeat(33);
        LatentAssemblyMiner.Result result = mine(List.of(
                new Input("VALID", "CCCOc1ccc(NCCC)cc1"),
                new Input("TOO_LARGE", thirtyThreeCarbonChain)
        ), LatentAssemblyMiner.ArmMode.TWO);

        assertEquals(2, result.stats().parsedMolecules());
        assertEquals(1, result.stats().rejectedProductsTooLarge());
        assertEquals(List.of("VALID"), result.decompositions().stream()
                .map(LatentAssemblyDecomposition::sourceMoleculeId).toList());

        assertThrows(IllegalArgumentException.class, () -> LatentAssemblyMiner.Options.builder()
                .input(tempDir.resolve("unused.tsv"))
                .maxProductNonHydrogenAtoms(33)
                .build());
    }

    @Test
    void writesSelfContainedGzippedJsonLines() throws Exception {
        LatentAssemblyDecomposition decomposition = mine(
                List.of(new Input("AB_JSON", "CCCOc1ccc(NCCC)cc1")),
                LatentAssemblyMiner.ArmMode.TWO).decompositions().getFirst();
        Path output = tempDir.resolve("assembly.jsonl.gz");
        try (LatentAssemblyJsonlWriter writer = new LatentAssemblyJsonlWriter(output)) {
            writer.write(decomposition);
        }

        String line;
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(Files.newInputStream(output)), StandardCharsets.UTF_8))) {
            line = reader.readLine();
            assertEquals(null, reader.readLine());
        }
        JSONObject json = new JSONObject(line);
        assertEquals(LatentAssemblyDecomposition.SCHEMA_VERSION, json.getString("schema_version"));
        assertEquals("AB_JSON", json.getString("source_molecule_id"));
        assertEquals(2, json.getJSONArray("arms").length());
        assertEquals(2, json.getJSONArray("assembly_views").length());
        assertEquals(2, json.getJSONArray("cut_bonds").length());
        assertTrue(json.getInt("product_non_hydrogen_atom_count") <= 32);
    }

    private LatentAssemblyMiner.Result mine(List<Input> inputs,
                                            LatentAssemblyMiner.ArmMode armMode) throws Exception {
        Path input = tempDir.resolve("input-" + armMode + ".tsv");
        List<String> lines = new ArrayList<>();
        lines.add("chembl_id\tidcode");
        for (Input value : inputs) {
            lines.add(value.id() + "\t" + idcode(value.smiles()));
        }
        Files.write(input, lines);
        return new LatentAssemblyMiner().mine(LatentAssemblyMiner.Options.builder()
                .input(input)
                .armMode(armMode)
                .build());
    }

    private static void assertPartitionAndViewMappings(LatentAssemblyDecomposition decomposition) {
        Set<Integer> partition = new HashSet<>(decomposition.scaffoldProductAtoms());
        decomposition.arms().forEach(arm -> {
            for (int atom : arm.productAtoms()) {
                assertTrue(partition.add(atom));
            }
        });
        assertEquals(decomposition.productAtomCount(), partition.size());

        StereoMolecule product = parseIdcode(decomposition.productIdcode());
        List<LatentAssemblyDecomposition.View> views = new ArrayList<>();
        views.add(decomposition.scaffoldView());
        views.addAll(decomposition.assemblyViews());
        for (LatentAssemblyDecomposition.View view : views) {
            StereoMolecule parsedView = parseIdcode(view.idcode());
            assertEquals(parsedView.getAtoms(), view.productAtomByViewAtom().size());
            for (int viewAtom = 0; viewAtom < parsedView.getAtoms(); viewAtom++) {
                int productAtom = view.productAtomByViewAtom().get(viewAtom);
                assertEquals(product.getAtomicNo(productAtom), parsedView.getAtomicNo(viewAtom));
            }
        }
    }

    private static String idcode(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule.getIDCode();
    }

    private static StereoMolecule parseIdcode(String idcode) {
        StereoMolecule molecule = new StereoMolecule();
        new IDCodeParser().parse(molecule, idcode);
        molecule.ensureHelperArrays(Molecule.cHelperRings);
        return molecule;
    }

    private record Input(String id, String smiles) {}
}
