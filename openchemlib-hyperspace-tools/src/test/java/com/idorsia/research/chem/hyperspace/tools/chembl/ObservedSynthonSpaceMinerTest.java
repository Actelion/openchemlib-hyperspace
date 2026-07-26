package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ObservedSynthonSpaceMinerTest {

    @TempDir
    Path tempDir;

    @Test
    void minesSmallObservedRawSpaceAndRoundTripsJson() throws Exception {
        Path input = tempDir.resolve("chembl.tsv");
        Files.write(input, List.of(
                "chembl_id\tidcode",
                "MOL_1\t" + idcode("CCOC"),
                "MOL_2\t" + idcode("CCCOC"),
                "MOL_3\t" + idcode("CCNCC")
        ));

        ObservedSynthonSpaceMiner.Options options = ObservedSynthonSpaceMiner.Options.builder()
                .input(input)
                .spaceName("test_chembl_observed")
                .maxMolecules(0)
                .topReactions(5)
                .minSetSize(1)
                .minSourceMolecules(1)
                .sampleProducts(3)
                .build();

        ObservedSynthonSpaceMiner.Result result = new ObservedSynthonSpaceMiner().mine(options);

        assertTrue(result.stats().processedMolecules() > 0);
        assertTrue(result.stats().acceptedSplits() > 0);
        assertFalse(result.rawSpace().getReactions().isEmpty());
        result.rawSpace().getReactions().values().forEach(reaction -> {
            assertFalse(reaction.getRawFragmentSets().isEmpty());
            reaction.getRawFragmentSets().values().forEach(set -> assertFalse(set.isEmpty()));
        });

        Path output = tempDir.resolve("raw.json.gz");
        RawSynthonSpaceIO.write(result.rawSpace(), output);
        RawSynthonSpace read = RawSynthonSpaceIO.read(output);
        assertEquals(result.rawSpace().getReactions().size(), read.getReactions().size());
        assertEquals("chembl-observed-synthon-space-v1",
                read.getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
    }

    private static String idcode(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule.getIDCode();
    }
}
