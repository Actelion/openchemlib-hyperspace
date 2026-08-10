package com.idorsia.research.chem.hyperspace.cleaning;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

class SynthonSpaceCleanerTest {

    @Test
    void maxSynthonAtomsDropsOversizedSynthonsBeforeCleaning() throws Exception {
        RawSynthon small1 = synthon("rxn", 0, "small-1", "[U]CC");
        RawSynthon large1 = synthon("rxn", 0, "large-1", "[U]CCCCCCCC");
        RawSynthon small2 = synthon("rxn", 0, "small-2", "[U]CN");
        RawSynthon large2 = synthon("rxn", 0, "large-2", "[U]CCCCCCCN");
        RawSynthon partner = synthon("rxn", 1, "partner", "[U]N");

        RawSynthonSpace raw = RawSynthonSpace.builder("test")
                .addRawFragments("rxn", 0, List.of(small1, large1, small2, large2))
                .addRawFragments("rxn", 1, List.of(partner))
                .addFragmentAttribute("rxn", "small-1", "source", "kept")
                .addFragmentAttribute("rxn", "large-1", "source", "dropped")
                .addFragmentAttribute("rxn", "small-2", "source", "kept")
                .build();

        SynthonCleaningOptions options = SynthonCleaningOptions.builder()
                .maxSynthonAtoms(4)
                .maxParallelism(3)
                .randomSeed(1L)
                .build();
        SynthonSpaceCleaner cleaner = new SynthonSpaceCleaner(new IdentityCleaner(), options);

        RawSynthonSpace cleaned = cleaner.clean(raw);
        RawSynthonSpace.ReactionData cleanedReaction = cleaned.getReactions().get("rxn");

        assertEquals("4", cleaned.getMetadata().get("cleaning.maxSynthonAtoms"));
        assertEquals(List.of("small-1", "small-2"), cleanedReaction.getRawFragmentSets().get(0).stream()
                .map(RawSynthon::getFragmentId)
                .toList());
        assertEquals(1, cleanedReaction.getRawFragmentSets().get(1).size());
        assertNull(cleanedReaction.findRawFragment("large-1"));
        assertNull(cleanedReaction.findRawFragment("large-2"));
        assertEquals("kept", cleanedReaction.getFragmentAttributes().get("small-1").get("source"));
        assertEquals("kept", cleanedReaction.getFragmentAttributes().get("small-2").get("source"));
        assertNull(cleanedReaction.getFragmentAttributes().get("large-1"));
    }

    private static RawSynthon synthon(String reactionId, int setIndex, String fragmentId, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return RawSynthon.fromMolecule(reactionId, setIndex, fragmentId, molecule);
    }

    private static final class IdentityCleaner implements StructureCleaner {
        @Override
        public List<StereoMolecule> cleanStructure(StereoMolecule molecule) {
            return List.of(new StereoMolecule(molecule));
        }
    }
}
