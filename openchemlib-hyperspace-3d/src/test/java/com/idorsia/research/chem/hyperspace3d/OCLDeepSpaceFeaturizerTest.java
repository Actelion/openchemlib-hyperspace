package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import java.util.List;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class OCLDeepSpaceFeaturizerTest {
    @Test void constructsAndPadsBenzeneDeterministically() throws Exception {
        StereoMolecule benzene = molecule("c1ccccc1");
        OCLDeepSpaceFeaturizer f = new OCLDeepSpaceFeaturizer();
        FeaturizationResult first = f.featurize(benzene);
        FeaturizationResult second = f.featurize(benzene);
        assertTrue(first.accepted());
        assertArrayEquals(first.atomFeatures(), second.atomFeatures());
        assertArrayEquals(first.pairFeatures(), second.pairFeatures());
        assertArrayEquals(first.atomMask(), second.atomMask());
        assertEquals(32 * 50, first.atomFeatures().length);
        assertEquals(32 * 32 * 36, first.pairFeatures().length);
        for (int i = 0; i < 32; i++) assertEquals(i < 6, first.atomMask()[i]);
        assertEquals(1f, first.pairFeatures()[0]);
        assertEquals(1f, first.pairFeatures()[13]);
        int paddedAtom = 6 * 50;
        for (int i = paddedAtom; i < first.atomFeatures().length; i++) assertEquals(0f, first.atomFeatures()[i]);
        DeepSpaceTensorBatch batch = new DeepSpaceTensorBatchBuilder(f).build(List.of(benzene, benzene));
        assertEquals(2, batch.batchSize());
    }

    @Test void rejectsOversizedMoleculeWithStableReason() throws Exception {
        FeaturizationResult result = new OCLDeepSpaceFeaturizer().featurize(
                molecule("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
        assertFalse(result.accepted());
        assertEquals("TOO_MANY_HEAVY_ATOMS", result.rejectionReason());
    }

    static StereoMolecule molecule(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return molecule;
    }
}
