package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class QueryFeatureHelper {


    /**
     *
     *
     * Removes the following restricting query features:
     *
     * Atoms:
     *
     *
     * @param mol
     */
    public static void removeNarrowingQueryFeatures(ExtendedMolecule mol) {

        if(!mol.isFragment()) {
            return;
        }

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            mol.setAtomQueryFeature(atom, Molecule.cAtomQFNarrowing, false);
            mol.setAtomQueryFeature(atom, Molecule.cAtomQFNoMoreNeighbours,false);
        }
        for (int bond=0; bond<mol.getAllBonds(); bond++) {
            mol.setBondQueryFeature(bond, Molecule.cBondQFNarrowing, false);
        }

    }

}
