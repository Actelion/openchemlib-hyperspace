package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

public class HyperspaceUtilsOCL {

    /**
     *
     * @param mi
     * @return
     */
    public static List<Integer> findNeighbors(StereoMolecule mi, int a) {
        List<Integer> neighbors = new ArrayList<>();
        //mi.ensureHelperArrays(Molecule.cHelperNeighbours);
        for(int za=0;za<mi.getAllConnAtoms(a);za++) {
            neighbors.add(mi.getConnAtom(a,za));
        }
        return neighbors;
    }

    public static BitSet findBondsInBetweenAtoms(StereoMolecule mol, BitSet atoms) {
        mol.ensureHelperArrays(Molecule.cHelperCIP);
        BitSet bonds = new BitSet(mol.getBonds());
        // 1. determine bonds:
        for(int zi=0;zi<mol.getBonds();zi++) {
            int baa = mol.getBondAtom(0,zi);
            int bab = mol.getBondAtom(1,zi);
            if( atoms.get(baa) && atoms.get(bab) ) {
                bonds.set(zi,true);
            }
        }
        return bonds;
    }

}
