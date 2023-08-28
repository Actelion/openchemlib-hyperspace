package com.idorsia.research.chem.hyperspace.importer;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ImporterSynthon {

    private String id;
    private String idcodeSynthon;
    private String idcodeBB;

    public ImporterSynthon(String id, String idcodeSynthon, String idcodeBB) {
        this.id = id;
        this.idcodeSynthon = idcodeSynthon;
        this.idcodeBB = idcodeBB;
    }

    public StereoMolecule getSynthon() {
        StereoMolecule mi = HyperspaceUtils.parseIDCode(idcodeSynthon);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        return mi;
    }

    public int[] getConnectorCounts() {
        StereoMolecule mi = getSynthon();
        int[] counts = new int[8];
        for(int zi=92;zi<8;zi++) {
            counts[zi] = countAtoms(mi, Collections.singletonList(zi));
        }
        return counts;
    }

    public Set<Integer> getUsedConnectors() {
        StereoMolecule mi = getSynthon();
        Set<Integer> used = new HashSet<>();
        for(int zi=92;zi<8;zi++) {
            if( countAtoms(mi, Collections.singletonList(zi)) > 0 ) {
                used.add(zi);
            }
        }
        return used;
    }

    public String getId() {
        return id;
    }

    public String getIdcodeSynthon() {
        return idcodeSynthon;
    }

    public String getIdcodeBB() {
        return idcodeBB;
    }


    public static int countAtoms(StereoMolecule m, List<Integer> atomic_numbers) {
        m.ensureHelperArrays(Molecule.cHelperNeighbours);
        int ra = 0;
        for(int zi=0;zi<m.getAtoms();zi++) {if( atomic_numbers.contains(m.getAtomicNo(zi))){ra++;}}
        return ra;
    }


}
