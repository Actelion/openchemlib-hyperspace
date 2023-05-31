package com.idorsia.research.chem.hyperspace.importer;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ConnectorRemap {

    public final Map<Integer,Integer> remap;
    //public final String remap_id;

    public ConnectorRemap() {
        this(new HashMap<>());
    }
    public ConnectorRemap(Map<Integer, Integer> remap) {
        this.remap = remap;
    }

    public static ImporterSynthonReaction remap(ConnectorRemap remap, ImporterSynthonReaction ra) {
        List<ImporterSynthonSet> isslist = ra.getSynthonSets();
        List<ImporterSynthonSet> isslist_remapped = new ArrayList<>();

        for(ImporterSynthonSet iss : isslist) {
            ArrayList<ImporterSynthon> synthons2 = new ArrayList<>();
            for( ImporterSynthon isi : iss.synthons ) {
                synthons2.add(remap(remap,isi));
            }
            isslist_remapped.add(new ImporterSynthonSet(iss.filePath,synthons2));
        }

        ImporterSynthonReaction ra2 = new ImporterSynthonReaction(isslist_remapped);
        return ra2;
    }

    public static ImporterSynthon remap(ConnectorRemap remap, ImporterSynthon is) {
        StereoMolecule mi = is.getSynthon();
        mi.ensureHelperArrays(Molecule.cHelperNeighbours);
        for(int zi=0;zi<mi.getAtoms();zi++) {
            int ai = mi.getAtomicNo(zi);
            if(remap.remap.containsKey( ai )) {
                mi.setAtomicNo(zi,remap.remap.get( ai ));
            }
        }
        ImporterSynthon is2 = new ImporterSynthon(is.getId(),mi.getIDCode(),is.getIdcodeBB());
        return is2;
    }

}
