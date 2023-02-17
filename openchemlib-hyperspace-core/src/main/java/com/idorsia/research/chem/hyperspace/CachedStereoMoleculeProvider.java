package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class CachedStereoMoleculeProvider {

    private Map<String, StereoMolecule> cache = new ConcurrentHashMap<>();

    public CachedStereoMoleculeProvider() {
    }

    Map<Thread,IDCodeParser> parsers = new HashMap<>();
    //private IDCodeParser icp = new IDCodeParser();

    private StereoMolecule parse(String idcode) {
        StereoMolecule nm = new StereoMolecule();
        Thread tc = Thread.currentThread();
        if(parsers.containsKey(tc)) {
            parsers.get(tc).parse(nm,idcode);
            return nm;
        }
        IDCodeParser nicp = new IDCodeParser();
        parsers.put(tc,nicp);
        nicp.parse(nm,idcode);
        return nm;
    }

    public StereoMolecule parseIDCode(String si) {
        if(cache.containsKey(si)){return new StereoMolecule(cache.get(si));}
        StereoMolecule nm = parse(si);
        cache.put(si,nm);
        return new StereoMolecule(nm);
    }

}
