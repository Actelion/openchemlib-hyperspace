package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class CachedStereoMoleculeProvider {

    private Map<String, StereoMolecule> cache = new ConcurrentHashMap<>();

    public CachedStereoMoleculeProvider() {
    }

    private final ThreadLocal<IDCodeParser> parser = ThreadLocal.withInitial(IDCodeParser::new);

    private StereoMolecule parse(String idcode) {
        StereoMolecule nm = new StereoMolecule();
        parser.get().parse(nm, idcode);
        return nm;
    }

    public StereoMolecule parseIDCode(String si) {
        if(cache.containsKey(si)){return new StereoMolecule(cache.get(si));}
        StereoMolecule nm = parse(si);
        cache.put(si,nm);
        return new StereoMolecule(nm);
    }

}
