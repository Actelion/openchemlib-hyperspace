package com.idorsia.research.chem.hyperspace.reparametrization;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * 1. the synthon sets are processed according to a SynthonSetReparametrization
 *
 */
public class RxnReparametrization {
    private Map<Integer,SynthonSetReparametrization> synthonReparametrization = new HashMap<>();

    public SynthonSetReparametrization getSynthonSetReparametrizatino(int synthonSet) {
        return synthonReparametrization.get(synthonSet);
    }

    /**
     * Uses the sam synthon set reparametrization for all synthon sets
     *
     * @param space
     * @param rxn_id
     * @param reparametrizatino
     */
    public RxnReparametrization(SynthonSpace space, String rxn_id, SynthonSetReparametrization reparametrizatino) {
        synthonReparametrization = new HashMap<>();
        for(Map.Entry<Integer, SynthonSpace.FragType> e1 : space.getFragTypes(rxn_id).entrySet()) {
            synthonReparametrization.put(e1.getKey(),reparametrizatino);
        }
    }

}
