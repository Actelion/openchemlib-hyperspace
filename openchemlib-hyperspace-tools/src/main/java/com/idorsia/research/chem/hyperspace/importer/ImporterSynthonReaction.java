package com.idorsia.research.chem.hyperspace.importer;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ImporterSynthonReaction {

    private List<ImporterSynthonSet> synthonSets;

    public ImporterSynthonReaction(List<ImporterSynthonSet> synthonSets) {
        this.synthonSets = synthonSets;
    }

    public Set<Integer> getOutConnectors() {
        int[] counts = new int[8];
        for(ImporterSynthonSet si : synthonSets) {
            for( int xi : si.getConnectors()) {
                counts[ xi - 92 ]++;
            }
        }
        Set<Integer> outc = new HashSet<>();
        for(int cci=0;cci<counts.length;cci++) {
            if(counts[cci]==1) { outc.add(counts[cci]);}
        }
        return outc;
    }

    public Set<Integer> getAllUsedConnectors() {
        Set<Integer> all_used = synthonSets.stream().map( xi -> xi.getConnectors() ).reduce( (x,y) -> {x.addAll(y); return x;} ).get();
        return all_used;
    }

    /**
     * Assumes: that synthon sets are correct, i.e. all have synthon in each set have same number of connectors.
     *
     * Checks that:
     * 1. all connectors appear 0,1 or 2 times.
     *
     * @throws Exception
     */
    public void checkValidity() throws Exception {
        int[] counts = new int[8];
        for(ImporterSynthonSet si : synthonSets) {
            for( int xi : si.getConnectors()) {
                counts[ xi - 92 ]++;
            }
        }

        // check 0,1,2:
        for(int cci=0;cci<counts.length;cci++) {
            int count_i = counts[cci];
            if(count_i==0 || count_i==1 || count_i==2 ) {
                // ok
            }
            else {
                System.out.println("[WARN] found error for connector: "+cci+" -> appears "+count_i+" times..");
            }
        }
    }

    public List<ImporterSynthonSet> getSynthonSets() {
        return this.synthonSets;
    }

}
