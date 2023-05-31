package com.idorsia.research.chem.hyperspace.importer;

import java.util.*;

public class CombinedMultiStepSynthonReaction implements MultiStepSynthonReaction {

    private List<MultiStepSynthonReaction> subReactions;

    public CombinedMultiStepSynthonReaction(ConcreteMultiStepSynthonReaction ri) {
        this.subReactions.add(ri);
    }

    private CombinedMultiStepSynthonReaction(List<MultiStepSynthonReaction> subReactions) {
        this.subReactions = subReactions;
    }

    @Override
    public boolean isConcrete() {
        return false;
    }

    @Override
    public Set<Integer> getOutConnectors() {
        int[] counts = new int[8];
        for(MultiStepSynthonReaction ri : subReactions) {
            for (ImporterSynthonSet si : ri.getAllSynthonSets()) {
                for (int xi : si.getConnectors()) {
                    counts[xi - 92]++;
                }
            }
        }
        Set<Integer> outc = new HashSet<>();
        for(int cci=0;cci<counts.length;cci++) {
            if(counts[cci]==1) { outc.add(counts[cci]);}
        }
        return outc;
    }

    @Override
    public Set<Integer> getAllUsedConnectors() {
        int[] counts = new int[8];
        for(MultiStepSynthonReaction ri : subReactions) {
            for (ImporterSynthonSet si : ri.getAllSynthonSets()) {
                for (int xi : si.getConnectors()) {
                    counts[xi - 92]++;
                }
            }
        }
        Set<Integer> used = new HashSet<>();
        for(int cci=0;cci<counts.length;cci++) {
            if(counts[cci]>0) { used.add(counts[cci]);}
        }
        return used;
    }

    @Override
    public List<ImporterSynthonSet> getAllSynthonSets() {
        List<ImporterSynthonSet> sets = new ArrayList<>();
        for(MultiStepSynthonReaction ri : subReactions) {
            sets.addAll(ri.getAllSynthonSets());
        }
        return sets;
    }

    @Override
    public void remapConnectors(ConnectorRemap remap) {
        for(MultiStepSynthonReaction cri : this.subReactions) {
            cri.remapConnectors(remap);
        }
    }

    /**
     * This function combines two multi step synthon reactions.
     *
     * It does the following things:
     * 1. compute both sets of out-connectors and check that they agree.
     * 2. find all internal connectors and make make them consistent.
     *
     *
     * So to create multistep synthon spaces:
     * Start at the leafs of the tree and combine synthon reactions upwards.
     *
     * This works as follows:
     * Assume we have A = a1/a2 that are one ConcreteMultistep and B=b1/b2 that is another ConcreteMultistep.
     * Now we process this as follows: we create CombinedMultistep from A, then we create CombinedMultistep from B,
     * then we use A.addMultiStepReaction(B).
     *
     * I.e. in general: start at the leafs of the tree and work upwards.
     *
     * @param ra
     * @param rb
     * @return
     */
    public static CombinedMultiStepSynthonReaction combineTwoMultiStepReactions(MultiStepSynthonReaction ra, MultiStepSynthonReaction rb) throws Exception {

        Set<Integer> out_a = ra.getOutConnectors();
        Set<Integer> out_b = rb.getOutConnectors();

        if(!out_a.equals(out_b)) {
            throw new Exception("Out Connectors do not agree..");
        }

        // collect all connectors:
        Set<Integer> all_a = ra.getAllUsedConnectors();
        Set<Integer> all_b = rb.getAllUsedConnectors();

        Set<Integer> internal_a = ra.getAllUsedConnectors(); internal_a.removeAll(out_a);
        Set<Integer> internal_b = rb.getAllUsedConnectors(); internal_b.removeAll(out_b);

        // check if we need remapping..
        if( Collections.disjoint(internal_a,internal_b) ) {
            // all good
        }
        else {
            // remap the ones in b
            List<Integer> available = new ArrayList<>();
            for(int zi=0;zi<8;zi++) {available.add(zi);}

            available.removeAll(out_a);
            available.removeAll(internal_a);

            // now the ones that are in b but not in available need remapping..
            Set<Integer> remap_needed = new HashSet<>( all_b );
            remap_needed.removeAll(available);

            List<Integer> still_available = new ArrayList<>(available);
            still_available.removeAll(all_b);

            // now we remap:
            Map<Integer,Integer> remap_b = new HashMap<>();
            int cnt_rm = 0;
            for(int ri_source : remap_needed) {
                remap_b.put(ri_source,still_available.get(0));
                still_available.remove(0);
                cnt_rm++;
            }
            System.out.println("[INFO] remap: "+remap_b.toString());

            // now we just have to apply the remap:
            rb.remapConnectors(new ConnectorRemap(remap_b));
        }

        List<MultiStepSynthonReaction> r_ci = new ArrayList<>();
        r_ci.add(ra);
        r_ci.add(rb);

        CombinedMultiStepSynthonReaction r_combined = new CombinedMultiStepSynthonReaction(r_ci);
        return r_combined;
    }

}
