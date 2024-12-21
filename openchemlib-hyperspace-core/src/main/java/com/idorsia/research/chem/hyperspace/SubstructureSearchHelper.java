package com.idorsia.research.chem.hyperspace;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class SubstructureSearchHelper {


    /**
     *
     * @param space
     * @param cdp
     * @param mi
     * @param threads
     * @param fillIncompleteMappings if true, then adds all synthons of the non-matched synthon sets into the
     *                               combinatorial hits.
     * @param omitRxnsWithHitsFromLowerSplitNumber if true, then if we find hits at a specific number of splits, we omit
     *                                             the synthon reaction in more complicated searches.
     * @return
     */
    public static List<SynthonSpace.CombinatorialHit> run_substructure_search_01(SynthonSpace space, CachedDescriptorProvider cdp, StereoMolecule mi, int threads, boolean fillIncompleteMappings, boolean omitRxnsWithHitsFromLowerSplitNumber) {

        List<SynthonSpace.CombinatorialHit> hits_no_splits = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_1s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_2s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_3s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_4s = new ArrayList<>();

        List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space,cdp,mi,threads);
        List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();

        Set<String> discovered_rxns = new HashSet<>();

        if(false) {
            if (!hits_bbs.isEmpty()) {
                BitSet fp_molecule = cdp.getFP_cached(mi);

                // 1. sort expanded hits by frag type.. (then we report all hits from the same frag type in the same combinatorial hit..)
                Map<SynthonSpace.FragType, List<SynthonSpace.ExpandedHit>> sorted_bb_hits = new HashMap<>();
                for (SynthonSpace.ExpandedHit hi : hits_bbs) {
                    Integer hit_frag = hi.hit_fragments.keySet().iterator().next();
                    SynthonSpace.FragId hit_frag_id = hi.hit_fragments.get(hit_frag);

                    SynthonSpace.FragType hit_ft = new SynthonSpace.FragType(hi.rxn, hit_frag_id.frag);
                    if (!sorted_bb_hits.containsKey(hit_ft)) {
                        sorted_bb_hits.put(hit_ft, new ArrayList<>());
                    }
                    sorted_bb_hits.get(hit_ft).add(hi);
                }

                for (SynthonSpace.FragType fti_a : sorted_bb_hits.keySet()) {

                    List<SynthonSpace.FragId> all_matching_fragids = new ArrayList<>();
                    for (SynthonSpace.ExpandedHit hi : sorted_bb_hits.get(fti_a)) {
                        all_matching_fragids.add(hi.hit_fragments.get(fti_a.frag));
                    }

                    //SynthonSpace.CombinatorialHit chi = new SynthonSpace.CombinatorialHit();

                    Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> frag_map = new HashMap<>();
                    Map<Integer, Pair<SynthonSpace.FragType, BitSet>> fp_mapping = new HashMap<>();

                    StereoMolecule sri_fragments[] = new StereoMolecule[space.getFragTypes(fti_a.rxn_id).size()];
                    StereoMolecule mi_a = new StereoMolecule(mi);
                    mi_a.ensureHelperArrays(Molecule.cHelperCIP);
                    for (int szi = 0; szi < sri_fragments.length; szi++) {
                        if (szi == 0) {
                            sri_fragments[szi] = mi_a;
                        } else {
                            sri_fragments[szi] = new StereoMolecule();
                        }
                    }
                    SynthonShredder.SplitResult sri = new SynthonShredder.SplitResult(sri_fragments, new ArrayList<>(), "", new ArrayList<>());

                    // loop over frag types and fill the required maps:
                    int other_splits_cnt = 1;
                    for (Integer ftii : space.getFragTypes(fti_a.rxn_id).keySet()) {
                        SynthonSpace.FragType fti = space.getFragTypes(fti_a.rxn_id).get(ftii);
                        if (ftii == fti_a.frag) {
                            fp_mapping.put(0, Pair.of(fti_a, fp_molecule));
                            frag_map.put(fti_a, all_matching_fragids);
                        }
//                    else {
//                        // put all fragments into combinatorial hit:
//                        frag_map.put( fti , new ArrayList<>( space.getSynthonSet(fti_a.rxn_id,ftii) ) );
//                        // put empty bitset into fp_mapping
//                        fp_mapping.put( other_splits_cnt , Pair.of( fti , new BitSet( space.getBits() ) ) );
//                        other_splits_cnt++;
//                    }
                    }
                    SynthonSpace.CombinatorialHit hit_bb_ci = new SynthonSpace.CombinatorialHit(fti_a.rxn_id, frag_map, sri, fp_mapping);
                    hits_bbs_as_combinatorial_hits.add(hit_bb_ci);
                }
            }
        }
        if(true) {
            List<SynthonSpace.CombinatorialHit> hits_from_bbs = new ArrayList<>();
            expandBuildingBlockExpandedHits(space, cdp, hits_bbs, mi, new SynthonSpace.CombinatorialHitReceiver() {
                @Override
                public synchronized void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                    hits_from_bbs.addAll(hi);
                }
            });
            hits_bbs_as_combinatorial_hits = hits_from_bbs;
        }




        discovered_rxns.addAll(hits_bbs_as_combinatorial_hits.stream().map(hi -> hi.rxn).collect(Collectors.toList()));

        hits_1s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,1,2,
                (omitRxnsWithHitsFromLowerSplitNumber?discovered_rxns:new HashSet<>()), 200,1);
        discovered_rxns.addAll(hits_1s.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
        hits_2s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,2,3,
                (omitRxnsWithHitsFromLowerSplitNumber?discovered_rxns:new HashSet<>()), 200,1);
        discovered_rxns.addAll(hits_2s.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
        hits_3s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,3,4,
                (omitRxnsWithHitsFromLowerSplitNumber?discovered_rxns:new HashSet<>()),200,1);
        discovered_rxns.addAll( hits_3s.stream().map( hi -> hi.rxn ).collect(Collectors.toList()) );
        //if(hits_3s.isEmpty()) {
        //    hits_4s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,4,4,200,1);
        //}

        System.out.println("Hits (bb): "+hits_bbs_as_combinatorial_hits.size());
        System.out.println("Hits (1s): "+hits_1s.size());
        System.out.println("Hits (2s): "+hits_2s.size());
        System.out.println("Hits (3s): "+hits_3s.size());
        System.out.println("Hits (4s): "+hits_4s.size());

        List<SynthonSpace.CombinatorialHit> all_hits = new ArrayList<>();
        all_hits.addAll(hits_bbs_as_combinatorial_hits);
        all_hits.addAll(hits_1s);
        all_hits.addAll(hits_2s);
        all_hits.addAll(hits_3s);
        all_hits.addAll(hits_4s);

        if(fillIncompleteMappings) {
            for(int zi=0;zi<all_hits.size();zi++) {
                SynthonSpace.CombinatorialHit ci = all_hits.get(zi);
                // check if all synthon sets of the synthon reaction are matched. If not, add the fragments to the combinatorial hit
                List<SynthonSpace.FragType> synthsets_matched = new ArrayList<>(ci.mapping.values().stream().map( xi -> (xi!=null)?xi.getLeft(): null).collect(Collectors.toList()));
                List<SynthonSpace.FragType> extra_synthsets = space.getFragTypes(ci.rxn).values().stream().collect(Collectors.toList());
                extra_synthsets.removeAll(synthsets_matched);
                if(!extra_synthsets.isEmpty()) {
                    for(SynthonSpace.FragType fti : extra_synthsets) {
                        System.out.println("Fill incomplete mapping for set: "+fti.toString());
                        ci.hit_fragments.put(fti,space.getSynthonSet(fti.rxn_id,fti.frag));
                    }
                }
            }
        }

        return all_hits;
    }


    public static void expandBuildingBlockExpandedHits(SynthonSpace space, CachedDescriptorProvider cdp, List<SynthonSpace.ExpandedHit> hits_bbs, StereoMolecule mi, SynthonSpace.CombinatorialHitReceiver receiver) {
        //List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();
        //Set<String> discovered_rxns = new HashSet<>();

        if(!hits_bbs.isEmpty()) {
            BitSet fp_molecule = cdp.getFP_cached(mi);

            // 1. sort expanded hits by frag type.. (then we report all hits from the same frag type in the same combinatorial hit..)
            Map<SynthonSpace.FragType,List<SynthonSpace.ExpandedHit>> sorted_bb_hits = new HashMap<>();
            for (SynthonSpace.ExpandedHit hi : hits_bbs) {
                Integer hit_frag = hi.hit_fragments.keySet().iterator().next();
                SynthonSpace.FragId hit_frag_id = hi.hit_fragments.get(hit_frag);

                SynthonSpace.FragType hit_ft = new SynthonSpace.FragType(hi.rxn, hit_frag_id.frag);
                if( ! sorted_bb_hits.containsKey(hit_ft) ) { sorted_bb_hits.put(hit_ft,new ArrayList<>()); }
                sorted_bb_hits.get(hit_ft).add(hi);
            }

            for (SynthonSpace.FragType fti_a : sorted_bb_hits.keySet()) {

                List<SynthonSpace.FragId> all_matching_fragids = new ArrayList<>();
                for(SynthonSpace.ExpandedHit hi : sorted_bb_hits.get(fti_a)) {
                    all_matching_fragids.add(hi.hit_fragments.get(fti_a.frag));
                }

                //SynthonSpace.CombinatorialHit chi = new SynthonSpace.CombinatorialHit();

                Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> frag_map = new HashMap<>();
                Map<Integer, Pair<SynthonSpace.FragType, BitSet>> fp_mapping = new HashMap<>();

                StereoMolecule sri_fragments[] = new StereoMolecule[space.getFragTypes(fti_a.rxn_id).size()];
                StereoMolecule mi_a = new StereoMolecule(mi); mi_a.ensureHelperArrays(Molecule.cHelperCIP);
                for(int szi=0;szi<sri_fragments.length;szi++) {
                    if(szi==0){ sri_fragments[szi] = mi_a; }
                    else{sri_fragments[szi] = new StereoMolecule();}
                }
                SynthonShredder.SplitResult sri = new SynthonShredder.SplitResult(sri_fragments,new ArrayList<>(),"",new ArrayList<>());

                // loop over frag types and fill the required maps:
                int other_splits_cnt = 1;
                for( Integer ftii : space.getFragTypes(fti_a.rxn_id).keySet() ) {
                    SynthonSpace.FragType fti = space.getFragTypes(fti_a.rxn_id).get(ftii);
                    if(ftii==fti_a.frag) {
                        fp_mapping.put(0,Pair.of(fti_a, fp_molecule));
                        frag_map.put(fti_a,all_matching_fragids);
                    }
                    else {
                        // put all fragments into combinatorial hit:
                        frag_map.put( fti , new ArrayList<>( space.getSynthonSet(fti_a.rxn_id,ftii) ) );
                        // put empty bitset into fp_mapping
                        fp_mapping.put( other_splits_cnt , Pair.of( fti , new BitSet( space.getBits() ) ) );
                        other_splits_cnt++;
                    }
                }
                SynthonSpace.CombinatorialHit hit_bb_ci = new SynthonSpace.CombinatorialHit(fti_a.rxn_id, frag_map, sri, fp_mapping);
                //hits_bbs_as_combinatorial_hits.add(hit_bb_ci);
                receiver.addCombinatorialHits(Collections.singletonList(hit_bb_ci),Collections.singletonList(null));
            }
        }
    }

    public static void expandIncompleteHits(SynthonSpace space, CachedDescriptorProvider cdp, List<SynthonSpace.CombinatorialHit> incomplete_hits) {
        //List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();
        //Set<String> discovered_rxns = new HashSet<>();

        if(!incomplete_hits.isEmpty()) {
            //BitSet fp_molecule = cdp.getFP_cached(mi);

            for (SynthonSpace.CombinatorialHit hit : incomplete_hits) {

                // 1. check if hit is incomplete
                List<SynthonSpace.FragType> fts = new ArrayList<>( space.getFragTypes( hit.rxn ).values());
                fts.removeAll( hit.hit_fragments.keySet() );
                if(!fts.isEmpty()) {
                    // OK we have to expand..
                    for(SynthonSpace.FragType fti : fts) {
                        hit.hit_fragments.put( fti , space.getSynthonSet(hit.rxn,fti.frag) );
                    }
                }
            }
        }
    }


    /**
     * Same as the other function, but expands bridged bonds and runs the queries one after the other.
     *
     * @param space
     * @param cdp
     * @param mi
     * @param threads
     * @param fillIncompleteMappings
     * @param omitRxnsWithHitsFromLowerSplitNumber
     * @param receiver
     */
    public static void run_substructure_search_streaming_01_withBridgedBondsExpansion(SynthonSpace space, CachedDescriptorProvider cdp, StereoMolecule mi, int threads, boolean fillIncompleteMappings, boolean omitRxnsWithHitsFromLowerSplitNumber, SynthonSpace.CombinatorialHitReceiver receiver) throws Exception {
        List<StereoMolecule> expandedQueries = expandBridgedSearches(mi);

        for(StereoMolecule qi : expandedQueries) {
            if(Thread.currentThread().isInterrupted()) {return;}
            System.out.println("[INFO] next bridge expansion query: "+qi.getIDCode());
            run_substructure_search_streaming_01(space,cdp,qi,threads,fillIncompleteMappings,omitRxnsWithHitsFromLowerSplitNumber,receiver);
        }
    }


    /**
     *
     * @param space
     * @param cdp
     * @param mi
     * @param threads
     * @param fillIncompleteMappings if true, then adds all synthons of the non-matched synthon sets into the
     *                               combinatorial hits.
     * @param omitRxnsWithHitsFromLowerSplitNumber if true, then if we find hits at a specific number of splits, we omit
     *                                             the synthon reaction in more complicated searches.
     * @return
     */
    public static void run_substructure_search_streaming_01(SynthonSpace space, CachedDescriptorProvider cdp, StereoMolecule mi, int threads, boolean fillIncompleteMappings, boolean omitRxnsWithHitsFromLowerSplitNumber, SynthonSpace.CombinatorialHitReceiver receiver) {

        // 1. bb hits:
        System.out.println("sss: bb screen");
        List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space,cdp,mi,threads);
        List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();

        Set<String> discovered_rxns = new HashSet<>();
        expandBuildingBlockExpandedHits(space, cdp, hits_bbs, mi, new SynthonSpace.CombinatorialHitReceiver() {
            @Override
            public synchronized void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                receiver.addCombinatorialHits(hi,stats);
                discovered_rxns.addAll(hi.stream().map(xi -> xi.rxn).collect(Collectors.toList()));
            }
        });

        // 2.
        System.out.println("sss: 1split -> start");
        space.findExpandedHits_withConnProximityMatching_streaming(space, cdp, mi, 1, 2,
                (omitRxnsWithHitsFromLowerSplitNumber ? discovered_rxns : new HashSet<>()), 1000, threads, new SynthonSpace.CombinatorialHitReceiver() {
                    @Override
                    public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                        if(fillIncompleteMappings) {
                            expandIncompleteHits(space, cdp, hi);
                        }
                        discovered_rxns.addAll(hi.stream().map(xi -> xi.rxn).collect(Collectors.toList()));
                        receiver.addCombinatorialHits(hi,stats);
                    }
                });
        System.out.println("sss: 1split -> done");

        System.out.println("sss: 2split -> start");
        if(Thread.currentThread().isInterrupted()) {
            return;
        }

        space.findExpandedHits_withConnProximityMatching_streaming(space, cdp, mi, 2, 3,
                (omitRxnsWithHitsFromLowerSplitNumber ? discovered_rxns : new HashSet<>()), 1000, threads, new SynthonSpace.CombinatorialHitReceiver() {
                    @Override
                    public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                        if(fillIncompleteMappings) {
                            expandIncompleteHits(space, cdp, hi);
                        }
                        discovered_rxns.addAll(hi.stream().map(xi -> xi.rxn).collect(Collectors.toList()));
                        receiver.addCombinatorialHits(hi,stats);
                    }
                });
        System.out.println("sss: 2split -> done");
//        space.findExpandedHits_withConnProximityMatching_streaming(space, cdp, mi, 2, 3,
//                (omitRxnsWithHitsFromLowerSplitNumber ? discovered_rxns : new HashSet<>()), 1000, threads, new SynthonSpace.CombinatorialHitReceiver() {
//                    @Override
//                    public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
//                        discovered_rxns.addAll(hi.stream().map(xi -> xi.rxn).collect(Collectors.toList()));
//                        receiver.addCombinatorialHits(hi,stats);
//                    }
//                });
//        System.out.println("sss: 2split -> done");
        System.out.println("sss: 3split -> start");
        if(Thread.currentThread().isInterrupted()) {return;}
        space.findExpandedHits_withConnProximityMatching_streaming(space, cdp, mi, 3, 4,
                (omitRxnsWithHitsFromLowerSplitNumber ? discovered_rxns : new HashSet<>()), 1000, threads, new SynthonSpace.CombinatorialHitReceiver() {
                    @Override
                    public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                        discovered_rxns.addAll(hi.stream().map(xi -> xi.rxn).collect(Collectors.toList()));
                        receiver.addCombinatorialHits(hi,stats);
                    }
                });
        System.out.println("sss: 3split -> stop");
    }

    /**
     * For up to three bonds that represent bridges, this function enumerates
     * all structures with explicit paths in between.
     *
     * @return
     */
    public static List<StereoMolecule> expandBridgedSearches(StereoMolecule query) throws Exception {
        // 1. search all bridge bonds: , l:bond,m:min:r:max
        List<Triple<Integer,Integer,Integer>> bridge_bonds = new ArrayList<>();
        query.ensureHelperArrays(Molecule.cHelperCIP);

        for(int zi=0;zi<query.getBonds();zi++) {
            //int bqf = query.getBondQueryFeatures(zi);
            if( query.isBondBridge(zi) ) {
                bridge_bonds.add(Triple.of(zi,query.getBondBridgeMinSize(zi),query.getBondBridgeMaxSize(zi)) );
            }
        }

        if(bridge_bonds.size()==0) { // just return the query molecule..
            List<StereoMolecule> mols = new ArrayList<>();
            mols.add(query);
            return mols;
        }

        if(bridge_bonds.size()>2) {throw new Exception("Too many bridge bonds, max supported: 2");}

        List<StereoMolecule> expandedQueries = new ArrayList<>();

        if(bridge_bonds.size()==1) {
            for (int za = bridge_bonds.get(0).getMiddle(); za <= bridge_bonds.get(0).getRight(); za++) {
                StereoMolecule expanded = expandBridgeBond(query,bridge_bonds.get(0).getLeft(),za);
                expandedQueries.add(expanded);
            }
        }
        if(bridge_bonds.size()==2) {
            for (int za = bridge_bonds.get(0).getMiddle(); za <= bridge_bonds.get(0).getRight(); za++) {
                for (int zb = bridge_bonds.get(1).getMiddle(); zb <= bridge_bonds.get(1).getRight(); zb++) {
                    StereoMolecule expanded_a = expandBridgeBond(query,bridge_bonds.get(0).getLeft(),za);
                    StereoMolecule expanded   = expandBridgeBond(expanded_a,bridge_bonds.get(1).getLeft(),zb);
                    expandedQueries.add(expanded);
                }
            }
        }

        return expandedQueries;
    }

    public static StereoMolecule expandBridgeBond(StereoMolecule mi_1, int bond, int num_bridge_atoms) {
        if(num_bridge_atoms<0) {return null;}
        int a1 = mi_1.getBondAtom(0,bond);
        int a2 = mi_1.getBondAtom(1,bond);


        StereoMolecule mi = new StereoMolecule(mi_1);
        mi.ensureHelperArrays(Molecule.cHelperCIP);

        // special handling for bridge bonds of size 0
        if(num_bridge_atoms==0) {
            List<Pair<Integer,Integer>> neighbors_a =new ArrayList<>();
            for(int ni =0;ni< mi.getConnAtoms(a1);ni++) {
                int nai = mi.getConnAtom(a1,ni);
                if(nai==a2){continue;} // skip the bridge neighbor..
                int nbondi = mi.getBondTypeSimple(mi.getBond(nai,a1));
                neighbors_a.add( Pair.of( nai , nbondi ));
            }

            // remove a, then connect b to all neighbors of a
            mi.markAtomForDeletion(a1);
            int amap[] = mi.deleteMarkedAtomsAndBonds();

            // insert bonds to b:
            for(Pair<Integer,Integer> pi : neighbors_a) {
                mi.addBond(amap[a2],amap[pi.getLeft()],pi.getRight());
            }
            mi.ensureHelperArrays(Molecule.cHelperCIP);
            return mi;
        }


        int qf_any_bond_type = Molecule.cBondTypeSingle | Molecule.cBondTypeDouble |
                Molecule.cBondTypeTriple | Molecule.cBondTypeDelocalized;

        mi.deleteBond(bond);
        int last = a1;
        for(int zi=0;zi<num_bridge_atoms;zi++) {
            int na = mi.addAtom(6);
            mi.setAtomQueryFeature(na,Molecule.cAtomQFAny,true);
            int nb = mi.addBond(last,na);
            mi.setBondQueryFeature(nb,qf_any_bond_type,true);
            last = na;
        }
        int nb = mi.addBond(last,a2);
        mi.setBondQueryFeature(nb,qf_any_bond_type,true);
        return mi;
    }

    public static void exportHitsToTSV_enumerated(List<SynthonSpace.CombinatorialHit> hits, BufferedWriter out, StereoMolecule reference, ProgressController pc) throws IOException {
        exportHitsToCSV_enumerated(hits,out,reference,pc,"\t");
    }

    public static void exportHitsToCSV_enumerated(List<SynthonSpace.CombinatorialHit> hits, BufferedWriter out, StereoMolecule reference, ProgressController pc, String delimiter) throws IOException {

        List<String> header = new ArrayList<>();
        header.add("Structure[idcode]");
        header.add("SynthonRxn");header.add("FragA[idcode]");header.add("FragB[idcode]");header.add("FragC[idcode]");;header.add("FragD[idcode]");
        header.add("FragA_id");header.add("FragB_id");header.add("FragC_id");header.add("FragD_id");
        out.write(String.join(delimiter,header)+"\n");
        //out.write("Structure[idcode],SynthonRxn,FragA[idcode],FragB[idcode],FragC[idcode],FragD[idcode],FragA_id,FragB_id,FragC_id,FragD_id\n");

        long hits_total = hits.stream().mapToLong(
                hi -> hi.hit_fragments.values().stream().mapToLong(vi -> vi.size() ).reduce((x,y) -> x*y ).getAsLong() ).sum();

        long hits_processed_cnt = 0;
        pc.startProgress("Exporting", 0,(int) hits_total );

        for(int zi=0;zi<hits.size();zi++) {
            SynthonSpace.CombinatorialHit hi = hits.get(zi);
            List<Pair<List<SynthonSpace.FragId>,StereoMolecule>> hits_enumerated = enumerateHits2( hi.hit_fragments );

            for(int zj=0;zj<hits_enumerated.size();zj++) {

                Pair<List<SynthonSpace.FragId>,StereoMolecule> ehi = hits_enumerated.get(zj);
                List<String> parts = new ArrayList<>();
                parts.add(ehi.getRight().getIDCode());
                parts.add(ehi.getLeft().get(0).rxn_id);
                List<String> frag_idcodes = new ArrayList<>();
                List<String> frag_ids     = new ArrayList<>();
                for(int zk=0;zk<4;zk++){
                    if(ehi.getLeft().size()<=zk) {
                        frag_idcodes.add("");
                        frag_ids.add("");
                    }
                    else {
                        frag_idcodes.add(ehi.getLeft().get(zk).idcode);
                        frag_ids.add(ehi.getLeft().get(zk).fragment_id);
                    }
                }
                parts.addAll(frag_idcodes);
                parts.addAll(frag_ids);

                out.write( String.join(delimiter,parts) + "\n" );
                hits_processed_cnt++;
                pc.updateProgress((int) hits_processed_cnt );
            }
        }
    }

    public static long countHits(SynthonSpace.CombinatorialHit hi) {
        return hi.hit_fragments.values().stream().mapToLong( vi -> vi.size() ).reduce( (x,y) -> x*y ).getAsLong();
    }

    public static long countHits(List<SynthonSpace.CombinatorialHit> hilist) {
        long hits = 0;
        for(int zi=0;zi<hilist.size();zi++) {
            SynthonSpace.CombinatorialHit hi = hilist.get(zi);
            hits += hi.hit_fragments.values().stream().mapToLong(vi -> vi.size()).reduce((x, y) -> x * y).getAsLong();
        }
        return hits;
    }


    /**
     *
     *
     * @param fragments
     *
     * @return
     */
    public static List<Pair<List<SynthonSpace.FragId>,StereoMolecule>> enumerateHits2(Map<SynthonSpace.FragType,List<SynthonSpace.FragId>> fragments) {

        Map<StereoMolecule, SynthonSpace.FragId> map_mol_to_fragid = new HashMap<>();

        IDCodeParser icp = new IDCodeParser();

        // create all assembly plans:
        List<List<StereoMolecule>> parts_to_assemble = new ArrayList<>();
        for(SynthonSpace.FragType fti : fragments.keySet().stream().sorted((x,y) -> Integer.compare(x.frag,y.frag)).collect(Collectors.toList())) {
            List<StereoMolecule> mlist = new ArrayList<>();
            for(SynthonSpace.FragId fid : fragments.get(fti)) {
                StereoMolecule si = new StereoMolecule();
                icp.parse(si,fid.idcode);
                mlist.add(si);
                map_mol_to_fragid.put(si,fid);
            }
            parts_to_assemble.add(mlist);
        }

        List<List<StereoMolecule>> partial_assemblies = new ArrayList<>();
        partial_assemblies.add(new ArrayList<>());

        for(List<StereoMolecule> pni : parts_to_assemble) {
            List<List<StereoMolecule>> next_partial_assemblies = new ArrayList<>();
            for(List<StereoMolecule> pai : partial_assemblies) {
                for(StereoMolecule pni_i : pni) {
                    List<StereoMolecule> pai_extended = new ArrayList<>(pai);
                    pai_extended.add(pni_i);
                    next_partial_assemblies.add(pai_extended);
                }
            }
            partial_assemblies = next_partial_assemblies;
        }

        List<Pair<List<SynthonSpace.FragId>,StereoMolecule>> assembled = new ArrayList<>();

        // then assemble all assemblies..
        for(List<StereoMolecule> asi : partial_assemblies) {
            StereoMolecule assembly_i = SynthonAssembler.assembleSynthons_faster(asi);
            List<SynthonSpace.FragId> assembly_in_fragids = new ArrayList<>();
            for(StereoMolecule asi_i : asi) { assembly_in_fragids.add(map_mol_to_fragid.get(asi_i));}
            assembled.add(Pair.of(assembly_in_fragids,assembly_i));
        }

        return assembled;
    }

    public static List<StereoMolecule> enumerateHits(SynthonSpace.CombinatorialHit hi) {

//        boolean all_frags_available = hi.hit_fragments.values().stream().
//                mapToInt( vi -> vi.size() ).allMatch( xi -> xi > 0 );
//        if(! all_frags_available) {
//            return new ArrayList<>();
//        }

        List<SynthonSpace.FragType> fragtypes = new ArrayList<>(hi.hit_fragments.keySet());

        List<List<StereoMolecule>> fms = new ArrayList<>();
        IDCodeParser icp = new IDCodeParser();

        for(SynthonSpace.FragType fti : fragtypes) {
            List<StereoMolecule> li = new ArrayList();
            for(SynthonSpace.FragId fid : hi.hit_fragments.get(fti)) {
                StereoMolecule mi = new StereoMolecule();
                icp.parse(mi,fid.idcode);
                li.add(mi);
            }
            fms.add(li);
        }

        List<StereoMolecule> assembled = fms.stream().reduce(
                (x,y) -> combine_bbs_xy(x,y)).get();
        return assembled;
//
//        if(fragtypes.size()==2) {
//            List<SynthonSpace.FragId> fa = hi.hit_fragments.get(fragtypes.get(0));
//            List<SynthonSpace.FragId> fb = hi.hit_fragments.get(fragtypes.get(1));
//            for(int za=0;za<fa.size();za++) {
//                for(int zb=0;za<fb.size();zb++) {
//                    List<StereoMolecule> mi = new StereoMolecule
//                            SynthonAssembler.
//                }
//            }
//        }

    }

    private static List<StereoMolecule> combine_bbs_xy(List<StereoMolecule> x, List<StereoMolecule> y) {

        List<StereoMolecule> assembled = new ArrayList<>();

        for(StereoMolecule xi : x) {
            for(StereoMolecule yi : y) {
                List<StereoMolecule> xy = new ArrayList<>();
                xy.add(xi);
                xy.add(yi);
                //assembled.add( SynthonAssembler.assembleSynthons_faster(xy) );
                assembled.add( SynthonAssembler.assembleSynthons(xy) );
            }
        }

        return assembled;
    }

    public static List<StereoMolecule> enumerateHits(List<SynthonSpace.CombinatorialHit> hits) {
        List<StereoMolecule> enumerated = new ArrayList<>();
        for(SynthonSpace.CombinatorialHit hi : hits) {
            enumerated.addAll( enumerateHits(hi) );
        }
        return enumerated;
    }


}
