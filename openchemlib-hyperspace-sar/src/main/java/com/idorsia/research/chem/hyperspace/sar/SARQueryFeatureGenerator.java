package com.idorsia.research.chem.hyperspace.sar;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonShredder;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import org.apache.commons.io.input.CountingInputStream;
import org.apache.commons.lang3.tuple.Pair;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class SARQueryFeatureGenerator {

    /**
     *
     * Computes all single splits, takes the larger part of the molecule,
     * fixes it as non-extendable, and sets the connection point to the missing part
     * as requiresNeighbors.
     *
     * The queries are sorted from large to small (by size of the remaining part)
     *
     * @return
     */
    public static List<SingleSplitResult> generateSingleSplitFragments(StereoMolecule m) {
        m.ensureHelperArrays(Molecule.cHelperCIP);

        List<SingleSplitResult> computed_fragments = new ArrayList<>();

        for(int zi=0;zi<m.getBonds();zi++) {
            StereoMolecule mi = new StereoMolecule(m);
            mi.ensureHelperArrays(Molecule.cHelperCIP);

            SynthonShredder.SplitResult sri = SynthonShredder.trySplit(mi,new int[]{zi},2);
            if(sri!=null) {
                StereoMolecule qi = new StereoMolecule(sri.fragments[0]); // zero is the larger (i.e. largest) fragment
                qi.ensureHelperArrays(Molecule.cHelperCIP);
                qi.setFragment(true);
                // now configure all parts as fixed, except for the atom connected to the uranium atom:
                boolean delete_flags[] = new boolean[qi.getAtoms()];
                for(int za=0;za<qi.getAtoms();za++) {
                    if(qi.getAtomicNo(za)>=88) {
                        delete_flags[za]=true;
                        // and set require substitution query feature for neighbor atom
                        if( qi.getConnAtoms(za) != 1) {
                            // this is an error..
                            System.out.println("[ERROR] encountered unexpected scaffold, neighbor count of "+za+" = "+qi.getNonHydrogenNeighbourCount(za)+" but must be 1");
                            break;
                            //List<StereoMolecule> fragments_empty = new ArrayList<>();
                            //return fragments_empty;
                            //return new RepresentativeCombinatorialHits("","",null,null);
                        }
                        else {
                            // set requires substitution query feature
                            //int atom_a = qi.getConnAtom(za,0);
                            //fragment.setAtomQueryFeature(atom_a,StereoMolecule.cAtomQFMoreNeighbours,true);
                        }
                    }
                    else {
                        // check if atom has >=88 neighbor, else set NoMoreNeighbors query feature
                        boolean connected_to_connector = false;
                        for(int zj=0;zj<qi.getConnAtoms(za);zj++) {
                            if(qi.getAtomicNo( qi.getConnAtom(za,zj))>=88) {
                                connected_to_connector = true;
                            }
                        }

                        if(connected_to_connector) {
                            qi.setAtomQueryFeature(za,StereoMolecule.cAtomQFMoreNeighbours,true);
                        }
                        else {
                            qi.setAtomQueryFeature(za,StereoMolecule.cAtomQFNoMoreNeighbours,true);
                            //fragment.setAtomQueryFeature(za,StereoMolecule.cAtomQFNoMoreNeighbours,false);
                        }

                    }
                }
                qi.deleteAtoms(delete_flags);
                qi.ensureHelperArrays(Molecule.cHelperCIP);
                computed_fragments.add(new SingleSplitResult(qi,sri.fragments[0],sri.fragments[1]));
            }
        }
        return computed_fragments;
    }

    public static class SingleSplitResult {
        public final StereoMolecule split_query;

        // this is the remaining split with exit vector
        public final StereoMolecule split_A;

        // this is the cut away split with exit vector
        public final StereoMolecule split_B;

        public SingleSplitResult(StereoMolecule split_query, StereoMolecule split_A, StereoMolecule split_B) {
            this.split_query = split_query;
            this.split_A = split_A;
            this.split_B = split_B;
        }
    }

    public static void main_a(String args[]) {
        IDCodeParser icp = new IDCodeParser();


        StereoMolecule mi = new StereoMolecule();
        icp.parse(mi,"fc\u007FIb@LRBU@d`bDfU^ve}nUUISpV}WnfjZjjfj`A`ae@@");

        List<SingleSplitResult> fragments = generateSingleSplitFragments(mi);

        System.out.println("Fragments[idcode]");
        for(SingleSplitResult fi : fragments) {
            System.out.println(fi.split_query.getIDCode());
        }

        System.out.println("\n mkay..");
    }

    public static void main(String args[]) {
        SynthonSpace space = null;
        try {
            //space = loadSynthonSpace("C:\\dev4\\git4\\hyperspace_core\\real_space_latest_FragFp.data");
            space = loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
        } catch (IOException e) {
            e.printStackTrace();
        }

        //DWARFileParser dwar = new DWARFileParser( "C:\\datasets\\tl_full_real_space_scaffold_enumeration_with_cmavailable_check_2_NEW_b_20210428_CLC_Filters_Flagged_A_ntot_over_1k_with_085_clustering.dwar" );
        //DWARFileParser dwar = new DWARFileParser( "C:\\datasets\\tl_full_real_space_scaffold_enumeration_with_cmavailable_check_2_NEW_b_20210428_CLC_Filters_Flagged_A_ntot_over_400_TriazolesRemoved_clustered_centroids.dwar" );
        List<StereoMolecule> structures = new ArrayList<>();

        CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");

        if(true) {
            String test_a = "edZTLD@E@iJn`IABPrJYQQQKKZJIRiQISIFNYEGW@puU@pQUTu@AUDTGb@@"; // this one should work excellently.. [not an idorsia compound]
            String test_b = "fcwpS@NJzstU`dDRFICHiBidhheJeDeMiYsQuTuUTu@DSQBR@@"; // this one should work excellently.. [not an idorsia compound]

            IDCodeParser icp = new IDCodeParser();
            StereoMolecule ma = new StereoMolecule();
            StereoMolecule mb = new StereoMolecule();
            icp.parse(ma, test_a);
            icp.parse(mb, test_b);

            structures.add(ma);
            structures.add(mb);
        }

        for(StereoMolecule si : structures) {
            List<SingleSplitResult> fragments_i = generateSingleSplitFragments(si);

            System.out.println("Fragment[idcode]");
            for(SingleSplitResult ffi : fragments_i) {
                System.out.println(ffi.split_query.getIDCode());
            }

            System.out.println("Search fragment hits:");
            Map<String,List<SynthonSpace.CombinatorialHit>> hits = new HashMap<>();
            Map<String,Map<SynthonSpace.CombinatorialHit,List<StereoMolecule>>> hits_enumerated = new HashMap<>();
            for(SingleSplitResult ssri : fragments_i) {
                StereoMolecule fi = ssri.split_query;
                List<SynthonSpace.CombinatorialHit> hits_i = run_substructure_search_01(space,cdp,fi,3);
                hits.put(fi.getIDCode(),hits_i);
                System.out.println("Search fragment hits: found "+ hits_i.size() + " for fragment: " + fi.getIDCode());

                hits_enumerated.put(fi.getIDCode(),new HashMap<>());

                for(SynthonSpace.CombinatorialHit hi : hits_i) {
                    List<StereoMolecule> enum_hits_i = new ArrayList<>();
                    if(hi.hit_fragments.size()==2) {
                        Iterator<List<SynthonSpace.FragId>> frags = hi.hit_fragments.values().iterator();
                        List<SynthonSpace.FragId>[] extracted_frags = new List[] { frags.next() , frags.next() };
                        for(SynthonSpace.FragId fa : extracted_frags[0]) {
                            for(SynthonSpace.FragId fb : extracted_frags[1]) {
                                List<StereoMolecule> to_assemble = new ArrayList<>();
                                to_assemble.add( parseIDCode( fa.idcode ) );
                                to_assemble.add( parseIDCode( fb.idcode ) );
                                StereoMolecule assembled_i = SynthonAssembler.assembleSynthons_faster(to_assemble);
                                enum_hits_i.add(assembled_i);
                                if(enum_hits_i.size()>2000) {break;}
                            }
                        }
                    }
                    if(hi.hit_fragments.size()==3) {
                        Iterator<List<SynthonSpace.FragId>> frags = hi.hit_fragments.values().iterator();
                        List<SynthonSpace.FragId>[] extracted_frags = new List[] { frags.next() , frags.next() , frags.next() };
                        for(SynthonSpace.FragId fa : extracted_frags[0]) {
                            for(SynthonSpace.FragId fb : extracted_frags[1]) {
                                for(SynthonSpace.FragId fc : extracted_frags[2]) {
                                    List<StereoMolecule> to_assemble = new ArrayList<>();
                                    to_assemble.add( parseIDCode( fa.idcode ) );
                                    to_assemble.add( parseIDCode( fb.idcode ) );
                                    to_assemble.add( parseIDCode( fc.idcode ) );
                                    StereoMolecule assembled_i = SynthonAssembler.assembleSynthons_faster(to_assemble);
                                    enum_hits_i.add(assembled_i);
                                    if(enum_hits_i.size()>2000) {break;}
                                }
                            }
                        }
                    }
                    hits_enumerated.get(fi.getIDCode()).put(hi,enum_hits_i);
                }
            }

            try {
                BufferedWriter out = new BufferedWriter(new FileWriter("out_test_sar_" + (System.currentTimeMillis() % 1000) + ".csv"));
                //System.out.println("\n\nSearchStructure[idcode],rxn,Structure[idcode]");
                out.write("SearchStructure[idcode],rxn,Structure[idcode]\n");
                for (String ssi : hits_enumerated.keySet()) {
                    for (SynthonSpace.CombinatorialHit chi : hits_enumerated.get(ssi).keySet()) {
                        for (StereoMolecule emi : hits_enumerated.get(ssi).get(chi)) {
                            //System.out.println(ssi + "," + chi.rxn + "," + emi.getIDCode());
                            out.write(ssi + "," + chi.rxn + "," + emi.getIDCode()+"\n");
                        }
                    }
                }
                out.flush();
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }

    private static StereoMolecule parseIDCode(SynthonSpace.FragId fragId) {
        return parseIDCode(fragId.idcode);
    }
    private static StereoMolecule parseIDCode(String idcode) {
        IDCodeParser icp = new IDCodeParser();
        StereoMolecule mi = new StereoMolecule();
        icp.parse(mi,idcode);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        return mi;
    }

    public static List<SynthonSpace.CombinatorialHit> run_substructure_search_01(SynthonSpace space, CachedDescriptorProvider cdp, StereoMolecule mi, int threads) {

        List<SynthonSpace.CombinatorialHit> hits_no_splits = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_1s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_2s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_3s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_4s = new ArrayList<>();

        List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space,cdp,mi,threads);
        List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();

        Set<String> discovered_rxns = new HashSet<>();

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
                hits_bbs_as_combinatorial_hits.add(hit_bb_ci);
            }
        }

        discovered_rxns.addAll( hits_bbs_as_combinatorial_hits.stream().map( hi -> hi.rxn ).collect(Collectors.toList()) );

        hits_1s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,1,2,discovered_rxns, 200,1);
        discovered_rxns.addAll( hits_1s.stream().map( hi -> hi.rxn ).collect(Collectors.toList()) );
        hits_2s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,2,3,discovered_rxns, 200,1);
        discovered_rxns.addAll( hits_2s.stream().map( hi -> hi.rxn ).collect(Collectors.toList()) );
        hits_3s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,3,4,discovered_rxns,200,1);
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

        return all_hits;
    }






    public static SynthonSpace loadSynthonSpace(String file) throws IOException {
        //try {
        File fi_file = new File(file);
        long file_size_in_bytes = fi_file.length();

        FileInputStream file_in = new FileInputStream(fi_file);
        CountingInputStream counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
        ObjectInputStream in = new ObjectInputStream(counting_in_stream);
        SynthonSpace space = null;

        Thread t_progress_monitor = new Thread() {
            @Override
            public void run() {
                while (true) {
                    try {
                        Thread.sleep(5000);
                    } catch (InterruptedException e) {
                        System.out.println("monitor done!");
                        break;
                    }
                    long bytes_read = counting_in_stream.getByteCount();
                    double ratio = (1.0 * bytes_read) / file_size_in_bytes;
                    System.out.println("read: " + String.format("%4.2f", (ratio * 100)) + " i.e. " + (bytes_read / 1024.0) + " kb / " + (file_size_in_bytes / 1024.0) + " kb");
                }
            }
        };
        t_progress_monitor.start();

        try {
            space = (SynthonSpace) in.readObject();
        } catch (OptionalDataException ex) {
            ex.printStackTrace();
        }
        catch (ClassNotFoundException ex) {
            ex.printStackTrace();
        }

        t_progress_monitor.interrupt();
        space.initAfterJavaDeserialization();

        return space;
    }



}
