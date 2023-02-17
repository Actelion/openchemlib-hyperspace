package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


import java.util.*;
import java.util.stream.Collectors;

public class SynthonSimilaritySpace extends SynthonSpace {

    //private int logLevel = 3;

    //private int logLevel_findCandidates = 1;
    //private int logLevel_hitExpansion   = 2;
     //private static int logLevel_findCandidates = 2;
    private static int logLevel_findCandidates = 0;
    private static int logLevel_hitExpansion   = 3;


    //private final int BITS = 512;




    public static void setLogLevel_B(int p_logLevel_hitExpansion, int p_logLevel_findCandidates) {
        logLevel_hitExpansion    = p_logLevel_hitExpansion;
        logLevel_findCandidates  = p_logLevel_findCandidates;
    }




//    public interface SimilarityHitCallback {
//        public void addSimilarityHit(SimilarityHit hi);
//    }

    public SynthonSimilaritySpace(SynthonSpace ss_space) {
        //this.cache_FPs = ss_space.cache_FPs;
        this.BITS = ss_space.BITS;
        this.bsts_by_connectors = ss_space.bsts_by_connectors;
        this.fragment_type_map  = ss_space.fragment_type_map;
        this.fragment_map_2 = ss_space.fragment_map_2;
        this.ffps_sorted_by_rxn_and_frag_BT = ss_space.ffps_sorted_by_rxn_and_frag_BT;
        this.fragment_map = ss_space.fragment_map;
        this.fragments_by_connectors = ss_space.fragments_by_connectors;
        this.rxns_by_connector_config = ss_space.rxns_by_connector_config;
        this.fragment_types_by_rxn_and_connector_config = ss_space.fragment_types_by_rxn_and_connector_config;
        this.setFP(ss_space.getDescriptorHandler_PerThread(),ss_space.getBits());
    }




    public SynthonSimilaritySpace getThis() {return this;}

    //private LSHOracle mANN = null;
    //private Map<Set<Integer>, LSHProvider> mANNs = new HashMap<>();
    private Map<BitSet, LSHProvider> mANNs = new HashMap<>();

    public LSHProvider getANN(BitSet ki) {
        return mANNs.get(ki);
    }

    //public Map<Set<Integer>, LSHProvider> getANNs() { return mANNs; }
    public Map<BitSet, LSHProvider> getANNs() { return mANNs; }

    private Map<String,Map<Integer, LSHProvider>> mANN_Rxns = new HashMap<>();



    public void initApproximateNN_Default() {
        //for(Set<Integer> ksi : this.fragments_by_connectors.keySet()){
        for(BitSet ksi : this.fragments_by_connectors.keySet()){
            LSHProvider mANN = LSHProvider.initDefault( new ArrayList<>(this.fragments_by_connectors.get(ksi).keySet()) );
            mANN_Rxns = new HashMap<>();
            for( String rxn : this.fragment_map.keySet() ) {
                if(!mANN_Rxns.containsKey(rxn)){mANN_Rxns.put(rxn,new HashMap<>());}
                for(Integer ki : this.fragment_map.get(rxn).keySet() ) {
                    mANN_Rxns.get(rxn).put( ki , LSHProvider.initDefault( new ArrayList<>(this.fragment_map.get(rxn).get(ki).stream().map(fi -> fi.fp ).collect(Collectors.toList())) ) );
                }
            }
            mANNs.put(ksi,mANN);
        }
    }


    public void initApproximateNN(int num_lsh_functions, int projection_bits) {
        initApproximateNN(num_lsh_functions, projection_bits , 1);
    }

    public void initApproximateNN(int num_lsh_functions, int projection_bits, int threads) {
        //mANN = LSHOracle.initDefault(num_lsh_functions, this.BITS, projection_bits, new ArrayList<>(this.fragments.keySet()) );
        int num_total_synthons     = (int) this.fragments_by_connectors.values().stream().flatMap( fi -> fi.values().stream() ).mapToInt(si -> si.size()).sum();
        int num_processed_synthons = 0;

        int num_total_synthon_sets = (int) this.fragment_map.values().stream().mapToInt(vi -> vi.keySet().size()).sum();
        int processed_synthon_sets = 0;

        System.out.print("\nInit main fragment ANNs ("+this.fragments_by_connectors.keySet().size()+") : ");
        //for(Set<Integer> ksi : this.fragments_by_connectors.keySet()) {
        for(BitSet ksi : this.fragments_by_connectors.keySet()) {
            LSHProvider mANN = LSHProvider.initDefault(num_lsh_functions, this.BITS, projection_bits, new ArrayList<>(this.fragments_by_connectors.get(ksi).keySet()),threads);
            mANNs.put(ksi,mANN);
            System.gc();
            System.out.print(".");
        }

        System.out.println("Done\nInit rxn synthon ANNs ("+num_total_synthon_sets+")\n\n");

        mANN_Rxns = new HashMap<>();
        for (String rxn : this.fragment_map.keySet()) {
            if (!mANN_Rxns.containsKey(rxn)) {
                mANN_Rxns.put(rxn, new HashMap<>());
            }
            for (Integer ki : this.fragment_map.get(rxn).keySet()) {
                processed_synthon_sets++;
                num_processed_synthons += this.fragment_map.get(rxn).get(ki).size();
                System.out.println("[INFO] " + rxn + " : Start , Progress: " + processed_synthon_sets + "/"+num_total_synthon_sets);
                System.out.println( "[INFO] "+ rxn +" : " + " -> Init ANN -> " +rxn+" : "+ki + " #synthons="+this.fragment_map.get(rxn).get(ki).size());
                if (this.fragment_map.get(rxn).get(ki).size() <= 800) {
                    System.out.println("[INFO] "+rxn+" :  -->  number of fragments: " + this.fragment_map.get(rxn).get(ki).size() + " -> no ANN required\n");
                    mANN_Rxns.get(rxn).put(ki, LSHProvider.initDefault(1, this.getBits(), 0,
                            new ArrayList<>(this.fragment_map.get(rxn).get(ki).stream().map(fi -> fi.fp).collect(Collectors.toList())) , 1 ));
                } else {
                    mANN_Rxns.get(rxn).put(ki, LSHProvider.initDefault(num_lsh_functions, this.getBits(), projection_bits,
                            new ArrayList<>(this.fragment_map.get(rxn).get(ki).stream().map(fi -> fi.fp).collect(Collectors.toList())) , threads) );
                }
            }
        }

        System.out.println("\nDone with Init of ANN datastructures!\n");


    }

    public static void runNNTest(Random r, LSHProvider mANN, int n_tests, int max_dist) {
        mANN.test_nn(r,n_tests,max_dist);
    }


//    public void sampleSimilar_default( Random r, StereoMolecule mol, int num_splits, int max_fragments, int thresh_initialHit_maxHD , int thresh_expansionHits_maxHD, int max_hits, List<SimilarityHit> out_sampled_similar ) {
//        sampleSimilar_default(r,mol,num_splits,max_fragments, ((Integer c) -> thresh_initialHit_maxHD) , ((Integer c) -> thresh_expansionHits_maxHD) , max_hits,out_sampled_similar);
//    }
//
//    public void sampleSimilar_default( Random r, StereoMolecule mol, int num_splits, int max_fragments, Function<Integer,Integer> thresh_initialHit_maxHD , Function<Integer,Integer> thresh_expansionHits_maxHD, int max_hits, List<SimilarityHit> out_sampled_similar ) {
//        for(int zs = 0; zs<num_splits;zs++) {
//            System.out.println("sample with "+zs+" splits:");
//            //sampleSimilar( r, mol, zs, max_fragments, thresh_initialHit_maxHD , thresh_expansionHits_maxHD, max_initial_hits, max_hits_per_extension, out_sampled_similar );
//            sampleSimilar( r, mol, zs, max_fragments, thresh_initialHit_maxHD , thresh_expansionHits_maxHD, max_hits, out_sampled_similar );
//        }
//    }





    /**
     * returns all lists of FTs that have FTs that are compatible
     *
     * @param frag_type_conn_configs
     * @param template
     * @return
     */
    public static List<List<FragType>> computeAllFragTypeOrderings( Map<Set<Integer>,List<FragType>> frag_type_conn_configs , List<Set<Integer>> template ) {

        if(template.size()==1) {
            if(frag_type_conn_configs.get(template.get(0)).size()==1) {
                // actually, we will always end up in this case, the other one is not possible..
                ArrayList<FragType> ft = new ArrayList<>();
                ft.add(frag_type_conn_configs.get(template.get(0)).get(0));
                ArrayList<List<FragType>> ft_la = new ArrayList<>();
                ft_la.add(ft);
                return ft_la;
            }
            else { // !! size must be two! (it is not possible, that more than two synthons have the same conn config..)
                if (frag_type_conn_configs.get(template.get(0)).size() != 2) {
                    System.out.println("ERROR.. HOW IS THIS POSSIBLE!?! More than two synthons have the same connector config!?");
                } else {
                    ArrayList<FragType> fta = new ArrayList<>();
                    ArrayList<FragType> ftb = new ArrayList<>();
                    fta.add(frag_type_conn_configs.get(template.get(0)).get(0));
                    ftb.add(frag_type_conn_configs.get(template.get(0)).get(1));
                    ArrayList<List<FragType>> ft_la = new ArrayList<>();
                    ft_la.add(fta);
                    ft_la.add(ftb);
                    return ft_la;
                }
            }
        }


        List<List<FragType>> new_results = new ArrayList<>();

//        // compute the recursion..
//        List<Set<Integer>> template_2 = new ArrayList<>(template);
//        template_2.remove(0);
//        List<List<FragType>> result_i = computeAllFragTypeOrderings(frag_type_conn_configs,template_2);

        // loop over the possibilities for this step:
        Set<Integer> si = template.get(0);

        for( FragType ft : frag_type_conn_configs.get(si) ) {
            // 1. compute the adjusted remaiing available configs:
            Map<Set<Integer>,List<FragType>> rem_frag_type_conn_configs = new HashMap<>(frag_type_conn_configs);
            rem_frag_type_conn_configs.get(si).remove(ft);

            // 2. compute the new list
            List<Set<Integer>> template_2 = new ArrayList<>(template);
            template_2.remove(0);

            // compute recursion:
            List<List<FragType>> result_i = computeAllFragTypeOrderings(frag_type_conn_configs,template_2);

            for(List<FragType> lfi : result_i) {
                lfi.add(0,ft);
                new_results.add(lfi);
            }
        }

        return new_results;
    }


//    public void sampleSimilar(Random r, StereoMolecule mol, int num_splits, int max_fragments, Function<Integer,Integer> thresh_initialHit_maxHD, Function<Integer,Integer> thresh_expansionHits_maxHD , int max_hits, List<SimilarityHit> out_sampled_similar ) {
//
//
//        // !! Most likely we will find the same initial fragments as candidates over and over again.
//        //    So we want to process them only once of course !!
//        //Set<FragId> processed_initial_frags = new HashSet<>();
//        // at least we maintain a set, such that we do not add duplicates to the out_sampled.
//        Set<SimilarityHit> detected_hits = new HashSet<>();
//
//
//        int nb = mol.getBonds();
//
//        List<int[]> splits  = null;
//        if(num_splits==0) {
//            splits = new ArrayList<>();
//            splits.add(new int[0]);
//        }
//        else {
//            splits = new ArrayList<>( CombinationGenerator.getAllOutOf(nb, num_splits).stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList()) );
//        }
//
//        // shuffle (then we get rather different ones maybe?)
//        Collections.shuffle(splits,r);
//
//        boolean found_enough_hits = false;
//
//        int split_cnt = 0; // this is just for measuring progress
//        for(int si_int[]:splits) {
//            split_cnt++;
//            if(found_enough_hits){break;}
//
//            SplitResult si = trySplit(mol,si_int,max_fragments);
//
//            if(si==null) {
//                continue;
//            }
//
//            List<SplitResult> split_vars = getAllSplitVariations(si,num_splits);
//
//            for(SplitResult svi : split_vars) {
//                if(found_enough_hits){break;}
//
//                // check for initial hit:
//                List<BitSet> hits_initial = new ArrayList<>();
//                BitSet fp_main = getFP(svi.fragments[0]);
//                Set<Integer> primary_connectorset = computeConnectorSet(svi.fragments[0]);
//
//                // is this enough?
//                int max_ann_hits = max_hits*8;
//                //mANN.findNearestNeighbors(r,fp_main,thresh_initialHit_maxHD.apply( fp_main.cardinality() ),max_ann_hits,hits_initial);
//                mANNs.get(primary_connectorset).findNearestNeighbors(r,fp_main,thresh_initialHit_maxHD.apply( fp_main.cardinality() ),max_ann_hits,hits_initial);
//                if(logLevel_findCandidates>0 ){
//                    if(hits_initial.size()>0) {
//                        System.out.println("Found " + hits_initial.size() + " hits!");
//                        if (hits_initial.size() > 0) {
//                            System.out.println("Query structure vs. Hit structure:");
//                            for (BitSet bi : hits_initial) {
//                                //for (FragId fid : fragments.get(bi)) {
//                                for (FragId fid : fragments_by_connectors.get(primary_connectorset).get(bi)) {
//                                    System.out.println(Hyperspace_old.idcodeToSmiles(svi.fragments[0].getIDCode()) + "." + Hyperspace_old.idcodeToSmiles(fid.idcode));
//                                }
//                            }
//                            System.out.println("distances: ");
//                            for (BitSet bi : hits_initial) {
//                                System.out.println("" + LSHOracle.hamming(bi, fp_main));
//                            }
//                        }
//                    }
//                }
//
//                // count the connectors (we need this then when we check for actual library fragments)
//                int frag_cnt_u = 0; int frag_cnt_np = 0; int frag_cnt_pu = 0; int frag_cnt_am = 0;
//                for(int zi=0;zi<svi.fragments[0].getAtoms();zi++) {
//                    frag_cnt_u  += (svi.fragments[0].getAtomicNo(zi)==92)?1:0;
//                    frag_cnt_np += (svi.fragments[0].getAtomicNo(zi)==93)?1:0;
//                    frag_cnt_pu += (svi.fragments[0].getAtomicNo(zi)==94)?1:0;
//                    frag_cnt_am += (svi.fragments[0].getAtomicNo(zi)==95)?1:0;
//                }
//
//                if(hits_initial.size()>0) {
//                    for(BitSet bhit : hits_initial) {
//
//                        //Set<FragId> frag_hits = fragments.get(bhit);
//                        Set<FragId> frag_hits = fragments_by_connectors.get(primary_connectorset).get(bhit);
//
//                        for (FragId fid : frag_hits) {
////                            if( processed_initial_frags.contains(fid) ) {
////                                continue;
////                            }
////                            processed_initial_frags.add(fid);
//
//                            if(found_enough_hits) {break;}
//
//                            // check that reaction is compatible with number of fragments, i.e. the reaction must have at least
//                            // as many fragments as the current split
//                            if (this.ffps_sorted_by_rxn_and_frag_BT.get(fid.rxn_id).keySet().size() < svi.fragments.length) {
//                                // no hit..
//                                continue;
//                            }
//
//                            IDCodeParser icp = new IDCodeParser();
//                            StereoMolecule fmi = new StereoMolecule();
//                            icp.parse(fmi, fid.idcode);
//
//                            // this check is redundant now, but we leave it in for sanity checking..
//                            if(true) {
//                                // !! The primary fragment MUST have the same number of connectors !!
//                                int lf_cnt_u = 0;
//                                int lf_cnt_np = 0;
//                                int lf_cnt_pu = 0;
//                                int lf_cnt_am = 0;
//                                for (int zi = 0; zi < fmi.getAtoms(); zi++) {
//                                    lf_cnt_u += (fmi.getAtomicNo(zi) == 92) ? 1 : 0;
//                                    lf_cnt_np += (fmi.getAtomicNo(zi) == 93) ? 1 : 0;
//                                    lf_cnt_pu += (fmi.getAtomicNo(zi) == 94) ? 1 : 0;
//                                    lf_cnt_am += (fmi.getAtomicNo(zi) == 95) ? 1 : 0;
//                                }
//                                if (frag_cnt_u != lf_cnt_u || frag_cnt_np != lf_cnt_np || frag_cnt_pu != lf_cnt_pu || frag_cnt_am != lf_cnt_am) {
//                                    throw new Error("We cannot end up here. Catastrohpic error in sampleSimilar(..)! Wrong connectors in primary fragment..");
//                                    //continue;
//                                }
//                            }
//                            //FragType ft = new FragType(fid.rxn_id, fid.frag);
//
//                            // ok! now check for some matching expansions.
//                            // We do this recursively, always taking the next smaller fragment (i.e. as in the canonical
//                            // order of the svi.fragments
//                            List<Integer> remaining_fragments = new ArrayList<>();
//                            remaining_fragments.addAll( this.fragment_map.get(fid.rxn_id).keySet() );
//                            remaining_fragments.remove(new Integer(fid.frag) );
//
//                            List<BitSet> remaining_fps = new ArrayList<>();
//                            for( int zi=1;zi<svi.fragments.length;zi++) {
//                                BitSet fp_fi = getFP(si.fragments[zi]);
//                                remaining_fps.add(fp_fi);
//                            }
//
//                            //expand_fragment_fp(fid,remaining_fragments,remaining_fps);
//                            List<StereoMolecule> remaining_molecules = new ArrayList<>(Arrays.asList(svi.fragments));
//                            remaining_molecules.remove(svi.fragments[0]);
//
//                            List<Map<Integer,List<FragId>>> expansion_results = new ArrayList<>();
//                            // should we leave this as a parameter?
//                            int max_hits_per_extension = 64;
//                            expand_fragment_fp(r, fid, remaining_fragments, remaining_molecules, remaining_fps , thresh_expansionHits_maxHD, max_hits_per_extension, new HashMap<>() , expansion_results);
//
//
//                            if(!expansion_results.isEmpty()) {
//                                System.out.println("Found expanded initial hits with similar structures!");
//                                System.out.println("Size =" + expansion_results.size());
//                                for(Map<Integer,List<FragId>> exp_ri : expansion_results ) {
//                                    SimilarityHit si_new = new SimilarityHit(fid,exp_ri);
//                                    if(detected_hits.contains(si_new)){
//                                        continue;
//                                    }
//                                    else {
//                                        detected_hits.add(si_new);
//                                        out_sampled_similar.add(si_new);
//                                    }
//                                }
//
//                                if (out_sampled_similar.size() > max_hits) {
//                                    System.out.println("Done, found enough similar!");
//                                    found_enough_hits = true;
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

//    private void expand_fragment_fp(Random r, FragId main_fragment, List<Integer> frags_remaining, List<StereoMolecule> query_frags_remaining, List<BitSet> fps_remaining , int max_dist, int max_results, Map<Integer,List<FragId>> current_frag_assignments , List<Map<Integer,List<FragId>>> all_frag_assignments ){
//        expand_fragment_fp(r,main_fragment,frags_remaining,query_frags_remaining,fps_remaining,
//                ((Integer c)-> max_dist) , max_results,current_frag_assignments,all_frag_assignments );
//    }
//
//    private void expand_fragment_fp(Random r, FragId main_fragment, List<Integer> frags_remaining, List<StereoMolecule> query_frags_remaining, List<BitSet> fps_remaining , Function<Integer,Integer> max_dist, int max_results, Map<Integer,List<FragId>> current_frag_assignments , List<Map<Integer,List<FragId>>> all_frag_assignments ) {
//
//        // analyze the current fragments connectors:
//        boolean fe_cnt_u,fe_cnt_np,fe_cnt_pu,fe_cnt_am;
//        fe_cnt_u=fe_cnt_np=fe_cnt_pu=fe_cnt_am = false;
//        for (int za = 0; za < query_frags_remaining.get(0).getAtoms(); za++) {
//            if (query_frags_remaining.get(0).getAtomicNo(za) == 92) {fe_cnt_u=true;}
//            if (query_frags_remaining.get(0).getAtomicNo(za) == 93) {fe_cnt_np=true;}
//            if (query_frags_remaining.get(0).getAtomicNo(za) == 94) {fe_cnt_pu=true;}
//            if (query_frags_remaining.get(0).getAtomicNo(za) == 95) {fe_cnt_am=true;}
//        }
//
//
//        if(frags_remaining.size()==1 && fps_remaining.size()==1) {
//            // check this, then end recursion
//            List<FragId> results_exp = new ArrayList<>();
//            //this.mANN_Rxns.get(main_fragment.rxn_id).get( frags_remaining.get(0) ).findNearestNeighbors(r,fps_remaining.get(0),max_dist,max_results, results);
//            findMatchingExpansions(r, main_fragment.rxn_id, main_fragment, query_frags_remaining.get(0), fps_remaining.get(0), fe_cnt_u, fe_cnt_np, fe_cnt_pu, max_dist.apply(fps_remaining.get(0).cardinality()), max_results, results_exp);
//
//            if(!results_exp.isEmpty()) {
//                HashMap assignments_a = new HashMap<>(current_frag_assignments);
//                assignments_a.put(frags_remaining.get(0),results_exp);
//                all_frag_assignments.add(assignments_a);
//            }
//            return;
//        }
//        //for(Integer fi : frags_remaining) {
//        // take the next frag-fi, and search for FragIds that we can assign.
//        Integer fi = frags_remaining.get(0);
//        for (BitSet fpi : fps_remaining) {
//            List<FragId> results_exp = new ArrayList<>();
//            //this.mANN_Rxns.get(main_fragment.rxn_id).get( frags_remaining.get(0) ).findNearestNeighbors(r,fps_remaining.get(0),max_dist,max_results, results);
//            findMatchingExpansions(r, main_fragment.rxn_id, main_fragment, query_frags_remaining.get(0), fps_remaining.get(0), fe_cnt_u, fe_cnt_np, fe_cnt_pu, max_dist.apply(fps_remaining.get(0).cardinality()), max_results, results_exp);
//            if (!results_exp.isEmpty()) {
//                HashMap assignments_a = new HashMap<>(current_frag_assignments);
//                assignments_a.put(frags_remaining.get(0), results_exp);
//                List<Integer> frags_remaining_a = new ArrayList<>(frags_remaining.subList(1, frags_remaining.size()));
//                List<BitSet> fps_remaining_a = new ArrayList<>(fps_remaining.subList(1, fps_remaining.size()));
//                List<StereoMolecule> query_frags_remaining_a = new ArrayList<>(query_frags_remaining.subList(1, fps_remaining.size()));
//                expand_fragment_fp(r, main_fragment, frags_remaining_a, query_frags_remaining_a, fps_remaining_a, max_dist, max_results, assignments_a, all_frag_assignments);
//            } else {
//            }
//        }
//        //}
//
//    }
//
//
//
//
//    /**
//     *
//     * @param rxn
//     * @return
//     */
//    public boolean findMatchingExpansions(Random r, String rxn, FragId initial_frag, StereoMolecule template, BitSet template_fp, boolean has_blue, boolean has_red, boolean has_orange, int max_hamming_dist, int max_hits, List<FragId> out_frags) {
//
//        Set<Integer> connectorset = new HashSet<>();
//        if(has_blue){   connectorset.add(92); }
//        if(has_red){    connectorset.add(93); }
//        if(has_orange){ connectorset.add(94); }
//
//        // 1. find the right fragment type:
//        int found_type = -1;
//        for(Integer ki : this.fragment_map.get(rxn).keySet()) {
//            if(ki.equals(initial_frag.frag)) {continue;}
//            boolean ft_ok =  ( this.fragment_map.get(rxn).get(ki).get(0).has_blue == has_blue ) &&
//                    ( this.fragment_map.get(rxn).get(ki).get(0).has_red == has_red) &&
//                    ( this.fragment_map.get(rxn).get(ki).get(0).has_orange == has_orange);
//            if(ft_ok) {
//                found_type = ki;
//                break;
//            }
//        }
//
//        // 2. look for similar fragments :)
//        if(found_type >= 0) {
//            // grab list of fragments, then either use oracle or just normal screen to find matching ones.
//            // just normal:
//            List<BitSet> hits = new ArrayList<>();
//            //this.mANN_Rxns.get(rxn).get(found_type).exactFindNearestNeighbors(template_fp,max_hamming_dist,max_hits,hits);
//            this.mANN_Rxns.get(rxn).get(found_type).findNearestNeighbors(r,template_fp,max_hamming_dist,max_hits,hits);
//
//            if(hits.isEmpty()) {
//                if(logLevel_hitExpansion>0) {
//                    System.out.println("ExpansionHit::findMatchingxpansions : did not find expansion");
//                }
//                return false;
//            }
//
//            if(logLevel_hitExpansion>0) {
//                System.out.println("ExpansionHits! "+hits.size());
//            }
//
//            Map<BitSet,Set<FragId>> fragments_i = this.fragments_by_connectors.get(connectorset);
//            final int f_found_type = found_type;
//            for(BitSet bsi : hits) {
//                //out_frags.addAll( this.fragments.get(bsi).stream().filter( fi -> fi.rxn_id.equals(rxn) && fi.frag==f_found_type ).collect(Collectors.toList()) );
//                out_frags.addAll( fragments_i.get(bsi).stream().filter( fi -> fi.rxn_id.equals(rxn) && fi.frag==f_found_type ).collect(Collectors.toList()) );
//            }
//
//            if(logLevel_hitExpansion>1) {
//                String template_idcode = template.getIDCode();
//                String template_smiles = Hyperspace_old.idcodeToSmiles(template_idcode);
//                for(FragId fid : out_frags) {
//                    System.out.println(template_smiles+"."+ Hyperspace_old.idcodeToSmiles(fid.idcode) );
//                }
//            }
//
//            return true;
//        }
//        else {
//            return false;
//        }
//    }


    public Map<String,Map<Integer, LSHProvider>> getReactionANNs() {
        return this.mANN_Rxns;
    }


}
