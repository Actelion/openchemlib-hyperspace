package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.StereoMolecule;

import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SynthonSpaceSimilarityExplorer {

    private SynthonSimilaritySpace space;
    private StereoMolecule            query;

    private CachedDescriptorProvider  cachedDescriptorProvider;


    private int conf_MaxInitialHitsPerSplit    = 128;
    private int conf_MaxInitialHitsScreened    = 1000000;
    private int conf_MaxInitialHitsConsidered  = 10000;
    private int[] conf_NumExtensionsPerCandidate = new int[] { 8 , 5 }; // [0]: in case of one extension, [1]: in case of two extensions, per extension.
    private Set<String> conf_AllowedRxns         = null;


    /**
     *
     *
     * num_initially_screened_hits is the number of hits that we screen in the first phase.
     * This number can be chosen really high, as the process for finding candidate hits is really fast. After
     * finding candidate hits, a simple scoring function is used to measure the prospected similarity of the
     * candidate hits.
     * The best num_initial_hits hits then enter the hit expansion phase, and each will get num_extensions_total
     * extended structures.
     *
     * if you want to allow all rxns, then just pass "null" as c_AllowedRxns
     *
     **/
    public SynthonSpaceSimilarityExplorer(SynthonSimilaritySpace space, StereoMolecule query,
                                          int c_MaxInitialHitsPerSplit, int c_MaxInitialHitsScreened,
                                          int c_MaxInitialHitsConsidered, int c_NumExtensionsPerCandidate[] ,
                                          List<String> c_allowedRxns) {
        this.space = space;
        this.query = query;

        this.setConf_MaxInitialHitsPerSplit(c_MaxInitialHitsPerSplit);
        this.setConf_MaxInitialHitsScreened(c_MaxInitialHitsScreened);
        this.setConf_MaxInitialHitsConsidered(c_MaxInitialHitsConsidered);
        this.setConf_NumExtensionsPerCandidate(c_NumExtensionsPerCandidate);

        if(c_allowedRxns == null) {
            this.conf_AllowedRxns = new HashSet<>(space.fragment_map.keySet());
        }
        else {
            this.conf_AllowedRxns = new HashSet<>(c_allowedRxns);
        }
        //this.cachedDescriptorProvider = new CachedDescriptorProvider(space);
        this.cachedDescriptorProvider = new CachedDescriptorProvider(space.getDescriptorHandler().getInfo().shortName);
    }



    private static int logLevel_phase_1 = 1;
    private static int logLevel_phase_2 = 1;

    public static void setLogLevel(int p_logLevel_phase_1, int p_logLevel_phase_2) {
        logLevel_phase_1 = p_logLevel_phase_1;
        logLevel_phase_2 = p_logLevel_phase_2;
    }



    /**
     * @param r
     * @param mol
     * @param max_num_splits
     * @param max_fragments
     * @param thresh_initialHit_maxHD
     * @param out_sampled_similar
     */
    public void sampleSimilar_new_default(Random r, StereoMolecule mol, int max_num_splits, int max_fragments, Function<Integer,Integer> thresh_initialHit_maxHD , List<SimilarityHit> out_sampled_similar) {
        sampleSimilar_new_default(r, mol, max_num_splits, max_fragments, thresh_initialHit_maxHD, out_sampled_similar, 4);
    }
    /**
     *
     *
     * @param r
     * @param mol
     * @param max_num_splits
     * @param max_fragments
     * @param out_sampled_similar
     */
    public void sampleSimilar_new_default(Random r, StereoMolecule mol, int max_num_splits, int max_fragments, Function<Integer,Integer> thresh_initialHit_maxHD , List<SimilarityHit> out_sampled_similar, int num_threads ) {
        for(int zs = 0; zs<=max_num_splits;zs++) {
            System.out.println("[SIM_MAIN] sample with "+zs+" splits:");
            //sampleSimilar( r, mol, zs, max_fragments, thresh_initialHit_maxHD , thresh_expansionHits_maxHD, max_initial_hits, max_hits_per_extension, out_sampled_similar );
            sampleSimilar_new( r, mol, zs, max_fragments, thresh_initialHit_maxHD, out_sampled_similar , num_threads );
        }
    }



//    public void sampleSimilar( Random r, StereoMolecule mol, int num_splits, int max_fragments, int thresh_initialHit_maxHD , int thresh_expansionHits_maxHD, int max_hits, List<SimilarityHit> out_sampled_similar ) {
//        sampleSimilar(r,mol,num_splits,max_fragments,
//                ((Integer c) -> thresh_initialHit_maxHD) , ((Integer c) -> thresh_expansionHits_maxHD) ,
//                max_hits,out_sampled_similar);
//    }

    /**
     *
     * @param r
     * @param mol
     * @param num_splits
     * @param max_fragments
     * @param thresh_initialHit_maxHD
     * @param out_sampled_similar
     * @param num_threads
     */
    // removed from params: int num_initial_hits_per_split, int num_extension_hits_total , int max_hits,
    public void sampleSimilar_new(Random r, StereoMolecule mol, int num_splits, int max_fragments, Function<Integer,Integer> thresh_initialHit_maxHD , List<SimilarityHit> out_sampled_similar, int num_threads ) {

        long ts_start = System.currentTimeMillis();

        ExecutorService main_pool = Executors.newFixedThreadPool(num_threads);

        // PHASE 1 : compute initial hits
        int nb = mol.getBonds();

        List<int[]> splits = null;
        if (num_splits == 0) {
            splits = new ArrayList<>();
            splits.add(new int[0]);
        } else {
            splits = new ArrayList<>(CombinationGenerator.getAllOutOf(nb, num_splits).stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList()));
        }

        // shuffle (then we get rather different ones maybe?)
        Collections.shuffle(splits, r);


        final List<PrimaryHit> primary_hits = Collections.synchronizedList(new ArrayList<>());
        boolean found_enough_hits = false;

        if(logLevel_phase_1>0) {
            System.out.println("[SIM_P1] : start");
        }

        List<Future> tasks_phase_1 = new ArrayList<>();
        for (int si_int[] : splits) {
            int num_connectors = num_splits;

            //sampleSimilar_new_examineSplit(si_int,mol,max_fragments,num_connectors,thresh_initialHit_maxHD,num_initial_hits_per_split,primary_hits,r);
            Runnable r_simil = new Runnable() {
                @Override
                public void run() {
                    sampleSimilar_new_examineSplit(si_int,mol,max_fragments,num_connectors,thresh_initialHit_maxHD,primary_hits,r);
                }
            };
            tasks_phase_1.add(main_pool.submit(r_simil));
        }

        // wait for pool:
        int cnt_p1 = 0;

        for(Future ft : tasks_phase_1) {
            try {
                ft.get();
                cnt_p1++;
                if(logLevel_phase_1>0) {
                    System.out.print(".");
                    if(cnt_p1%60==0){System.out.println();}
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        long ts_p1_done = System.currentTimeMillis();

        if(logLevel_phase_1>-1) {
            System.out.println("\n[SIM_P1] done! Primary Hits: "+primary_hits.size() +
                    " Time P1= "+( ts_p1_done - ts_start ) );
        }

        // maybe output some statistics about the primary hits?

        // here we should have different options for the sorting that we apply (also consider fragment size / complexity?)
        primary_hits.sort( (x,y) -> -Double.compare( x.primary_tanimoto_sim , y.primary_tanimoto_sim ) );

        // CUT OFF:
        final List<PrimaryHit> considered_primary_hits = primary_hits.subList(0, Math.min(primary_hits.size(),this.conf_MaxInitialHitsConsidered) );

        if(logLevel_phase_1>-1) {
            if (!primary_hits.isEmpty()) {
                System.out.println("[SIM_P1] best tanimoto_sim= " + primary_hits.get(0).primary_tanimoto_sim);
            }
        }
        // PHASE 2 : expand initial hits

        List<SimilarityHit> similarity_hits = Collections.synchronizedList( new ArrayList<>() );

        // sort by rxn:
        Map<String,List<PrimaryHit>> sorted_hits = new HashMap<>();

        //for(PrimaryHit hi : primary_hits ) {
        for(PrimaryHit hi : considered_primary_hits ) {
            if(!sorted_hits.containsKey(hi.primary_hit.rxn_id)) {
                sorted_hits.put(hi.primary_hit.rxn_id , new ArrayList<>());
            }
            sorted_hits.get(hi.primary_hit.rxn_id).add(hi);
        }

        // add tasks in inner loop..
        List<Future> tasks_p2 = new ArrayList<>();

        for( String rxn : sorted_hits.keySet() ) {
            //System.out.println("Expand Hits Rxn "+rxn + " -> " + sorted_hits.get(rxn).size() + " primary hits");

            for(int zi=0;zi<sorted_hits.get(rxn).size();zi++) {
                PrimaryHit hi = sorted_hits.get(rxn).get(zi);
                //sampleSimilar_new_expandHit(hi,num_extension_hits_total,similarity_hits);
                Runnable r_hxp = new Runnable() {
                    @Override
                    public void run() {
                        sampleSimilar_new_expandHit(hi,similarity_hits);
                    }
                };
                tasks_p2.add(main_pool.submit(r_hxp));
            }
        }

        // wait for pool:
        int cnt_p2 = 0;

        if(logLevel_phase_2>-1) {
            System.out.println("[SIM_P2] start! (with "+tasks_p2.size()+" hits)");
        }

        for(Future ft : tasks_p2) {
            try {
                ft.get();
                cnt_p2++;
                if(logLevel_phase_2>0) {
                    System.out.print(".");
                    if(cnt_p2%60==0){System.out.println();}
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        long ts_p2_done = System.currentTimeMillis();

        if(logLevel_phase_2>-1) {
            System.out.println("[SIM_P2] done! Expanded hits: "+similarity_hits.size() +
                    " Time P2= "+( ts_p2_done - ts_p1_done ) +
                    " Time Total = " + (ts_p2_done-ts_start) );
        }

        // required?
        main_pool.shutdown();

        out_sampled_similar.addAll(similarity_hits);
    }

    private void sampleSimilar_new_expandHit(PrimaryHit hi, List<SimilarityHit> similarity_hits) {
        String rxn = hi.primary_hit.rxn_id;

        if(logLevel_phase_2>1) {
            System.out.println("Expand Hit Rxn " + rxn + " : " + hi.primary_hit_template.toString());
        }

        //ArrayList<BitSet> sets = new ArrayList<>( hi.extension_fragment_types.values() );
        ArrayList<SynthonSpace.FragType> frags = new ArrayList<>( hi.extension_fragment_types.keySet() );
        ArrayList<BitSet> sets    = new ArrayList<>( frags.stream().map( fi -> hi.extension_fragment_types.get(fi) ).collect(Collectors.toList()) );


        int num_extensions = frags.size();
        //int frags_per_extension = (num_extensions==1)?num_extension_hits_total: (int) Math.ceil( Math.sqrt( num_extension_hits_total ) );
        int frags_per_extension = this.conf_NumExtensionsPerCandidate[num_extensions-1];

        // all frags for extension i will go into i'th position
        List<Set<SynthonSpace.FragId>> fragments_for_expansion = new ArrayList<>();

        for(int ze=0;ze<num_extensions;ze++) {
            SynthonSpace.FragType fti        = frags.get(ze);
            BitSet    bi         = hi.extension_fragment_types.get( fti );

            //Set<Integer>  conni  = space.getConnectorSetForFragType(fti);
            BitSet  conni  = space.getConnectorSetForFragType(fti);

            if(logLevel_phase_2>1) {
                System.out.println("Expand extension " + ze + " : " + bi);
            }

            List<BitSet> extension_hits = new ArrayList<>();
            //space.mANN_Rxns.get(rxn).get(fti.frag).findAllNearestNeighbors(bi,frags_per_extension,extension_hits);
            space.getReactionANNs().get(rxn).get(fti.frag).findAllNearestNeighbors(bi,frags_per_extension,extension_hits);

            // find the right fragment(s).. (multiple is extremely unlikely, but you'll never know..)
            Set<SynthonSpace.FragId> fragments_i = new HashSet<>();
            for(BitSet bsi : extension_hits) {
                fragments_i.addAll( space.fragment_map_2.get(fti.rxn_id).get(fti.frag).get(bsi) );
            }
            fragments_for_expansion.add(fragments_i);
        }


        //List<int[]> assemblies = new ArrayList<>();
        // now we manually check for 1 / 2 extensions and start assembling our hits :)
        Map<Integer,List<SynthonSpace.FragId>> assembly_plan = new HashMap<>();
        if (num_extensions == 1) {
            if(logLevel_phase_2>1){ System.out.println("Start pre-assembling " + fragments_for_expansion.get(0).size() + " molecules" );}
            //assembly_plan.put( hi.extension_fragment_types.get(sets.get(0)).frag , new ArrayList<>( fragments_for_expansion.get(0) ) );
            assembly_plan.put( frags.get(0).frag , new ArrayList<>( fragments_for_expansion.get(0) ) );
        }

        if (num_extensions == 2) {
            if(logLevel_phase_2>1){ System.out.println("Start pre-assembling " + fragments_for_expansion.get(0).size() + " x " + fragments_for_expansion.get(1).size() + " molecules" ); }
            //assembly_plan.put( hi.extension_fragment_types.get(sets.get(0)).frag , new ArrayList<>( fragments_for_expansion.get(0) ) );
            //assembly_plan.put( hi.extension_fragment_types.get(sets.get(1)).frag , new ArrayList<>( fragments_for_expansion.get(1) ) );
            assembly_plan.put( frags.get(0).frag , new ArrayList<>( fragments_for_expansion.get(0) ) );
            assembly_plan.put( frags.get(1).frag , new ArrayList<>( fragments_for_expansion.get(1) ) );
        }

        // Assemble :)
        SimilarityHit sim_hit = new SimilarityHit( space, hi.primary_hit , assembly_plan);
        similarity_hits.add(sim_hit);
    }

    private void sampleSimilar_new_examineSplit(int si_int[], StereoMolecule mol, int max_fragments, int num_connectors, Function<Integer,Integer> thresh_initialHit_maxHD, List<PrimaryHit> primary_hits, Random r) {

        boolean found_enough_hits = primary_hits.size() >= this.conf_MaxInitialHitsScreened;

        int num_initial_hits_per_split = this.conf_MaxInitialHitsPerSplit;

        if (found_enough_hits) {
            //break;
            return;
        }

        SynthonShredder.SplitResult si = SynthonShredder.trySplit(mol, si_int, max_fragments);

        if (si == null) {
            //continue;
            return;
        }

        List<SynthonShredder.SplitResult> split_vars = SynthonShredder.getAllSplitVariations(si, num_connectors);
        //List<SplitResult> split_vars = getAllSplitVariations(si, 4);

        for (SynthonShredder.SplitResult svi : split_vars) {

            found_enough_hits = primary_hits.size()>=this.conf_MaxInitialHitsScreened;

            if (found_enough_hits) {
                break;
            }

            // check for initial hit:
            List<BitSet> hits_initial = new ArrayList<>();
            //BitSet fp_main = space.getFP_cached(svi.fragments[0]);
            BitSet fp_main = cachedDescriptorProvider.getFP_cached(svi.fragments[0]);


            //Set<Integer> primary_connectorset     = SynthonSpace.computeConnectorSet(svi.fragments[0]);
            BitSet primary_connectorset     = SynthonSpace.computeConnectorBitSet(svi.fragments[0]);

            // is this enough?
            //int max_ann_hits = max_hits * 8;
            //mANN.findNearestNeighbors(r,fp_main,thresh_initialHit_maxHD.apply( fp_main.cardinality() ),max_ann_hits,hits_initial);

            //mANNs.get(primary_connectorset).findNearestNeighbors(r, fp_main, thresh_initialHit_maxHD.apply(fp_main.cardinality()), max_ann_hits, hits_initial);
            //mANNs.get(primary_connectorset).findAllNearestNeighbors(fp_main, num_initial_hits_per_split , hits_initial );
            //if(!mANNs.containsKey(primary_connectorset)) {
            if(!space.getANNs().containsKey(primary_connectorset)) {
                continue;
            }
            //space.getANNs().get(primary_connectorset).findNearestNeighbors(r, fp_main, thresh_initialHit_maxHD.apply(fp_main.cardinality()), num_initial_hits_per_split, hits_initial);
            //int num_initial_hits_per_split = 1024;
            //space.getANNs().get(primary_connectorset).findNearestNeighbors(r, fp_main, thresh_initialHit_maxHD.apply(fp_main.cardinality()),num_initial_hits_per_split , hits_initial);
            List<LSHProvider.Hit2> ann_hits = new ArrayList<>();
            space.getANNs().get(primary_connectorset).findAllNearestNeighbors(fp_main, thresh_initialHit_maxHD.apply(fp_main.cardinality()),num_initial_hits_per_split , ann_hits);

            if(true) {
                List<LSHProvider.Hit2> hits_exact = new ArrayList<>();
                space.getANNs().get(primary_connectorset).exactFindNearestNeighbors2(fp_main, thresh_initialHit_maxHD.apply(fp_main.cardinality()), 1000, hits_exact);
                System.out.println("exact hits: "+hits_exact.size() + " hits ann: " + ann_hits.size() );
                ann_hits = hits_exact;
            }
            if(ann_hits.size()>0) {
                hits_initial = ann_hits.stream().map(hi -> hi.bs).collect(Collectors.toList());
            }
            else {
                hits_initial = new ArrayList<>();
            }

            // now expand, then group by rxn, and check if rxn is compatible with cut set
            Map<String,List<SynthonSpace.FragId>> frags_by_rxn = new HashMap<>();

            for(BitSet bhit : hits_initial) {
                Set<SynthonSpace.FragId> frag_hits = space.fragments_by_connectors.get(primary_connectorset).get(bhit);
                for(SynthonSpace.FragId fid : frag_hits) {
                    if(!space.rxns_by_connector_config.containsKey(svi.connector_config)) {
                        //System.out.println("unknown connector config.. skip");
                        continue;
                    }
                    if(!conf_AllowedRxns.contains(fid.rxn_id)) {
                        continue;
                    }
                    if( space.rxns_by_connector_config.get(svi.connector_config).contains(fid.rxn_id) ) {
                        if(!frags_by_rxn.containsKey(fid.rxn_id)){ frags_by_rxn.put(fid.rxn_id,new ArrayList<>()); }
                        frags_by_rxn.get(fid.rxn_id).add(fid);
                    }
                }
            }

            // now for each rxn: figure out the mapping from cut fragments to fragType
            for( String rxn : frags_by_rxn.keySet() ) {
                for (SynthonSpace.FragId phit : frags_by_rxn.get(rxn)) {

                    SynthonSpace.FragType primary_fragtype    = new SynthonSpace.FragType(rxn,phit.frag);

                    // compute the extension options (by removing the primary connector type once):
                    //Map<Set<Integer>, List<SynthonSpace.FragType>> extension_options = new HashMap<>();// space.fragment_types_by_rxn_and_connector_config.get(rxn) );
                    //Map<Set<Integer>, List<SynthonSpace.FragType>> extension_options_template = space.fragment_types_by_rxn_and_connector_config.get(rxn);
                    Map<BitSet, List<SynthonSpace.FragType>> extension_options = new HashMap<>();// space.fragment_types_by_rxn_and_connector_config.get(rxn) );
                    Map<BitSet, List<SynthonSpace.FragType>> extension_options_template = space.fragment_types_by_rxn_and_connector_config.get(rxn);

                    // remove the phit fragtype:
                    //for (Set<Integer> ki : extension_options_template.keySet()) {
                    for (BitSet ki : extension_options_template.keySet()) {
                        ArrayList<SynthonSpace.FragType> frag_types_ki = new ArrayList<>( extension_options_template.get(ki) );
                        frag_types_ki.remove(primary_fragtype);
                        if(frag_types_ki.size()>0) {
                            extension_options.put(ki, frag_types_ki);
                        }
                    }


                    boolean unambig_mapping = extension_options.values().stream().mapToInt( li -> li.size() ).allMatch( siz -> siz==1 );

                    if (unambig_mapping) {
                        // âll unique connector configs, that's easy..

                        // BitSet should rather be FragType!!
                        //Map<BitSet, FragType> map_extensions = new HashMap<>();
                        //Map<BitSet, StereoMolecule> map_extension_templates = new HashMap<>();
                        Map<SynthonSpace.FragType, BitSet> map_extensions = new HashMap<>();
                        Map<SynthonSpace.FragType, StereoMolecule> map_extension_templates = new HashMap<>();

                        for (int zi = 1; zi < svi.fragments.length; zi++) {
                            //Set<Integer> connectors_fi = SynthonSpace.computeConnectorSet(svi.fragments[zi]);
                            BitSet connectors_fi = SynthonSpace.computeConnectorBitSet(svi.fragments[zi]);

                            //BitSet fp_fi = space.getFP_cached(svi.fragments[zi]);
                            BitSet fp_fi = cachedDescriptorProvider.getFP_cached(svi.fragments[zi]);

                            // find fragment type: (is unique, therefore just get the first)
                            //SynthonSpace.FragType fti = space.fragment_types_by_rxn_and_connector_config.get(rxn).get(connectors_fi).get(0); //.get(0);
                            SynthonSpace.FragType fti = extension_options.get(connectors_fi).get(0); //.get(0);

//                            map_extensions.put(fp_fi, fti);
//                            map_extension_templates.put(fp_fi, svi.fragments[zi]);
                            map_extensions.put(fti, fp_fi);
                            map_extension_templates.put(fti, svi.fragments[zi]);
                        }

                        double primary_tanimoto_sim = LSHProvider.tanimoto_similarity(fp_main,phit.fp);
                        primary_hits.add(new SynthonSpaceSimilarityExplorer.PrimaryHit(phit, svi.fragments[0],
                                map_extension_templates, map_extensions,primary_tanimoto_sim));
//                        for (SynthonSpace.FragId fid : frags_by_rxn.get(rxn)) {
//                            primary_hits.add(new FragmentSpaceSimilarityExplorer.PrimaryHit(fid, svi.fragments[0], map_extension_templates, map_extensions));
//                        }
                    }

                    else {
                        System.out.println("Hmm.. not yet implemented..");
                    }

                }
            }

//                Map<Set<Integer>,List<SynthonSpace.FragType>> frag_options = space.fragment_types_by_rxn_and_connector_config.get(rxn);
//
//                // NOW: check (aka verify..) the unambiguous mapping:
//                boolean unambig_mapping = false;
//
//                // hmm.. actually, we really have to consider it as ambiguous if there are two similar fragments
//                for( Set<Integer> ki : frag_options.keySet()) {
//                    //int s_initial = ( primary_connectorset.equals(ki) )?1:0;
//                    int s_initial = 0;
//                    if( frag_options.get(ki).size() - s_initial <= 1 ) {
//                        unambig_mapping = true;
//                        break;
//                    }
//                }
//
//                //if( frag_options.values().stream().mapToInt( li -> li.size() ).allMatch( siz -> siz==1 ) ) {
//                if(unambig_mapping) {
//                    // âll unique connector configs, that's easy..
//
//                    // BitSet should rather be FragType!!
//                    //Map<BitSet, FragType> map_extensions = new HashMap<>();
//                    //Map<BitSet, StereoMolecule> map_extension_templates = new HashMap<>();
//                    Map<SynthonSpace.FragType,BitSet> map_extensions = new HashMap<>();
//                    Map<SynthonSpace.FragType,StereoMolecule> map_extension_templates = new HashMap<>();
//
//                    for (int zi = 1; zi < svi.fragments.length; zi++) {
//                        Set<Integer> connectors_fi = SynthonSpace.computeConnectorSet(svi.fragments[zi]);
//
//                        //BitSet fp_fi = space.getFP_cached(svi.fragments[zi]);
//                        BitSet fp_fi = cachedDescriptorProvider.getFP_cached(svi.fragments[zi]);
//
//                        // find fragment type: (is unique, therefore just get the first)
//                        SynthonSpace.FragType fti = space.fragment_types_by_rxn_and_connector_config.get(rxn).get(connectors_fi).get(0); //.get(0);
//
////                            map_extensions.put(fp_fi, fti);
////                            map_extension_templates.put(fp_fi, svi.fragments[zi]);
//                        map_extensions.put(fti,fp_fi);
//                        map_extension_templates.put(fti,svi.fragments[zi]);
//                    }
//
//                    for (SynthonSpace.FragId fid : frags_by_rxn.get(rxn)) {
//                        primary_hits.add(new FragmentSpaceSimilarityExplorer.PrimaryHit(fid, svi.fragments[0], map_extension_templates, map_extensions));
//                    }
//                }
//                else {
//                    // shared connector configs. Here it is a bit more complicated..
//                    // actually, it is not that complicated: after setting one as primary, everything should be well defined..
//                    List<Set<Integer>> template_list = new ArrayList<>();
//
////                    for(int zx = 1; zx<svi.fragments.length ; zx++ ) {
////                        Set<Integer> cs_ext_a = SynthonSpace.computeConnectorSet(svi.fragments[zx]);
////                        template_list.add(cs_ext_a);
////                    }
//
//                    // loop over the fid that we pick for the initial hit
//                    for(SynthonSpace.FragId fid : frags_by_rxn.get(rxn)) {
//                        SynthonSpace.FragType fid_fragtype = new SynthonSpace.FragType(fid.rxn_id,fid.frag);
//
//                        // frag types (important, remove the primary hit frag from the possibles!)
//                        Map<Set<Integer>,List<SynthonSpace.FragType>> frag_options_i = new HashMap<>(frag_options);
//
//                        // create new HashMap without the fid..
//                        //frag_options_i.values().stream().forEach(li -> li.remove(fid)); //
//                        for(Set<Integer> sxi : frag_options_i.keySet()) {
//                            frag_options_i.get(sxi).remove( fid_fragtype );
//                        }
//
//                        // now, basically all should be clear, there cannot be ambiguous options..
//                        boolean now_unambiguous = false;
//                        now_unambiguous = frag_options_i.values().stream().mapToInt( foi -> foi.size() ).allMatch( pi -> pi==1 );
//                        if(now_unambiguous) {
//                            Map<SynthonSpace.FragType, BitSet> map_extensions = new HashMap<>();
//                            Map<SynthonSpace.FragType, StereoMolecule> map_extension_templates = new HashMap<>();
//                            for (int zx = 1; zx < svi.fragments.length ; zx++ ) {
//                                Set<Integer> connectors_fi = SynthonSpace.computeConnectorSet(svi.fragments[zx]);
//                                SynthonSpace.FragType fti = frag_options_i.get(connectors_fi).get(0); // only one element left, so we just grab it.
//                                StereoMolecule frag_i = svi.fragments[zx];
//                                //map_extensions.put(fti, space.getFP_cached(frag_i));
//                                map_extensions.put(fti, cachedDescriptorProvider.getFP_cached(frag_i));
//                                map_extension_templates.put(fti, frag_i);
//                            }
//                            primary_hits.add(new FragmentSpaceSimilarityExplorer.PrimaryHit(fid, svi.fragments[0], map_extension_templates, map_extensions));
//                        }
//                        else {
//                            System.out.println("How can we end up here?? This is not good :(");
//                        }
////                        List<List<SynthonSpace.FragType>> all_assignments = SynthonSimilaritySpace.computeAllFragTypeOrderings(frag_options_i, template_list);
////                        for (List<SynthonSpace.FragType> assi : all_assignments) {
////                            Map<SynthonSpace.FragType, BitSet> map_extensions = new HashMap<>();
////                            Map<SynthonSpace.FragType, StereoMolecule> map_extension_templates = new HashMap<>();
////
////                            for (int zx = 0; zx < assi.size(); zx++) {
////                                SynthonSpace.FragType fti = assi.get(zx);
////                                StereoMolecule frag_i = svi.fragments[zx + 1];
////                                //map_extensions.put(fti, space.getFP_cached(frag_i));
////                                map_extensions.put(fti, cachedDescriptorProvider.getFP_cached(frag_i));
////                                map_extension_templates.put(fti, frag_i);
////                            }
////
////                            // maybe like this?
////                            primary_hits.add(new FragmentSpaceSimilarityExplorer.PrimaryHit(fid, svi.fragments[0], map_extension_templates, map_extensions));
////                        }
//                    }
//                }


//            }
        }
    }

    public CachedDescriptorProvider getCachedDescriptorProvider() {
        return this.cachedDescriptorProvider;
    }

    /**
     * Settings:
     */
    public int getConf_MaxInitialHitsPerSplit() {
        return conf_MaxInitialHitsPerSplit;
    }

    public void setConf_MaxInitialHitsPerSplit(int conf_MaxInitialHitsPerSplit) {
        this.conf_MaxInitialHitsPerSplit = conf_MaxInitialHitsPerSplit;
    }

    public int getConf_MaxInitialHitsScreened() {
        return conf_MaxInitialHitsScreened;
    }

    public void setConf_MaxInitialHitsScreened(int conf_MaxInitialHitsScreened) {
        this.conf_MaxInitialHitsScreened = conf_MaxInitialHitsScreened;
    }

    public int getConf_MaxInitialHitsConsidered() {
        return conf_MaxInitialHitsConsidered;
    }

    public void setConf_MaxInitialHitsConsidered(int conf_MaxInitialHitsConsidered) {
        this.conf_MaxInitialHitsConsidered = conf_MaxInitialHitsConsidered;
    }

    public int[] getConf_NumExtensionsPerCandidate() {
        return conf_NumExtensionsPerCandidate;
    }

    public void setConf_NumExtensionsPerCandidate(int[] conf_NumExtensionsPerCandidate) {
        this.conf_NumExtensionsPerCandidate = conf_NumExtensionsPerCandidate;
    }

    /**
     * Contains all information about a primary hit, such that it is
     * ready for efficient expansion
     *
     */
    public static class PrimaryHit {
        public final SynthonSpace.FragId primary_hit;
        public final StereoMolecule primary_hit_template;
        //public final Map<BitSet,FragType> extension_fragment_types;
        //public final Map<BitSet,StereoMolecule> extension_fragment_templates;
        public final Map<SynthonSpace.FragType,BitSet> extension_fragment_types;
        public final Map<SynthonSpace.FragType,StereoMolecule> extension_fragment_templates;

        double primary_tanimoto_sim=0;

        public PrimaryHit(SynthonSpace.FragId primary_hit, StereoMolecule primary_hit_template, Map<SynthonSpace.FragType,StereoMolecule> extension_fragment_templates , Map<SynthonSpace.FragType,BitSet> extension_fragment_types, double primary_tanimoto_sim) {
            this.primary_hit = primary_hit;
            this.primary_hit_template = primary_hit_template;
            this.extension_fragment_types = extension_fragment_types;
            this.extension_fragment_templates = extension_fragment_templates;
            this.primary_tanimoto_sim = primary_tanimoto_sim;
        }

    }




    public static class SimilarityHit {
        public final SynthonSpace.FragId initial;
        public final Map<Integer,List<SynthonSpace.FragId>> extensions;
        public final SynthonSimilaritySpace space;

        public SimilarityHit(SynthonSimilaritySpace space, SynthonSpace.FragId initial, Map<Integer,List<SynthonSpace.FragId>> extensions) {
            this.space = space;
            this.initial = initial;
            this.extensions = extensions;
        }





        public List<FinalizedSimilarityHit> assembleAllStructures(CachedDescriptorProvider cached_dh, CachedStereoMoleculeProvider icp, StereoMolecule reference) {
            Set<Set<SynthonSpace.FragId>> already_assembled_sets = new HashSet<>();
            return assembleAllStructures(cached_dh,icp,already_assembled_sets,reference);
        }

        /**
         * Takes the reference to compute the hamming distance
         *
         * @param reference
         * @return
         */
        public List<FinalizedSimilarityHit> assembleAllStructures(CachedDescriptorProvider cached_dh, CachedStereoMoleculeProvider icp, Set<Set<SynthonSpace.FragId>> already_assembled_sets, StereoMolecule reference) {
            //BitSet fp_ref = space.getFP(reference);
            //BitSet fp_ref = this.space.getFP_cached(reference);
            BitSet fp_ref = cached_dh.getFP_cached(reference);

            Set<Set<SynthonSpace.FragId>> assembled_frag_id_sets   = already_assembled_sets; //new HashSet<>();

            if(logLevel_phase_2 > 1) {
                System.out.println(":assembleAllStructure: size: ( " + extensions.values().stream().map(vi -> "" + vi.size()).collect(Collectors.joining(",")) + ")");
            }

            List<FinalizedSimilarityHit> assembled    = new ArrayList<>();
            // hack, we assume there are at most two extension hits. (else we would do a recursion here..)
            List<Integer> extension_keys = new ArrayList<>(extensions.keySet());
            for( SynthonSpace.FragId fe_a : extensions.get(extension_keys.get(0)) ) {
                List<SynthonSpace.FragId> to_assemble = new ArrayList<>();
                to_assemble.add(initial);
                to_assemble.add(fe_a);
                if(extension_keys.size()==2) {
                    for (SynthonSpace.FragId fe_b : extensions.get(extension_keys.get(1))) {
                        List<SynthonSpace.FragId> to_assemble_2 = new ArrayList<>(to_assemble);
                        to_assemble_2.add(fe_b);
                        Set<SynthonSpace.FragId> to_assemble_2_set = new HashSet<>(to_assemble_2);
                        if(assembled_frag_id_sets.contains(to_assemble_2_set)) {
                            // already assembled this one..
                            continue;
                        }
                        else {
                            assembled_frag_id_sets.add(to_assemble_2_set);
                        }
                        StereoMolecule mi = SynthonSpaceExplorer.assembleBuildingBlocks(icp, to_assemble_2);
                        BitSet fp_mi = cached_dh.getFP_cached(mi);
                        int hamming_distance_i = LSHProvider.hamming(fp_ref,fp_mi);
                        Map<Integer, SynthonSpace.FragId> frag_assignments = new HashMap<>(); to_assemble_2.stream().forEach(fi -> frag_assignments.put( fi.frag , fi) );
                        FinalizedSimilarityHit fhi = new FinalizedSimilarityHit(this.space,initial.rxn_id,frag_assignments,new ArrayList<>(),mi.getIDCode(),hamming_distance_i);
                        assembled.add(fhi);
                    }
                }
                else if(extension_keys.size()==1){
                    Set<SynthonSpace.FragId> to_assemble_1_set = new HashSet<>(to_assemble);
                    if(assembled_frag_id_sets.contains(to_assemble_1_set)) {
                        // already assembled this one..
                        continue;
                    }
                    else {
                        assembled_frag_id_sets.add(to_assemble_1_set);
                    }
                    StereoMolecule mi = SynthonSpaceExplorer.assembleBuildingBlocks(icp, to_assemble);
                    BitSet fp_mi = cached_dh.getFP_cached(mi);
                    int hamming_distance_i = LSHProvider.hamming(fp_ref,fp_mi);
                    Map<Integer, SynthonSpace.FragId> frag_assignments = new HashMap<>(); to_assemble.stream().forEach(fi -> frag_assignments.put( fi.frag , fi) );
                    FinalizedSimilarityHit fhi = new FinalizedSimilarityHit(this.space,initial.rxn_id,frag_assignments,new ArrayList<>(),mi.getIDCode(),hamming_distance_i);

                    assembled.add(fhi);
                }
            }
            return assembled;
        }

        public int hashCode() {
            return this.extensions.hashCode() + this.initial.hashCode();
        }
        public boolean equals(Object o) {
            if(!(o instanceof SimilarityHit)){return false;}
            if(!((SimilarityHit) o).initial.equals(this.initial)){return false;}
            return ((SimilarityHit) o).extensions.equals(this.extensions);
        }


    }



    public static class FinalizedSimilarityHit {
        public final SynthonSpace space;
        public final String rxn;
        public final Map<Integer, SynthonSpace.FragId> matched_fragments;
        public final List<Integer> unmatched_fragments;

        public final String idcode_assembled;
        public final int hamming_distance;

        public static final String idcode_u;
        public static final String idcode_np;
        public static final String idcode_pu;
        public static final String idcode_am;

        static {
            StereoMolecule sm_u   = new StereoMolecule();
            StereoMolecule sm_np  = new StereoMolecule();
            StereoMolecule sm_pu  = new StereoMolecule();
            StereoMolecule sm_am  = new StereoMolecule();
            sm_u.addAtom(92);
            sm_np.addAtom(93);
            sm_pu.addAtom(94);
            sm_am.addAtom(95);

            idcode_u    = sm_u.getIDCode();
            idcode_np   = sm_np.getIDCode();
            idcode_pu   = sm_pu.getIDCode();
            idcode_am   = sm_am.getIDCode();
        }

        public FinalizedSimilarityHit(SynthonSpace space, String rxn, Map<Integer, SynthonSpace.FragId> matched_fragments, List<Integer> unmatched_fragments, String idcode_assembled, int hamming_distance) {
            this.space = space;
            this.rxn = rxn;
            this.matched_fragments = matched_fragments;
            this.unmatched_fragments = unmatched_fragments;
            this.idcode_assembled = idcode_assembled;
            this.hamming_distance = hamming_distance;
        }

        public long getNumberOfStructures() {
            long count = 1;
            for(Integer ki : unmatched_fragments) {
                count *= this.space.ffps_sorted_by_rxn_and_frag_BT.get(rxn).get(ki).root.countAll();
            }
            return count;
        }

        public String toString() {
            List<String> matched = new ArrayList<>();
            for(Integer ki : matched_fragments.keySet().stream().sorted().collect(Collectors.toList()) ) {
                matched.add( matched_fragments.get(ki).fragment_id );
            }
            for(Integer ki : unmatched_fragments.stream().sorted().collect(Collectors.toList()) ) {
                matched.add( "UNMATCHED["+ki+"]");
            }

            String str = this.rxn + " : "+String.join(" + ",matched) + "  ["+this.getNumberOfStructures()+"]";
            return str;
        }


        public String toCSV(int mol_format) {
            return toCSV(mol_format,",");
        }

        /**
         * mol_format: 1 -> idcode
         * mol_format: 2 -> smiles
         * @param mol_format
         * @return
         */
        public String toCSV(int mol_format, String delimiter) {
            List<String> csv = new ArrayList<>();

            String idc = this.idcode_assembled;
            csv.add( (mol_format==1)?idc: HyperspaceUtils.idcodeToSmiles(idc) );

            csv.add(this.rxn.replaceAll(delimiter,"_"));

            List<Integer> keys =  this.matched_fragments.keySet().stream().sorted().collect(Collectors.toList());

            for(int zi=0;zi<=2;zi++) {
                if( ! (zi<keys.size())  ){
                    if(mol_format==1){csv.add(idcode_pu);}else{csv.add("[Pu]");}
                }
                else {
                    int i = keys.get(zi);
                    String fidc = this.matched_fragments.get(i).idcode;
                    csv.add( (mol_format==1)?fidc: HyperspaceUtils.idcodeToSmiles(fidc) );
                }
//                if(this.matched_fragments.containsKey(i)) {
//                    String fidc = this.matched_fragments.get(i).idcode;
//                    csv.add( (mol_format==1)?fidc: Hyperspace_old.idcodeToSmiles(fidc) );
//                }
//                else if(this.unmatched_fragments.contains(i)) {
//                    //csv.add("unmatched");
//                    if(mol_format==1){csv.add(idcode_u);}else{csv.add("[U]");}
//                }
//                else {
//                    //csv.add("none");
//                    //csv.add("[Pu]");
//                    if(mol_format==1){csv.add(idcode_pu);}else{csv.add("[Pu]");}
//                }
            }

            //csv.add( ""+this.getNumberOfStructures() );

            //for(int i=1;i<=3;i++) {
            for(int zi=0;zi<=2;zi++) {
                if(zi<keys.size()) {
                    int i = keys.get(zi);
                    boolean is_in = this.matched_fragments.containsKey(i);
                    csv.add((is_in) ? this.matched_fragments.get(i).fragment_id.replaceAll(delimiter, "_") : "na");
                }
                else {
                    csv.add( "na" );
                }
            }

            csv.add(""+this.hamming_distance);

            return String.join( delimiter , csv );
        }

    }

}
