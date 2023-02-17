package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.*;
import com.actelion.research.chem.descriptor.FingerPrintGenerator;

import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class SynthonSpaceExplorer {

    private static int       logLevel = 1;


    private static int       logLevel_finalize = 1;

    private boolean   logRejected = true;


    private SynthonSpace space;

    private StereoMolecule query;
    private BitSet         query_fp;

    private CachedDescriptorProvider  cachedDescriptorProvider;

    public static void setLogLevel(int p_logLevel, int p_logLevel_finalize, int p_logLevel_assemble ){
        logLevel= p_logLevel;
        logLevel_finalize = p_logLevel_finalize;
        SynthonAssembler.logLevel_assemble = p_logLevel_assemble;
    }


    public SynthonSpaceExplorer(SynthonSpace space, StereoMolecule query) {
        this.space = space;
        //this.query = query;
        this.setQuery(query);

        //this.cachedDescriptorProvider = new CachedDescriptorProvider(space);
        this.cachedDescriptorProvider = new CachedDescriptorProvider(space.getDescriptorHandler().getInfo().shortName);
    }

    public void setQuery(StereoMolecule query) {
        this.query = query;
        this.query.setFragment(true);
        this.query_fp = this.space.getFP(query);
    }

    private Map<Integer,List<SynthonSpace.ExpandedHit>> expanded_hit_candidates = new HashMap<>();

    private Map<Integer,Map<String,List<FinalizedHit>>> hits_finalized = new HashMap<>();

    public List<SynthonSpace.ExpandedHit> getExpandedHitCandidates(int splits) { return this.expanded_hit_candidates.get(splits); }



    public void screen_default(int max_splits, int max_fragments, int max_extension_hits_per_fragment) {
        screen_default(max_splits,max_fragments,max_extension_hits_per_fragment,4);
    }


    /**
     * The parameter max_extension_hits_per_fragment limits the number of extension hits that are generated before
     * we go into assembling and finalization. So in cases, where during the final check after finalization ( which
     * is mainly for stereo config things) assembled molecules are rejected, it may be necessary to increase this
     * number (but these cases should be very rare).
     *
     * @param max_splits
     * @param max_fragments
     * @param max_extension_hits_per_fragment
     */
    public void screen_default(int max_splits, int max_fragments, int max_extension_hits_per_fragment, int num_threads) {

        //int max_initial_hits = 256;
        int max_initial_hits = 4096;

        CachedStereoMoleculeProvider csmp = new CachedStereoMoleculeProvider();

        StereoMolecule sm = new StereoMolecule( this.query );
        sm.ensureHelperArrays(StereoMolecule.cHelperCIP);

        // here we collect the reactions with hits
        // reactions with hits will not be considered in subsequent scanning with more splits
        Set<String> discovered_rxns = new HashSet<>();

        // Screen Building Blocks
        System.out.println("[SS_MAIN] Screen building blocks");

        //List<SynthonSpace.ExpandedHit> hit_candidates_0 = space.screen_building_blocks(new StereoMolecule(sm));
        List<SynthonSpace.ExpandedHit> hit_candidates_0 = SynthonSpace.screen_building_blocks(space,cachedDescriptorProvider, new StereoMolecule(sm), num_threads);
        List<FinalizedHit> hits_bb_rejected         = new ArrayList<>();
        List<FinalizedHit> hits_bb                  = finalizeExpandedHits(csmp, hit_candidates_0, hits_bb_rejected, true, num_threads);
        discovered_rxns.addAll( hits_bb.stream().map( bbi -> bbi.rxn ).collect(Collectors.toSet()) );

        long num_structures = hits_bb.stream().mapToLong( bi -> bi.getNumberOfStructures() ).sum();
        storeResults(0,hits_bb);

        System.out.println("[SS_MAIN] BB-Hits in reactions (omit in next complexity level): "+discovered_rxns.stream().map(ri -> "["+ri+"]").reduce( (x,y) -> x+","+y ) );

        // Screen With Splits
        for(int splits = 1;splits<=max_splits;splits++) {

            System.out.println("[SS_MAIN] Screen with splits="+splits);
            //List<SynthonSpace.Hit> hit_candidates_1 = space.screen(sm, splits, discovered_rxns, max_fragments, max_initial_hits);
            List<SynthonSpace.InitialHit> hit_candidates_1 = space.screen( cachedDescriptorProvider,  sm, splits, discovered_rxns, max_fragments, max_initial_hits, num_threads);

            System.out.println("Candidate hits: "+hit_candidates_1.size());

            if(logLevel > 1) {
                List<String> hi_primary_frags = new ArrayList<>();
                for(SynthonSpace.InitialHit hi : hit_candidates_1) {
                    hi_primary_frags.add(hi.primary_fragment_id.idcode);
                }
                if(false) {
                    HyperspaceUtils.DebugOutput.plotMolecules("Hits", hi_primary_frags.stream().toArray(String[]::new), 8, 8);
                    try {
                        Thread.sleep(10000);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                }
            }

            //Map<SynthonSpace.Hit, List<SynthonSpace.ExpandedHit>> expanded_hit_candidates_sorted_1 = space.expandHits(hit_candidates_1, max_extension_hits_per_fragment);
            Map<SynthonSpace.InitialHit, List<SynthonSpace.ExpandedHit>> expanded_hit_candidates_sorted_1 = space.expandHits_New(hit_candidates_1, max_extension_hits_per_fragment, num_threads);
            List<SynthonSpace.ExpandedHit> expanded_hit_candidates_1 = expanded_hit_candidates_sorted_1.keySet().stream().flatMap(ki -> expanded_hit_candidates_sorted_1.get(ki).stream()).distinct().collect(Collectors.toList());


            if(logLevel>0) {
                if (!expanded_hit_candidates_1.isEmpty()) {
                    for (SynthonSpace.ExpandedHit ehi : expanded_hit_candidates_1) {
                        List<String> idcodes = new ArrayList<>();
                        idcodes.add(query.getIDCode());
                        ehi.hit_fragments.values().stream().map(fi -> fi.idcode).forEach(fi -> idcodes.add(fi));
                        if(false) {
                            HyperspaceUtils.DebugOutput.plotMolecules("X", idcodes.stream().toArray(String[]::new), 1, 4);
                        }
                    }
                }
            }

            this.expanded_hit_candidates.put(splits,expanded_hit_candidates_1);
            System.out.println("Candidate Hits from splits=" + splits + "  -> " + expanded_hit_candidates_1.size());


            if(logLevel>0) {
                System.out.println("Expanded hit candidates: "+expanded_hit_candidates_1.size());
                if(logLevel>1) {
                    for (SynthonSpace.ExpandedHit hi : expanded_hit_candidates_1) {
                        //System.out.println("Query:");
                        //System.out.println( idcodeToSmiles(query.getIDCode()) );
                        System.out.println("InitialHit: " + hi.hit_fragments.values().stream().map(fi -> SynthonAssembler.idcodeToSmiles(fi.idcode)).reduce((x, y) -> x + "." + y).orElse("ERROR").toString());
                        ;
                    }
                }
            }

            boolean include_expanded_molecules = true;//expanded_hit_candidates_1.size() <= 2000;

            List<FinalizedHit> hits_1 = new ArrayList<>();
            List<FinalizedHit> hits_1_rejected = Collections.synchronizedList(new ArrayList<>() );
            hits_1 = finalizeExpandedHits(csmp, expanded_hit_candidates_1, hits_1_rejected, include_expanded_molecules,num_threads);

            if(true) {
                List<String> rej_structures = new ArrayList<>();
                for(int zi=0;zi<Math.min(64,hits_1_rejected.size());zi++) {
                    rej_structures.add(hits_1_rejected.get(zi).idcode_assembled);
                }
                if(!rej_structures.isEmpty()) {
                    //DebugOutput.plotMolecules("rejected", rej_structures.stream().toArray(String[]::new), 8, 8);
                    System.out.println("----------------------- Rejected in finalization! Structures: --------------------------");
                    for(String rjs : rej_structures) {
                        System.out.println(rjs);
                    }
                    System.out.println("----------------------- Rejected in finalization! END-------- --------------------------");
                }
            }


            discovered_rxns.addAll(hits_1.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
            num_structures = hits_1.stream().map( hi -> hi.getNumberOfStructures()).reduce((x,y) -> x + y).orElse(0L);

            System.out.println("Confirmed Hits from splits="+ splits + "  -> " + hits_1.size() + " -> # structures = " + num_structures + " -> rejected in Finalization: "+hits_1_rejected.size());
            discovered_rxns.addAll(hits_1.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
            System.out.println("[SS_MAIN] Hit in reactions (omit in next complexity level): "+discovered_rxns.stream().map(ri -> "["+ri+"]").reduce( (x,y) -> x+","+y ) );
            storeResults(splits,hits_1);
        }


        // add all to the results
        //storeResults(1,hits_1);
        //storeResults(2,hits_2);
        //storeResults(3,hits_3);

    }



    private void storeResults(int splits, List<FinalizedHit> hits ) {
        Map<String,List<FinalizedHit>> sorted = new HashMap<>();
        for(FinalizedHit fhi : hits) {
            if(!sorted.containsKey(fhi.rxn)){sorted.put(fhi.rxn,new ArrayList<>());}
            sorted.get(fhi.rxn).add(fhi);
        }
        this.hits_finalized.put(splits,sorted);
    }

    public List<FinalizedHit> getAllResults() {
        List<FinalizedHit> all_hits = new ArrayList<>();
        for(Integer ki:this.hits_finalized.keySet()) {
            for(String rxn:this.hits_finalized.get(ki).keySet()) {
                all_hits.addAll(this.hits_finalized.get(ki).get(rxn));
            }
        }
        return all_hits;
    }

    private Map<String,long[]> hashed_indexes = new HashMap<>();

    public List<FinalizedHit> finalizeExpandedHits(CachedStereoMoleculeProvider csmp, List<SynthonSpace.ExpandedHit> candidates , List<FinalizedHit> out_rejected_hits , boolean include_assembled_molecule, int threads){

        //CachedStereoMoleculeProvider csmp = new CachedStereoMoleculeProvider();
        List<FinalizedHit> finalized_hits = Collections.synchronizedList(new ArrayList<>());

        if(logLevel_finalize>0) {
            System.out.println("[FIN] finalize hits:\n");
        }

        ExecutorService main_pool = Executors.newFixedThreadPool(threads);
        List<Future> tasks_finalize = new ArrayList<>();

        for(SynthonSpace.ExpandedHit hi : candidates) {
            //finalizeExpandedHits_finalizeExpandedHit(hi,finalized_hits,csmp,candidates,out_rejected_hits,include_assembled_molecule);
            Runnable r_fhi = new Runnable() {
                @Override
                public void run() {
                    finalizeExpandedHits_finalizeExpandedHit(hi,finalized_hits,csmp,candidates,out_rejected_hits,include_assembled_molecule);
                }
            };
            tasks_finalize.add( main_pool.submit(r_fhi) );
        }

        // wait for pool:
        int cnt_fin = 0;

        for(Future ft : tasks_finalize) {
            try {
                ft.get();
                cnt_fin++;
                if(logLevel_finalize>0) {
                    System.out.print(".");
                    if(cnt_fin%60==0){System.out.println();}
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
                System.out.println(e.getMessage());
            }
        }

        // needed?
        main_pool.shutdown();

        if(logLevel_finalize>0) {
            System.out.println("[FIN] finalize hits: Done!");
        }

        return finalized_hits;
    }

    private static String empty_idcode = (new StereoMolecule()).getIDCode();

    private void finalizeExpandedHits_finalizeExpandedHit(SynthonSpace.ExpandedHit hi, List<FinalizedHit> finalized_hits, CachedStereoMoleculeProvider csmp, List<SynthonSpace.ExpandedHit> candidates , List<FinalizedHit> out_rejected_hits , boolean include_assembled_molecule) {
        List<String> frag_idcodes = hi.hit_fragments.values().stream().map( fi -> fi.idcode ).collect(Collectors.toList());

        List<Integer> unmatched = new ArrayList<>(space.ffps_sorted_by_rxn_and_frag_BT.get(hi.rxn).keySet());
        unmatched.removeAll(hi.hit_fragments.keySet());

        if(include_assembled_molecule) {

            //StereoMolecule assembled = assembleBuildingBlocks(frag_idcodes);
            if(logLevel_finalize>1){
                System.out.println("FragIDs: "+ hi.hit_fragments.values().stream().map(vi->vi.fragment_id).reduce((x,y)->x + " , "+y) );
            }
            StereoMolecule assembled = SynthonSpaceExplorer.assembleBuildingBlocks( csmp, new ArrayList<>(hi.hit_fragments.values()) );

            if (logLevel_finalize >= 2) {
                System.out.println("query:     " + SynthonAssembler.idcodeToSmiles(query.getIDCode()));
                System.out.println("assembled: " + SynthonAssembler.idcodeToSmiles(assembled.getIDCode()));
                System.out.println("query:     " + query.getIDCode());
                System.out.println("assembled: " + assembled.getIDCode());
                System.out.println("Perform final substructure search");

                List<String> idcodes_to_plot = new ArrayList<>();
                idcodes_to_plot.add(query.getIDCode());
                idcodes_to_plot.add(assembled.getIDCode());
                //DebugOutput.plotMolecules( "LC:"+assembled.getIDCode() , idcodes_to_plot.stream().toArray(String[]::new),1,2 );
            }

            // !!! Else the sss wont do anything !!!
            StereoMolecule query_fragment = new StereoMolecule(query);
            query_fragment.setFragment(true);
            query_fragment.ensureHelperArrays(StereoMolecule.cHelperCIP);

            String idcode_query = query_fragment.getIDCode();
            String idcode_assembled = assembled.getIDCode();

            //SSSearcher searcher = new SSSearcher(SSSearcher.cDefaultMatchMode);
            SSSearcherWithIndex searcher = new SSSearcherWithIndex();

            if (!hashed_indexes.containsKey(idcode_query)) {
                hashed_indexes.put(idcode_query, searcher.createLongIndex(query));
            }
            if (!hashed_indexes.containsKey(idcode_assembled)) {
                hashed_indexes.put(idcode_assembled, searcher.createLongIndex(assembled));
            }

            long index_query[] = hashed_indexes.get(idcode_query);
            long index_assembled[] = hashed_indexes.get(idcode_assembled);

            //SmilesParser smp = new SmilesParser();
            //StereoMolecule que_a = new StereoMolecule();
            //StereoMolecule ass_a = new StereoMolecule();
            boolean validate_final_match = true;

            FinalizedHit fhit = new FinalizedHit(space, hi.rxn, hi.hit_fragments, unmatched, idcode_assembled);

            if (validate_final_match) {
                assembled.setFragment(false);
                assembled.ensureHelperArrays(StereoMolecule.cHelperCIP);

                //BitSet fp_assembled = space.getFP(assembled);
                BitSet fp_assembled = cachedDescriptorProvider.getFP_cached(assembled);

                //BitSet fp_assembled = space.getFP_cached(assembled);
                //BitSet b2 = (BitSet) fp_assembled.clone();
                //b2.or(query_fp);

                BitSet ti = (BitSet) fp_assembled.clone();
                ti.or(query_fp);
                if(!ti.equals(fp_assembled)) {
                    if(logLevel_finalize>1) {
                        System.out.println("Rejected by FP:");
                        System.out.println( HyperspaceUtils.idcodeToSmiles(query.getIDCode()) );
                        System.out.println( HyperspaceUtils.idcodeToSmiles(assembled.getIDCode()) );
                    }

//                    if(false) { // good if you want to test the sub-bla property of some fp. But for HyperFP it works.
//                        searcher.setMolecule(assembled, index_assembled);
//                        searcher.setFragment(query_fragment, index_query);
//                        boolean matched = searcher.isFragmentInMolecule();
//                        if(matched) {
//                            System.out.println("\n\n\n\n\n");
//                            System.out.println("------------------------------------");
//                            System.out.println("CATASTROPHIC ERROR: for");
//                            System.out.println( HyperspaceUtils.idcodeToSmiles(query.getIDCode()) );
//                            System.out.println( HyperspaceUtils.idcodeToSmiles(assembled.getIDCode()) );
//                            String[] fp_a = FragmentFingerPrintGenerator.getHyperFP(1024).getFingerPrintStrings(query);
//                            String[] fp_b = FragmentFingerPrintGenerator.getHyperFP(1024).getFingerPrintStrings(assembled);
//                            List<String> fps_a = new ArrayList<>(Arrays.asList(fp_a));
//                            fps_a.removeAll( Arrays.asList(fp_b) );
//                            System.out.println("------------------------------------");
//                            System.out.println("\n\n\n\n\n");
//                        }
//                    }

                    out_rejected_hits.add(fhit);
                }
                else {
                    if(logLevel_finalize>1) {
                        System.out.println("PASSED FP CHECK!");
                        System.out.println( HyperspaceUtils.idcodeToSmiles(query.getIDCode())+"."+ HyperspaceUtils.idcodeToSmiles(assembled.getIDCode()));
                        System.out.println( query.getIDCode()+" <-> "+ assembled.getIDCode() );
                    }
                    try {
                        searcher.setMolecule(assembled, index_assembled);
                        searcher.setFragment(query_fragment, index_query);

                        boolean matched = searcher.isFragmentInMolecule();


                        if (matched) {
                            if (logLevel_finalize >= 2) {
                                System.out.println("Verified Hit: "+ SynthonAssembler.idcodeToSmiles(assembled.getIDCode()));
                            }
                            finalized_hits.add(fhit);
                        } else {
                            if (logLevel_finalize >= 1 || logRejected) {
                                System.out.println("Hit rejected: "+ SynthonAssembler.idcodeToSmiles(query_fragment.getIDCode())+"."+ SynthonAssembler.idcodeToSmiles(assembled.getIDCode()));
                            }
                            out_rejected_hits.add(fhit);
                        }
                    } catch (Exception ex) {
                        System.out.println("error during final substructure matching");
                        ex.printStackTrace();
                    }
                }
            } else {
                finalized_hits.add(fhit);
            }
        }
        else {
            FinalizedHit fhit = new FinalizedHit(space, hi.rxn, hi.hit_fragments, unmatched, empty_idcode);
            finalized_hits.add(fhit);
        }
    }

    public static Map<String,StereoMolecule> mAssembledMoleculesCache = new HashMap<>();



    public static StereoMolecule assembleBuildingBlocks(List<SynthonSpace.FragId> frag_ids) {
        return assembleBuildingBlocks(new CachedStereoMoleculeProvider(),frag_ids);
    }

    /**
     * Assembles the building blocks at uranium atoms with similar color.
     *
     * NOTE! This will throw an error for wrongly annotated building block sets.
     * There seem to be mistakes in the building blocks m274860a_1.smi (not the actual version of Enamine REAL space btw.)
     *
     * @param frag_ids
     * @return
     */
    public static StereoMolecule assembleBuildingBlocks(CachedStereoMoleculeProvider icp , List<SynthonSpace.FragId> frag_ids) {


        // 1. check that we do not have more than two of every connector (this must never happen..)
        int connector_counts[] = new int[SynthonSpace.MAX_CONNECTORS];
        for(SynthonSpace.FragId fid : frag_ids) {
            fid.connectors.stream().forEach( zi -> connector_counts[zi]++ );
        }
        for(int zi=0;zi<SynthonSpace.MAX_CONNECTORS;zi++) {
            if(connector_counts[zi]>2) {
                throw new Error("Cannot assemble, more than 2 connectors found of type: "+zi);
            }
        }
//
//        int connectors_u = 0;
//        int connectors_np = 0;
//        int connectors_pu = 0;
//        int connectors_am = 0;
//        for (SynthonSpace.FragId fid : frag_ids) {
//            connectors_u += (fid.has_blue ? 1 : 0);
//            connectors_np += (fid.has_red ? 1 : 0);
//            connectors_pu += (fid.has_orange ? 1 : 0);
//            connectors_am += (fid.has_green ? 1 : 0);
//        }
//        if (connectors_u > 2) {
//            throw new Error("Cannot assemble, more than 2 U connectors found..");
//        }
//        if (connectors_np > 2) {
//            throw new Error("Cannot assemble, more than 2 Np connectors found..");
//        }
//        if (connectors_pu > 2) {
//            throw new Error("Cannot assemble, more than 2 Pu connectors found..");
//        }
//        if (connectors_am > 2) {
//            throw new Error("Cannot assemble, more than 2 Am connectors found..");
//        }

        //Collections.sort(parts);
        //String assembly_hash_key = parts.stream().reduce("",(x,y) -> x + "_x_" + y);
        frag_ids.sort((x, y) -> (x.fragment_id.compareTo(y.fragment_id)));

        //IDCodeParser p = new IDCodeParser();

        List<StereoMolecule> m_parts = new ArrayList<>();
        for (SynthonSpace.FragId fid : frag_ids) {
            if (SynthonAssembler.logLevel_assemble > 0) {
                System.out.println("Process: " + fid.fragment_id + "frag: ");
            }
            //StereoMolecule mi = new StereoMolecule();
            //p.parse(mi, fid.idcode);
            StereoMolecule mi = icp.parseIDCode(fid.idcode);
//            // color uranium atoms:
//            if(fid.uranium_pos_blue>=0) {
//                if( mi.getAtomicNo(fid.uranium_pos_blue)!=92 ) {
//                    System.out.println("Hyperspace input data corrupted.. :(");
//                }
//                mi.setAtomColor(fid.uranium_pos_blue,StereoMolecule.cAtomColorBlue);
//            }
//            if(fid.uranium_pos_red>=0) {
//                if( mi.getAtomicNo(fid.uranium_pos_red)!=92 ) {
//                    System.out.println("Hyperspace input data corrupted.. :(");
//                }
//                mi.setAtomColor(fid.uranium_pos_red,StereoMolecule.cAtomColorRed);
//            }
//            if(fid.uranium_pos_orange>=0) {
//                if( mi.getAtomicNo(fid.uranium_pos_orange)!=92 ) {
//                    System.out.println("Hyperspace input data corrupted.. :(");
//                }
//                mi.setAtomColor(fid.uranium_pos_orange,StereoMolecule.cAtomColorOrange);
//            }
            m_parts.add(mi);
        }


        //StereoMolecule assembled = assembleSynthons(m_parts);
        StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(m_parts);
        return assembled;
    }


    public CachedDescriptorProvider getCachedDescriptorProvider() {
        return cachedDescriptorProvider;
    }


    public static class FinalizedHit {
        public final SynthonSpace space;
        public final String rxn;
        public final Map<Integer, SynthonSpace.FragId> matched_fragments;
        public final List<Integer> unmatched_fragments;

        public final String idcode_assembled;

        public FinalizedHit(SynthonSpace space, String rxn, Map<Integer, SynthonSpace.FragId> matched_fragments, List<Integer> unmatched_fragments, String idcode_assembled) {
            this.space = space;
            this.rxn = rxn;
            this.matched_fragments = matched_fragments;
            this.unmatched_fragments = unmatched_fragments;
            this.idcode_assembled = idcode_assembled;
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


        /**
         * mol_format: 1 -> idcode
         * mol_format: 2 -> smiles
         * @param mol_format
         * @return comma separated: 1. assembled mol 2. synthon rxn 3 - 5: 
         *
         */
        public String toCSV(int mol_format) {
            return toCSV(mol_format,",");
        }

        /**
         * mol_format: 1 -> idcode
         * mol_format: 2 -> smiles
         * @param mol_format
         * @return comma separated: 1. assembled mol 2. synthon rxn 3 - 5:
         *
         */
        public String toCSV(int mol_format,String separator) {
            List<String> csv = new ArrayList<>();

            String idc = this.idcode_assembled;
            csv.add( (mol_format==1)?idc: SynthonAssembler.idcodeToSmiles(idc) );

            csv.add(this.rxn.replaceAll(separator,"_"));

            //for(int i=1;i<=3;i++) {
            // collect ALL keys, matched and unmatched..
            List<Integer> sorted_keys = new ArrayList<>(this.matched_fragments.keySet());
            sorted_keys.addAll(this.unmatched_fragments);
            Collections.sort(sorted_keys);

            for(int ii=1;ii<=3;ii++) {
                int i = (sorted_keys.size()>=ii)?sorted_keys.get(ii-1).intValue():1000000;
                if(this.matched_fragments.containsKey(i)) {
                    String fidc = this.matched_fragments.get(i).idcode;
                    csv.add( (mol_format==1)?fidc: SynthonAssembler.idcodeToSmiles(fidc) );
                }
                else if(this.unmatched_fragments.contains(i)) {
                    //csv.add("unmatched");
                    //csv.add("[U]");
                    if(mol_format==1) {
                        csv.add(SynthonSpaceSimilarityExplorer.FinalizedSimilarityHit.idcode_u);
                    }
                    else {
                        csv.add("[U]");
                    }
                }
                else {
                    //csv.add("none");
                    //csv.add("[Pu]");
                    if(mol_format==1) {
                        csv.add(SynthonSpaceSimilarityExplorer.FinalizedSimilarityHit.idcode_pu);
                    }
                    else {
                        csv.add("[Pu]");
                    }
                }
            }

            for(int ii=1;ii<=3;ii++) {
                int i = (sorted_keys.size() >= ii) ? sorted_keys.get(ii - 1).intValue() : 1000000;
                if (this.matched_fragments.containsKey(i)) {
                    String fidc = this.matched_fragments.get(i).fragment_id;
                    csv.add(fidc);
                } else if (this.unmatched_fragments.contains(i)) {
                    csv.add("unmatched");
                }
                else {
                    csv.add("na");
                }
            }

//            for(int i=1;i<=3;i++) {
//                boolean is_in =this.matched_fragments.containsKey(i);
//                csv.add( (is_in)?this.matched_fragments.get(i).fragment_id:"na");
//            }

            csv.add( ""+this.getNumberOfStructures() );

            //return String.join( "," , csv );
            return String.join( separator , csv );
        }

    }


    /**
     * Checks if B is subset of A
     *
     * @param a
     * @param b
     * @return
     */
    public static boolean testSubset(BitSet a, BitSet b) {
        BitSet a2 = (BitSet) a.clone();
        a2.or(b);
        return a2.equals(a);
    }

    public static void main(String args[]) {
        sss_test_01();
    }

    public static void sss_test_01() {

        String ssa = "CC(C=C(c(cc1)ccc1[N+]([O-])=O)[U])[U]";
        String ssb = "CC(/C=C(/c(cc1)ccc1[N+]([O-])=O)\\[U])[U]";


        SmilesParser sp = new SmilesParser();
        StereoMolecule sa = new StereoMolecule();
        StereoMolecule sb = new StereoMolecule();
        try {
            sp.parse(sa,ssa);
            sp.parse(sb,ssb);
        } catch (Exception e) {
            e.printStackTrace();
        }

        sa.setFragment(true);

        SSSearcher sss = new SSSearcher(SSSearcher.cIndexMatchMode);

        sss.setMolecule(sb);
        sss.setFragment(sa);

        FingerPrintGenerator fp = new FingerPrintGenerator();
        BitSet fp_a = fp.getFingerprint(sa);
        BitSet fp_b = fp.getFingerprint(sb);

        System.out.println("FP-Subset: "+testSubset(fp_b,fp_a));

        boolean match = sss.isFragmentInMolecule();
        System.out.println("Matching: "+match);
    }

}
