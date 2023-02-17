package com.idorsia.research.chem.hyperspace;

import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerSkeletonSpheres;
import com.actelion.research.gui.JEditableStructureView;
import org.apache.commons.lang3.tuple.Pair;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;


/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace
 */

public class SynthonSimilaritySpaceExplorer3 {


    private static int LOG_LEVEL_A = 2;


    public static void main(String args[]) {

        JFrame f = new JFrame();
        JEditableStructureView jsv = new JEditableStructureView();
        JButton jsearch = new JButton("Search");
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(jsv,BorderLayout.CENTER);
        f.getContentPane().add(jsearch,BorderLayout.NORTH);
        f.setSize(400,400);
        f.setVisible(true);


        SynthonSimilaritySpace3 sa = createSpace_A(false);

        CachedDescriptorProvider cdp    = new CachedDescriptorProvider("FragFp");
        SimilaritySearchConfig3 conf = new SimilaritySearchConfig3(10,3,3,0.85,8000);

        jsearch.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                StereoMolecule mi = new StereoMolecule(jsv.getMolecule());
                mi.ensureHelperArrays(Molecule.cHelperCIP);

                Thread ti = new Thread() {
                    @Override
                    public void run() {
                        SynthonSimilaritySpaceExplorer3 explorer = new SynthonSimilaritySpaceExplorer3(sa,cdp);
                        List<Pair<FinalSimilarityResult,Double>> results = explorer.runSearch(mi,conf,(x)->{});

                        BitSet mi_descriptor_ffp = BitSet.valueOf((new DescriptorHandlerLongFFP512()).createDescriptor(mi));
                        Random r = new Random();
                        String filename_i = ("sim_results_"+ (Math.abs(r.nextInt())%10000)+".csv");
                        writeResultsToCSV(filename_i,mi_descriptor_ffp,results);
                    }
                };
                ti.start();
            }
        });


        if(false) {

            StereoMolecule mi = HyperspaceUtils.parseIDCode("");

            System.out.println("Mkay, ready to go");


            SynthonSimilaritySpaceExplorer3 explorer = new SynthonSimilaritySpaceExplorer3(sa, cdp);
            explorer.runSearch(mi, conf, (x)->{});

            List<Pair<MoleculeAssembly, Double>> hits_a = findHits(mi, 1, sa, cdp, conf);
            BitSet mi_descriptor = cdp.getFP_cached(mi);
            byte[] mi_dh_skelspheres = (new DescriptorHandlerSkeletonSpheres()).createDescriptor(mi);

            System.out.println("mkay, ready to assemble");
            //List<Pair<StereoMolecule,Pair<MoleculeAssembly,Double>>> hits_a_expanded = hits_a.parallelStream().map(
            //        hi -> Pair.of( hi.getLeft().getAssembled() , hi)        ).collect(Collectors.toList());
            List<FinalSimilarityResult> hits_a_expanded = hits_a.parallelStream().map(
                    hi -> new FinalSimilarityResult(hi.getLeft(), hi.getRight(), hi.getLeft().getAssembled())).collect(Collectors.toList());

            System.out.println("mkay, everything assembled");
            System.out.println("mkay, compute descriptors A");
            hits_a_expanded.parallelStream().forEach(pi -> pi.computeFFP());
            System.out.println("mkay, compute descriptors B");
            hits_a_expanded.parallelStream().forEach(pi -> pi.computeSkelSpheres());

            System.out.println("mkay, computed descriptors");
            List<Pair<FinalSimilarityResult, Double>> hits_a_expanded_scored = hits_a_expanded.stream().map(
                    hi -> Pair.of(hi, (double) (new DescriptorHandlerSkeletonSpheres()).getSimilarity(mi_dh_skelspheres, hi.dh_skelspheres))).collect(Collectors.toList());
            hits_a_expanded_scored.sort((x, y) -> -Double.compare(x.getRight(), y.getRight()));

            // output:
            System.out.println("Results:\n\n");
            Random r = new Random();
            BufferedWriter data_out = null;
            try {
                data_out = new BufferedWriter(new FileWriter(("sim_results_" + (Math.abs(r.nextInt()) % 10000) + ".csv")));
                //System.out.println("Hit[idcode],Sim,Rxn,SimA2");
                data_out.write("Hit[idcode],SimSkelSpheres,SimFragFp,Rxn,SimA2" + "\n");
                for (Pair<FinalSimilarityResult, Double> asi : hits_a_expanded_scored) {
                    List<String> result_parts = new ArrayList<>();
                    result_parts.add(asi.getLeft().mol.getIDCode());
                    result_parts.add("" + asi.getRight());
                    result_parts.add("" + LSHProvider.tanimoto_similarity(mi_descriptor, asi.getLeft().dh_ffp));
                    result_parts.add(asi.getLeft().assembly.hit.initial_hit.frag.rxn_id);
                    result_parts.add("" + asi.getLeft().assembly.getScore_A(2));
                    //System.out.println(String.join(",",result_parts));
                    data_out.write(String.join(",", result_parts) + "\n");
                }
                data_out.flush();
                data_out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

            //SimilaritySearchConfig3 conf3 = new SimilaritySearchConfig3(6,3,3,0.8);

            List<Pair<MoleculeAssembly, Double>> hits_a_2 = findHits(mi, 1, sa, cdp, conf);
            List<Pair<MoleculeAssembly, Double>> hits_b = findHits(mi, 2, sa, cdp, conf);
            List<Pair<MoleculeAssembly, Double>> hits_c = findHits(mi, 3, sa, cdp, conf);
        }
    }

    public static SynthonSimilaritySpace3 createSpace_A(boolean small) {
        SynthonSpace sspace = null;
        try {
            if(!small) {
                sspace = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
            }
            else {
                sspace = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/real_space_latest_SMALL_FragFp.data");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        sspace.initAfterJavaDeserialization();
        SynthonSimilaritySpace3 space3 = new SynthonSimilaritySpace3(sspace);
        space3.initFastSimilaritySearchers(10,(x->{}));

        // try to serialize..
        try {
            ObjectOutputStream out_space = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream("space_test_out.data")));
            out_space.writeObject(space3);
            out_space.flush();
            out_space.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


        return space3;
    }


    public static class SimilaritySearchConfig3 {
        public final int    num_threads;

        public final int    max_splits;
        public final int    max_fragments;
        public final double initial_tanimoto;

        public final SynthonSimilaritySpace3.TopoConstraints topo_constraints;

        public final int max_hits_per_split_level;

        public SimilaritySearchConfig3(int threads,
                                       int max_splits, int max_fragments,
                                       double initial_tanimoto,
                                       int max_hits_per_split_level) {
            this.num_threads = threads;

            this.max_splits = max_splits;
            this.max_fragments = max_fragments;
            this.initial_tanimoto = initial_tanimoto;

            // load default
            this.topo_constraints = new SynthonSimilaritySpace3.TopoConstraints();
            this.max_hits_per_split_level = max_hits_per_split_level;
        }
    }



    private SynthonSimilaritySpace3   space3;
    private CachedDescriptorProvider  cdp;


    public SynthonSimilaritySpaceExplorer3(SynthonSimilaritySpace3  space, CachedDescriptorProvider cdp) {
        this.space3 = space;
        this.cdp    = cdp;
    }

    private static int LOG_LEVEL_SEARCH = 3;

    public List<Pair<FinalSimilarityResult,Double>> runSearch(StereoMolecule query, SimilaritySearchConfig3 config, Consumer<Double> progress) {

        if(LOG_LEVEL_SEARCH > 0) {System.out.println("[INFO] Query: "+query.getIDCode());}

        BitSet query_descriptor = cdp.getFP_cached(query);
        byte[] query_dh_skelspheres = (new DescriptorHandlerSkeletonSpheres()).createDescriptor(query);


        Map<Integer,List<Pair<FinalSimilarityResult,Double>>> map_all_hits = new HashMap<>();

        progress.accept(0d);

        for(int splits = 1; splits <= config.max_splits; splits++) {
            if(LOG_LEVEL_SEARCH > 0) {System.out.println("[INFO] Start Splits: "+splits);}

            List<Pair<MoleculeAssembly,Double>> hits_a   = findHits(query,splits,this.space3,cdp,config);

            if(LOG_LEVEL_SEARCH > 1) {
                System.out.println("[INFO] Splits: "+splits+" -> done -> hits: "+hits_a.size()+" -> continue with: "+config.max_hits_per_split_level);
                System.out.println("[INFO] Splits: "+splits+" -> start assembling");
            }

            hits_a = new ArrayList<>( hits_a.subList(0,Math.min(hits_a.size(),config.max_hits_per_split_level)) );

            //List<FinalSimilarityResult> hits_a_expanded = hits_a.parallelStream().map(
            //        hi -> new FinalSimilarityResult(hi.getLeft(), hi.getRight(), hi.getLeft().getAssembled())).collect(Collectors.toList());

            ExecutorService executor_A = Executors.newFixedThreadPool(config.num_threads);
            List<Callable<FinalSimilarityResult>> tasks_Assembly = new ArrayList<Callable<FinalSimilarityResult>>();
            for (final Pair<MoleculeAssembly,Double> hi : hits_a) {
                Callable<FinalSimilarityResult> c = new Callable<FinalSimilarityResult>() {
                    @Override
                    public FinalSimilarityResult call() throws Exception {
                        return new FinalSimilarityResult(hi.getLeft(), hi.getRight(), hi.getLeft().getAssembled());
                    }
                };
                tasks_Assembly.add(c);
            }
            List<FinalSimilarityResult> hits_a_assembled = new ArrayList<>();
            try {
                List<Future<FinalSimilarityResult>> results = executor_A.invokeAll(tasks_Assembly);
                for (Future<FinalSimilarityResult> fr : results) {
                    hits_a_assembled.add( fr.get() );
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }

            progress.accept( (double)(( 2.0*(splits-1) + 1 )  / (2.0*(config.max_splits-1))) );

            if(LOG_LEVEL_SEARCH > 1) {
                System.out.println("[INFO] Splits: "+splits+" -> all assembled");
                System.out.println("[INFO] Splits: "+splits+" -> start descriptor calculation");
            }

            //hits_a_expanded.parallelStream().forEach( pi -> pi.computeFFP() );
            //hits_a_expanded.parallelStream().forEach( pi -> pi.computeSkelSpheres() );

            List<Callable<Integer>> tasks_DescriptorCalc = new ArrayList<Callable<Integer>>();
            for(FinalSimilarityResult fhsi : hits_a_assembled) {
                Callable<Integer> c = new Callable<Integer>() {
                    @Override
                    public Integer call() throws Exception {
                        fhsi.computeFFP();
                        fhsi.computeSkelSpheres();
                        return 1;
                    }
                };
                tasks_DescriptorCalc.add(c);
            }
            try {
                List<Future<Integer>> results = executor_A.invokeAll(tasks_DescriptorCalc);
                for (Future<Integer> fr : results) {
                    fr.get();
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }

            if(LOG_LEVEL_SEARCH > 1) {
                System.out.println("[INFO] Splits: " + splits + " -> all descriptors computed");
            }

            progress.accept( (double)(( 2.0*(splits-1) + 2 )  / (2.0*(config.max_splits-1))) );

            List<Pair<FinalSimilarityResult,Double>> hits_a_expanded_scored = hits_a_assembled.stream().map(
                    hi -> Pair.of(hi, (double)(new DescriptorHandlerSkeletonSpheres()).getSimilarity(query_dh_skelspheres,hi.dh_skelspheres) ) ).collect(Collectors.toList());
            hits_a_expanded_scored.sort( (x,y) -> - Double.compare( x.getRight() , y.getRight()) );
            executor_A.shutdown();
            map_all_hits.put(splits,hits_a_expanded_scored);
        }

        // return list of expanded hits with scores:
        List<Pair<FinalSimilarityResult,Double>> all_hits = new ArrayList<>();
        for(Integer ki : map_all_hits.keySet()) {
            all_hits.addAll(map_all_hits.get(ki));
        }
        all_hits.sort( (x,y) -> -Double.compare(x.getRight(),y.getRight()) );
        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] All done, return "+all_hits.size()+" hits");
        }
        return all_hits;
    }

    /**
     *
     *
     * @param mi
     * @param num_splits
     * @param space3
     * @param cdp
     *
     * @return list containing candidate assemblies
     */
    public static List<Pair<MoleculeAssembly,Double>> findHits(StereoMolecule mi, int num_splits, SynthonSimilaritySpace3 space3, CachedDescriptorProvider cdp, SimilaritySearchConfig3 config ) {

        double initial_tanimoto = config.initial_tanimoto;

        SynthonSimilaritySpace3.TopoConstraints default_constraints = config.topo_constraints;  //new SynthonSimilaritySpace3.TopoConstraints();

        ExecutorService executor_A = Executors.newFixedThreadPool(config.num_threads);

        int nb = mi.getBonds();

        // Protecting access to the concurrent hashmaps:
        // Synchronization via these locks must be used when checking / adding new objects!
        Object LOCK_frag_counts           = new Object();
        Object LOCK_compatible_topo_infos = new Object();
        Object LOCK_compatible_fragtypes  = new Object();

        Map<Pair<SynthonShredder.SplitResult,Integer>,String> frag_idcodes = new ConcurrentHashMap<>();
        Map<String,Integer> frag_counts = new ConcurrentHashMap<>();
        Map<String,Set<SynthonSimilaritySpace3.TopoInfo>> compatible_topo_infos = new ConcurrentHashMap<>();
        Map<String,Set<SynthonSpace.FragType>> compatible_fragtypes             = new ConcurrentHashMap<>();

        // maps frag idcodes to their computed topoinfo.
        Map<String, SynthonSimilaritySpace3.TopoInfo> topo_cache = new HashMap<>();

        // 1. compute all splits
        if(LOG_LEVEL_SEARCH > 1) {System.out.println("[INFO] Start computing splits");}

        List<int[]> splits  = null;
        if(num_splits==0) {
            splits = new ArrayList<>();
            splits.add(new int[0]);
        }
        else {
            List<int[]> combi_list = CombinationGenerator.getAllOutOf(nb, num_splits);
            // !! returns null if b > a..
            if(combi_list==null) { return new ArrayList<>(); }
            splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());
        }

        List<SynthonShredder.SplitResult> labeled_splits = Collections.synchronizedList(new ArrayList<>());

        List<Callable<Integer>> tasks_trySplit = new ArrayList<>();
        for(int zi=0;zi<splits.size();zi++) {
            final int[] split_i = splits.get(zi);
            Callable<Integer> trySplit = new Callable<Integer>(){
                @Override
                public Integer call() throws Exception {
                    SynthonShredder.SplitResult sri = SynthonShredder.trySplit(mi, split_i, 3);
                    if (sri == null) {
                        return 0;
                    }
                    List<SynthonShredder.SplitResult> splits_v = SynthonShredder.getAllSplitVariations(sri, num_splits);
                    for (SynthonShredder.SplitResult svi : splits_v) {
                        labeled_splits.add(svi);
                    }
                    return 1;
                }
            };
            tasks_trySplit.add(trySplit);
        }

        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_trySplit);
            for (Future<Integer> fr : results) {
                //hits_a_assembled.add( fr.get() );
                fr.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] Computing splits -> done, Labeled splits: "+labeled_splits.size());
            System.out.println("[INFO] Start computing topo data");
        }

        // 2. compute all topo data:
        List<Callable<Integer>> tasks_computeTopo = new ArrayList<>();
        for(int zi=0; zi<labeled_splits.size(); zi++)
        {
            final SynthonShredder.SplitResult svi = labeled_splits.get(zi);
            Callable<Integer> task_computeTopo = new Callable<Integer>(){
                @Override
                public Integer call() throws Exception {
                    for(int fzi = 0; fzi<svi.fragments.length;fzi++) {
                        StereoMolecule ffi = svi.fragments[fzi];
                        String idc_i = ffi.getIDCode();
                        frag_idcodes.put(Pair.of(svi,fzi),idc_i);

                        synchronized(LOCK_frag_counts) {
                            if (!frag_counts.containsKey(idc_i)) {
                                frag_counts.put(idc_i, 1);
                            } else {
                                frag_counts.put(idc_i, frag_counts.get(idc_i) + 1);
                            }
                        }
                        synchronized(LOCK_compatible_topo_infos) {
                            if (!compatible_topo_infos.containsKey(idc_i)) {
                                SynthonSimilaritySpace3.TopoInfo tinfi = SynthonSimilaritySpace3.computeTopoInfoForMolecule(ffi);
                                topo_cache.put(idc_i, tinfi);
                                compatible_topo_infos.put(idc_i, new HashSet<>(default_constraints.createCompatibleTopoConstraints(tinfi)));
                            }
                        }

                        synchronized(LOCK_compatible_fragtypes) {
                            if (!compatible_fragtypes.containsKey(idc_i)) {
                                compatible_fragtypes.put(idc_i, new HashSet<>());
                                for (SynthonSimilaritySpace3.TopoInfo tinf : compatible_topo_infos.get(idc_i)) {
                                    if (space3.mFragmentToposInv.containsKey(tinf)) {
                                        compatible_fragtypes.get(idc_i).addAll(space3.mFragmentToposInv.get(tinf));
                                    }
                                }
                            }
                        }
                        //if(!compatible_fragtypes.containsKey(idc_i)) {compatible_fragtypes.put(idc_i,new HashSet<>());}
                        //for(  )
                    }
                    return 1;
                }
            };
            tasks_computeTopo.add(task_computeTopo);
        }

        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_computeTopo);
            for (Future<Integer> fr : results) {
                //hits_a_assembled.add( fr.get() );
                fr.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] Computing topo data -> done");
            System.out.println("[INFO] Start finding valid split->synthon maps");
        }


        // 3. loop over all splits / reactions / split -> fragtype mappings and determine all mappings that agree with topo

        // this will contain all valid splits with their mappings. Mapping int value indicates index of splitresult.fragment array.
        List<Pair<SynthonShredder.SplitResult,Map<Integer,SynthonSpace.FragType>>> all_maps = Collections.synchronizedList(new ArrayList<>());

        List<Callable<Integer>> tasks_findValidMaps = new ArrayList<>();

        for(int zi=0;zi<labeled_splits.size();zi++) {
            final SynthonShredder.SplitResult sri = labeled_splits.get(zi);
            Callable<Integer> task_findValidMap = new Callable<Integer>() {
                @Override
                public Integer call() throws Exception {
                    List<String> idcs = Arrays.stream(sri.fragments).map(fi -> fi.getIDCode()).collect(Collectors.toList());

                    // 1. check if any reactand has no matching frag types
                    List<Set<SynthonSpace.FragType>> compatible_frag_types = idcs.stream().map(idi -> compatible_fragtypes.get(idi)).collect(Collectors.toList());
                    if (compatible_frag_types.stream().anyMatch(cfi -> (cfi == null) ? true : cfi.isEmpty())) {
                        // then we are done :)
                        return 0;
                    }

                    // 2. check for reaction ids that appear in any:
                    List<String> rxns_in_any = compatible_frag_types.stream().flatMap(fi -> fi.stream().map(xi -> xi.rxn_id)).distinct().collect(Collectors.toList());

                    if (LOG_LEVEL_A >= 2) {
                        System.out.println("sr: " + sri.toString());
                        System.out.println("rxns: " + rxns_in_any);
                    }

                    // 3. loop over these reactions:
                    for (int ri = 0; ri < rxns_in_any.size(); ri++) {

                        String rxi = rxns_in_any.get(ri);

                        if (space3.fragment_map.get(rxi).keySet().size() != sri.fragments.length) {
                            // then we are done :) we dont consider partial matches here
                            continue;
                        }

                        if (LOG_LEVEL_A >= 2) {
                            System.out.println("rxn " + rxi + " -> 1");
                        }

                        // collect all rxn frag types
                        List<List<SynthonSpace.FragType>> compatible_frag_types_rxn =
                                compatible_frag_types.stream().map(ftis -> ftis.stream().filter(fft -> fft.rxn_id.equals(rxi)).collect(Collectors.toList())).collect(Collectors.toList());


                        // check if rxn appears in all compatible frag types:
                        if (compatible_frag_types_rxn.stream().anyMatch(fxi -> fxi.isEmpty())) {
                            // then we are done :)
                            continue;
                        }

                        if (LOG_LEVEL_A >= 2) {
                            System.out.println("rxn " + rxi + " -> 2");
                        }

                        // check if all different frag types of the rxn appear:
                        List<SynthonSpace.FragType> available_frag_types = compatible_frag_types_rxn.stream().flatMap(ffi -> ffi.stream()).distinct().collect(Collectors.toList());
                        if (available_frag_types.size() != space3.fragment_map.get(rxi).keySet().size()) {
                            // then we are done :)
                            continue;
                        }

                        if (LOG_LEVEL_A >= 2) {
                            System.out.println("rxn " + rxi + " -> 3");
                        }

                        // now we have to consider all mappings from fragments to fts
                        //if(sri.fragments[0].getIDCode().equals("fi{pR@HNrRu^`NdfYu][V]mhrRdyjZ@JZjijB@@")) {
                        //    System.out.println("mkay..");
                        //}

                        List<Map<Integer, SynthonSpace.FragType>> all_maps_i = computeAllPossibleSynthonMaps(available_frag_types, compatible_frag_types_rxn);

                        // and now we have to actually start searching for fragments that look similar ;)
                        for (Map<Integer, SynthonSpace.FragType> mapi : all_maps_i) {
                            all_maps.add(Pair.of(sri, mapi));
                        }

                        //return 1;
                    }
                    return 1;
                }
            };
            tasks_findValidMaps.add(task_findValidMap);
        }

        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_findValidMaps);
            for (Future<Integer> fr : results) {
                //hits_a_assembled.add( fr.get() );
                fr.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] valid maps -> done, totla srs/maps: "+all_maps.size());
            System.out.println("[INFO] Start computing initial hits");
        }

        //System.out.println("Total srs/maps: "+all_maps.size());



        // Now we have all maps and we search for the most similar initial hits
        //List<Pair<Pair<SynthonShredder.SplitResult,Map<Integer, SynthonSpace.FragType>>,Pair<SynthonSpace.FragId,Double>>> initial_hits = new ArrayList<>();
        List<SimilarityHit3> initial_hits = Collections.synchronizedList( new ArrayList<>() );

        List<Callable<Integer>> tasks_computeInitialHits = new ArrayList<>();

        for(Pair<SynthonShredder.SplitResult,Map<Integer, SynthonSpace.FragType>> mapi : all_maps) {
            final SynthonShredder.SplitResult sri = mapi.getLeft();
            final Map<Integer, SynthonSpace.FragType> mapped_fts = mapi.getRight();

            final SynthonSpace.FragType initial_ft = mapped_fts.get(0); // initial hit fragment (largest)

            final Map<BitSet,List<SynthonSpace.FragId>> all_frags_i = space3.fragment_map_2.get(initial_ft.rxn_id).get(initial_ft.frag);
            final Set<SynthonSimilaritySpace3.TopoInfo> topos = compatible_topo_infos.get( frag_idcodes.get( Pair.of( sri , 0) ) );

//            if(initial_ft.rxn_id.equals("m_11bbc")) {
//                System.out.println("mkay..");
//                System.out.println(HyperspaceUtils.idcodeToSmiles(frag_idcodes.get( Pair.of( sri , 0) )));
//                System.out.println(HyperspaceUtils.idcodeToSmiles(frag_idcodes.get( Pair.of( sri , 1) )));
//            }

            Callable<Integer> task_findInitialHits = new Callable<Integer>(){
                @Override
                public Integer call() throws Exception {
                    for( SynthonSimilaritySpace3.TopoInfo tinfi : topos ) {
                        int cnt_initial_hits_tinfi = 0;
                        FastSimilaritySearcher fssi = space3.mFastSimilaritySearchers.get(initial_ft).get(tinfi);
                        if(fssi != null) {
                            BitSet bsi = cdp.getFP_cached(sri.fragments[0]);
                            // NOTE: here we assume that the initial fragment is ALWAYS bigger than the "dont test similarity" threshold..
                            int cardinality_i = bsi.cardinality();
                            int hd_tani = (int) Math.floor( cardinality_i * ( (1.0/initial_tanimoto) - 1.0 ) );
                            List<LSHProvider.Hit2> hits = fssi.mANN.exactFindNearestNeighbors2(bsi,hd_tani,1000);
                            for(LSHProvider.Hit2 hi : hits) {
                                double tani_i = LSHProvider.tanimoto_similarity(bsi,hi.bs);
                                if(tani_i>=initial_tanimoto) {
                                    for (SynthonSpace.FragId fidi : all_frags_i.get(hi.bs)) {
                                        // add initial hit
                                        //initial_hits.add(Pair.of(mapi, Pair.of(fidi, tani_i)));
                                        SynthonSimilaritySpace3.TopoInfo split_topo = topo_cache.get(frag_idcodes.get( Pair.of( sri , 0)));
                                        SimilarFragmentHit sfh = new SimilarFragmentHit(sri,0,split_topo,fidi,tinfi,tani_i);
                                        SimilarityHit3 sim_hit = new SimilarityHit3(sri,mapped_fts,sfh);
                                        initial_hits.add(sim_hit);
                                        cnt_initial_hits_tinfi++;
                                    }
                                }
                            }
                        }
                        //System.out.println("initial hits: "+((cnt_initial_hits_tinfi==0)?"NONE":""+cnt_initial_hits_tinfi));
                    }
                    return 1;
                }
            };

            tasks_computeInitialHits.add(task_findInitialHits);
        }



        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_computeInitialHits);
            for (Future<Integer> fr : results) {
                //hits_a_assembled.add( fr.get() );
                fr.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] Finding initial hits done, hits: "+initial_hits.size());
            System.out.println("[INFO] Start expanding initial hits");
        }



        // now expand initial hits
        int expansions_per_frag_2s = 20;
        int expansions_per_frag_3s = 6;

        initial_hits.sort( (x,y) -> - Double.compare( x.initial_hit.tanimoto , y.initial_hit.tanimoto ) );


        List<Callable<Integer>> tasks_expandHits = new ArrayList<>();

        for(int zi=0;zi<initial_hits.size();zi++) {
            //Pair<Pair<SynthonShredder.SplitResult,Map<Integer, SynthonSpace.FragType>>,Pair<SynthonSpace.FragId,Double>> hit_i = initial_hits.get(zi);
            SimilarityHit3 hit_i = initial_hits.get(zi);
            SynthonShredder.SplitResult sri = hit_i.split;
            int num_expansions_per_frag = expansions_per_frag_2s;
            if(sri.fragments.length>=3) { num_expansions_per_frag = expansions_per_frag_3s;}
            int final_num_expansions_per_frag = num_expansions_per_frag;

            Callable<Integer> task_expandHit = new Callable<Integer>(){
                @Override
                public Integer call() throws Exception {

                    // we have to expand all split result frags except frag 0:
                    //Map<Integer,List<Pair<SynthonSpace.FragId,Double>>> expansions = new ArrayList<>();
                    for(int frag = 1 ; frag < sri.fragments.length ; frag++) {
                        SynthonSpace.FragType fti = hit_i.map.get(frag);
                        BitSet bsi = cdp.getFP_cached(sri.fragments[frag]);

                        SynthonSimilaritySpace3.TopoInfo split_topo = topo_cache.get(frag_idcodes.get( Pair.of( sri , frag)));

                        Map<BitSet,List<SynthonSpace.FragId>> all_frags_i = space3.fragment_map_2.get(fti.rxn_id).get(fti.frag);

                        //List<LSHProvider.Hit2> hits_frag = new ArrayList<>();
                        // find all hit2 objects:
                        List<SimilarFragmentHit> frag_all_hits = new ArrayList<>(); // this should rather be a priority queue..

                        Set<SynthonSimilaritySpace3.TopoInfo> topos = compatible_topo_infos.get( frag_idcodes.get( Pair.of( sri , frag) ) );
                        for( SynthonSimilaritySpace3.TopoInfo tinfi : topos ) {
                            FastSimilaritySearcher fssi = space3.mFastSimilaritySearchers.get(fti).get(tinfi);
                            if(fssi!=null) {
                                List<LSHProvider.Hit2> hits_frag_i = fssi.mANN.exactFindKNearestNeighbors2_Tanimoto(bsi,final_num_expansions_per_frag);
                                for(int zh = 0; zh<hits_frag_i.size() ;zh++) {
                                    LSHProvider.Hit2 hi = hits_frag_i.get(zh);
                                    for (SynthonSpace.FragId fidi : all_frags_i.get(hi.bs)) {
                                        SimilarFragmentHit sfh = new SimilarFragmentHit(sri, frag, split_topo, fidi, tinfi, hits_frag_i.get(zh).tanimoto_sim);
                                        frag_all_hits.add(sfh);
                                    }
                                }
                            }
                        }

                        frag_all_hits.sort( (x,y) -> - Double.compare(x.tanimoto,y.tanimoto) );
                        hit_i.addExtensionHits(frag,frag_all_hits);

                        // clean up hits and construct
                        //List<LSHProvider.Hit2> hits_to_use_sorted = all_hits.stream().distinct().sorted( (x, y) -> -Double.compare(x.tanimoto_sim,y.tanimoto_sim) ).collect(Collectors.toList());
                        //for(int zh = 0; zh<num_expansions_per_frag ;zh++) {
                        //LSHProvider.Hit2 hi2 = hits_to_use_sorted.get(zh);
                        //    SynthonSimilaritySpace3.TopoInfo split_topo = topo_cache.get(frag_idcodes.get( Pair.of( sri , frag)));
                        //    SimilarFragmentHit sfh = new SimilarFragmentHit(sri,frag,split_topo,fidi,tinfi,tani_i);
                        //}

                    }

                    return 1;
                }
            };
            tasks_expandHits.add(task_expandHit);
        }


        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_expandHits);
            for (Future<Integer> fr : results) {
                //hits_a_assembled.add( fr.get() );
                fr.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if(LOG_LEVEL_SEARCH > 1) {
            System.out.println("[INFO] Done expanding hits!");
        }

        // now we have nicely expanded hits :)
        //System.out.println("mkay..");


        // now decide what to assemble..
        double score_exponent = 2.0;
        Map<String,Pair<MoleculeAssembly,Double>> assemblies = new HashMap<>();


        for(int zi=0;zi<initial_hits.size();zi++) {
            SimilarityHit3 hi = initial_hits.get(zi);
            List<MoleculeAssembly> assemblies_i  = hi.getAssemblies();
            for(MoleculeAssembly masi : assemblies_i) {
                String hash_masi  = masi.getFragHash();
                double score_masi = masi.getScore_A(score_exponent);
                if(assemblies.containsKey(hash_masi)) {
                    if(assemblies.get(hash_masi).getRight() < score_masi) {
                        assemblies.put(hash_masi,Pair.of(masi,score_masi));
                    }
                }
                else {
                    assemblies.put(hash_masi,Pair.of(masi,score_masi));
                }
            }
        }

        List<Pair<MoleculeAssembly,Double>> unique_assemblies = new ArrayList<>( assemblies.values() );
        unique_assemblies.sort( (x,y) -> -Double.compare( x.getRight() , y.getRight() ) );

        executor_A.shutdown();

        System.out.println("mkay.. :) assembly candidates: "+assemblies.size());
        return unique_assemblies;
    }

    /**
     *
     *
     * @param all_frag_types
     * @param compatible_frag_types
     * @return
     */
    public static List<Map<Integer, SynthonSpace.FragType>> computeAllPossibleSynthonMaps(
            List<SynthonSpace.FragType> all_frag_types,
            List<List<SynthonSpace.FragType>> compatible_frag_types) {

        List<Map<Integer,SynthonSpace.FragType>> all_possible_maps = new ArrayList<>();

        List<Integer> indeces = new ArrayList<>(); for(int zi=0;zi<all_frag_types.size();zi++){indeces.add(zi);}
        List<List<Integer>> all_maps = SynthonSpace.all_permutations(indeces);

        // try all of these maps:
        for(List<Integer> map_i : all_maps) {
            // try to map the fragments to the frag types:
            boolean map_ok = true;
            for( int zi=0;zi<map_i.size();zi++ ) {
                SynthonSpace.FragType ft_i = all_frag_types.get( map_i.get( zi ) );
                if(compatible_frag_types.get(zi).contains(ft_i)) {
                    // ok :)
                }
                else {
                    map_ok = false;
                    break;
                }
            }
            if(map_ok) {
                // add map
                Map<Integer,SynthonSpace.FragType> map_new = new HashMap<>();
                for( int zi=0;zi<map_i.size();zi++ ) {
                    SynthonSpace.FragType ft_i = all_frag_types.get( map_i.get( zi ) );
                    //map_new.put(ft_i,zi);
                    map_new.put(zi,ft_i);
                }
                all_possible_maps.add(map_new);
            }
        }

        return all_possible_maps;
    }

    public static class SimilarityHit3 implements Serializable {
        SynthonShredder.SplitResult               split;
        Map<Integer, SynthonSpace.FragType>       map;
        SimilarFragmentHit                        initial_hit;
        Map<Integer,List<SimilarFragmentHit>>     extension_hits;
        public SimilarityHit3(SynthonShredder.SplitResult sri,
                              Map<Integer, SynthonSpace.FragType> map,
                              SimilarFragmentHit initial) {
            this.split = sri;
            this.map   = map;
            this.initial_hit = initial;
        }

        public SimilarFragmentHit getInitialHit() {
            return this.initial_hit;
        }

        public void addExtensionHits(int frag, List<SimilarFragmentHit> hits) {
            if(this.extension_hits==null) {this.extension_hits = new HashMap<>();}
            if(!this.extension_hits.containsKey(frag)) { this.extension_hits.put(frag,hits); }
            else{ this.extension_hits.get(frag).addAll( hits );}
            // sort now? maybe not
        }

        public List<MoleculeAssembly> getAssemblies() {
            List<MoleculeAssembly> assemblies = new ArrayList<>();
            List<Integer> extension_keys = new ArrayList<>(extension_hits.keySet());
            if(extension_keys.size()==1) {
                for (SimilarFragmentHit hi : this.extension_hits.get(extension_keys.get(0)) ) {
                    List<SimilarFragmentHit> parts = new ArrayList<>();
                    parts.add(initial_hit);
                    parts.add(hi);
                    MoleculeAssembly asi = new MoleculeAssembly(this,parts);
                    assemblies.add(asi);
                }
            }
            else if(extension_keys.size()==2) {
                for (SimilarFragmentHit hia : this.extension_hits.get(extension_keys.get(0)) ) {
                    for (SimilarFragmentHit hib : this.extension_hits.get(extension_keys.get(1)) ) {
                        List<SimilarFragmentHit> parts = new ArrayList<>();
                        parts.add(initial_hit);
                        parts.add(hia);
                        parts.add(hib);
                        MoleculeAssembly asi = new MoleculeAssembly(this,parts);
                        assemblies.add(asi);
                    }
                }
            }
            return assemblies;
        }
    }
    public static class SimilarFragmentHit implements Serializable {
        public final SynthonShredder.SplitResult       split;
        public final int                               split_idx;
        public final SynthonSimilaritySpace3.TopoInfo  split_topo;
        public final SynthonSpace.FragId               frag;
        public final SynthonSimilaritySpace3.TopoInfo  frag_topo;
        public final double tanimoto;
        public SimilarFragmentHit(SynthonShredder.SplitResult split, int split_idx,
                                  SynthonSimilaritySpace3.TopoInfo split_topo,
                                  SynthonSpace.FragId frag, SynthonSimilaritySpace3.TopoInfo frag_topo,
                                  double tanimoto) {
            this.split       = split;
            this.split_idx   = split_idx;
            this.split_topo  = split_topo;
            this.frag        = frag;
            this.frag_topo   = frag_topo;
            this.tanimoto    = tanimoto;
        }
    }

    /**
     *
     */
    public static class MoleculeAssembly implements Serializable {
        SimilarityHit3 hit;
        /**
         * fragments, in canonical order, based on fragment idcodes.
         */
        List<SimilarFragmentHit> parts;

        public MoleculeAssembly(SimilarityHit3 hit, List<SimilarFragmentHit> parts) {
            this.hit = hit;
            this.parts = new ArrayList<>(parts);
            this.parts.sort( (x,y) -> (x.frag.idcode.compareTo( y.frag.idcode )) );
        }

        public List<SimilarFragmentHit> getFragmentHits() {
            return this.parts;
        }

        public SimilarityHit3 getSimilarityHit3() {
            return this.hit;
        }

        /**
         * Returns hash for 1. rxn, and 2. the set of fragment idcodes (independent
         * of mapping), i.e different frag hash guarantees that at least one fragment
         * is different.
         * However, it does obviously not guarantee that the final molecule is
         * different, becuse there might be different sets of fragments to
         * assemble the same molecule.
         *
         * @return
         */
        public String getFragHash() {
            List<String> idcds = parts.stream().map( pi -> pi.frag.idcode ).collect(Collectors.toList());
            idcds.sort((x,y) -> x.compareTo(y));
            List<String> hash_parts = new ArrayList<>();
            hash_parts.add(hit.map.get(0).rxn_id);
            hash_parts.addAll(idcds);
            return String.join( "##" , hash_parts );
        }

        /**
         * fragment similarities are scored based on exp of size
         *
         * @param exp
         * @return
         */
        public double getScore_A(double exp) {
            int size_tot = -1;
            try {
                size_tot = this.parts.stream().mapToInt(pi -> pi.split_topo.size).sum();
            }
            catch(Exception ex) {
                System.out.println("mkay..");
                return 0.0;
            }
            double[] mix = new double[parts.size()];
            double   mix_sum = 0;
            for(int zi=0;zi<mix.length;zi++) {
                double mi = Math.pow( (1.0*parts.get(zi).frag_topo.size) , exp) / size_tot;
                mix[zi] =  mi;
                mix_sum += mi;
            }

            double score = 0;
            for(int zi=0;zi<mix.length;zi++) {
                score += parts.get(zi).tanimoto * ( mix[zi] / mix_sum );
            }
            return score;
        }

        public StereoMolecule getAssembled() {
            List<StereoMolecule> mols = this.parts.stream().map( pi -> HyperspaceUtils.parseIDCode(pi.frag.idcode)).collect(Collectors.toList());
            StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(mols);
            assembled.ensureHelperArrays(Molecule.cHelperCIP);
            return assembled;
        }
    }

    public static class FinalSimilarityResult implements Serializable {
        public final MoleculeAssembly assembly;
        public final double           assembly_score;
        public final StereoMolecule   mol;

        private BitSet   dh_ffp;
        private byte[]   dh_skelspheres;

        public FinalSimilarityResult(MoleculeAssembly asm, double asm_score, StereoMolecule mi) {
            this.assembly        = asm;
            this.assembly_score  = asm_score;
            this.mol             = mi;
        }

        public BitSet getFFP() {return this.dh_ffp;}
        public byte[] getSkelSpheres() {return this.dh_skelspheres;}

        public void computeFFP() { this.dh_ffp =  BitSet.valueOf( (new DescriptorHandlerLongFFP512()).createDescriptor(this.mol) ); }
        public void computeSkelSpheres() { this.dh_skelspheres = (new DescriptorHandlerSkeletonSpheres()).createDescriptor(this.mol); }
    }

    public static void writeResultsToCSV(String filename, BitSet mi_descriptor_ffp, List<Pair<FinalSimilarityResult,Double>> data) {
        BufferedWriter data_out = null;
        try {
            data_out = new BufferedWriter(new FileWriter(filename));
            //System.out.println("Hit[idcode],Sim,Rxn,SimA2");
            data_out.write("Hit[idcode],SimSkelSpheres,SimFragFp,Rxn,SimA2"+"\n");
            for(Pair<FinalSimilarityResult,Double> asi : data) {
                List<String> result_parts = new ArrayList<>();
                result_parts.add(asi.getLeft().mol.getIDCode());
                result_parts.add(""+asi.getRight());
                result_parts.add(""+LSHProvider.tanimoto_similarity(mi_descriptor_ffp,asi.getLeft().dh_ffp));
                result_parts.add(asi.getLeft().assembly.hit.initial_hit.frag.rxn_id);
                result_parts.add(""+asi.getLeft().assembly.getScore_A(2));
                //System.out.println(String.join(",",result_parts));
                data_out.write(String.join(",",result_parts)+"\n");
            }
            data_out.flush();
            data_out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
