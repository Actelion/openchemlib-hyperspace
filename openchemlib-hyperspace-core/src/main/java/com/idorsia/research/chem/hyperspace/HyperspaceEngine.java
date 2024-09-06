package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Function;

public class HyperspaceEngine {

    public static final class SubstructureSearchSpace {
        public final String name;
        public final String version;
        public final String provider;
        public final String file;
        public final SynthonSpace space;
        public  SubstructureSearchSpace(String name, String version, String provider, String file, SynthonSpace space) {
            this.name = name; this.version=version; this.provider=provider; this.file = file; this.space = space;
        }
    }

    public static final class SimilaritySearchSpace {
        public final String name;
        public final String version;
        public final String provider;
        public final String file;
        public final SynthonSimilaritySpace space;
        public  SimilaritySearchSpace(String name, String version, String provider, String file, SynthonSimilaritySpace space) {
            this.name = name; this.version=version; this.provider=provider; this.file = file; this.space = space;
        }
    }

    Map<String,SubstructureSearchSpace> mSpaces_SSS        = new HashMap<>();
    Map<String,SimilaritySearchSpace>   mSpaces_Similarity = new HashMap<>();

    public void loadSpace_SS(String space_description_file) throws Exception {
        File fi_prop = new File(space_description_file);
        loadSpace_SS(new FileInputStream(fi_prop));
    }

    public void loadSpace_SS(InputStream in_prop) throws Exception {

        Properties pi = new Properties();
        pi.load( in_prop );

        if(!pi.getProperty("type").equals("sss")) {
            throw new Exception("Space description file must be for type \"sss\"");
        }

        String name     = pi.getProperty("name");
        String version  = pi.getProperty("version");
        String provider = pi.getProperty("provider");
        String file     = pi.getProperty("file");

        System.out.println("\nInitialize SSS Space: " + name + " from file "+file);
        System.out.println("Name          = " + name );
        System.out.println("Version       = " + version);
        System.out.println("Provider      = " + provider);

        System.out.print("\n\nLoad data..");
        ObjectInputStream in = new ObjectInputStream( new FileInputStream( file ) );
        SynthonSpace space = (SynthonSpace) in.readObject();
        space.initAfterJavaDeserialization();
        System.out.print(" done!\n\n");

        if(mSpaces_SSS.containsKey(name)) {
            String name_old = name;
            name = name + "["+(new Random()).nextInt()+"]";
            System.out.println("[WARN] sss space name "+name_old+" appears multiple times! -> rename to ");
        }
        mSpaces_SSS.put(name, new SubstructureSearchSpace(name,version,provider,file,space) );
        System.out.println("Initialize SSS Space: " + name + " from file "+file + "  --> DONE!");

        System.out.println("Statistics for Synthon Space:");
        printSpaceStats(space);
        System.out.println();

        System.out.println();
    }

    public void loadSpace_Similarity(String space_description_file) throws Exception {
        File fi_prop = new File(space_description_file);
        loadSpace_Similarity(new FileInputStream(fi_prop));
    }

    public void loadSpace_Similarity(InputStream space_in) throws Exception {

        Properties pi = new Properties();
        pi.load( space_in );

        if(!pi.getProperty("type").equals("similarity")) {
            throw new Exception("Space description file must be for type \"similarity\"");
        }

        String name     = pi.getProperty("name");
        String version  = pi.getProperty("version");
        String provider = pi.getProperty("provider");
        String file     = pi.getProperty("file");

        int lsh_functions = Integer.parseInt(pi.getProperty("lsh_num_functions"));
        int lsh_bits      = Integer.parseInt(pi.getProperty("lsh_projection_bits"));


        System.out.println("Initialize Similarity Space: " + name + " from file "+file);
        System.out.println("Name          = " + name );
        System.out.println("Version       = " + version);
        System.out.println("Provider      = " + provider);
        System.out.println();
        System.out.println("LSH_Functions = " + lsh_functions);
        System.out.println("LSH_Bits      = " + lsh_bits);
        System.out.println();

        System.out.print("\nLoad data..");
        ObjectInputStream in = new ObjectInputStream( new FileInputStream( file ) );
        SynthonSpace ss_space            = (SynthonSpace) in.readObject();
        ss_space.initAfterJavaDeserialization();
        System.out.print(" done!\n\n");

        System.out.println("Init approximate nearest neighbor datastructures");

        SynthonSimilaritySpace space   = new SynthonSimilaritySpace(ss_space);
        System.out.println("Init ANNs");
        space.initApproximateNN(lsh_functions,lsh_bits,8);
        System.out.println("Init ANNs -> done!");

        if(mSpaces_Similarity.containsKey(name)) {
            String name_old = name;
            name = name + "["+(new Random()).nextInt()+"]";
            System.out.println("[WARN] similarity space name "+name_old+" appears multiple times! -> rename to ");
        }
        mSpaces_Similarity.put(name, new SimilaritySearchSpace(name,version,provider,file,space) );
        System.out.println("Initialize Similarity Space: " + name + " from file "+file + "  --> DONE!\n");

        System.out.println("Statistics for Synthon Space:");
        printSpaceStats(space);
        System.out.println();

    }

    public static void setLogLevelToQuiet() {
        SynthonSpace.setLogLevel(0,0);
        SynthonSimilaritySpace.setLogLevel_B(0,0);
        SynthonSpaceExplorer.setLogLevel(0,0,0);
        SynthonSpaceSimilarityExplorer.setLogLevel(0,0);
    }

    /**
     * This should be threadsafe.
     *
     */
    public List<SynthonSpaceExplorer.FinalizedHit> findSubstructures(String space_name , StereoMolecule mi ,
                                                                     int max_num_full_extensions, int threads) {
        SubstructureSearchSpace ss_space = mSpaces_SSS.get(space_name);
        SynthonSpace space = ss_space.space;

        SynthonSpaceExplorer explorer = new SynthonSpaceExplorer(space,mi);
        explorer.screen_default(3, 3, max_num_full_extensions, threads);
        return explorer.getAllResults();
    }

    /**
     * This should be threadsafe.
     *
     * @param space_name
     * @param mi
     * @param min_tanimoto_sim_initial
     * @param primary_hits_screened
     * @param primary_hits_considered
     * @param num_hit_expansions
     * @param allowed_rxns pass null if you don't want to exclude any rxns.
     * @param threads
     * @return
     */
    public List<SynthonSpaceSimilarityExplorer.FinalizedSimilarityHit> findSimilarStructures(String space_name , StereoMolecule mi,
                                                                                             double min_tanimoto_sim_initial,
                                                                                             int primary_hits_screened,
                                                                                             int primary_hits_considered,
                                                                                             int num_hit_expansions,
                                                                                             List<String> allowed_rxns,
                                                                                             int threads) {
        SimilaritySearchSpace sim_space = mSpaces_Similarity.get(space_name);
        SynthonSimilaritySpace space = sim_space.space;

        //SynthonSpaceSimilarityExplorer explorer = new SynthonSpaceSimilarityExplorer(space,mi);

        // translate for two expansions to "per expansion" value:
        int hit_expansion_values[] = new int[]{num_hit_expansions, (int) Math.ceil( Math.sqrt( num_hit_expansions ))};

        SynthonSpaceSimilarityExplorer explorer = new SynthonSpaceSimilarityExplorer(space,mi,
                128,primary_hits_screened,
                primary_hits_considered,hit_expansion_values,
                allowed_rxns);
        //explorer.screen_default(3, 3, );

        Function<Integer,Integer> f_initial_dist = (Integer bits) -> Integer.valueOf( (int) Math.ceil(  ((double)bits.intValue()) * (1.0-min_tanimoto_sim_initial) ) + 2 );

        Random r = new Random();
        List<SynthonSpaceSimilarityExplorer.SimilarityHit> hits = new ArrayList<>();

//        explorer.sampleSimilar_new_default(r,mi,3,3,f_initial_dist,1024,
//                num_hit_expansions,8000,hits,threads);
        explorer.sampleSimilar_new_default(r,mi,3,3,f_initial_dist,
                hits,threads);


        // finalize.. (aka assemble)
        CachedDescriptorProvider cdh = new CachedDescriptorProvider(space.getDescriptorHandler().getInfo().shortName);

        CachedStereoMoleculeProvider csmp = new CachedStereoMoleculeProvider();
        List<SynthonSpaceSimilarityExplorer.FinalizedSimilarityHit> finalized_hits = Collections.synchronizedList(new ArrayList<>());
        Set<Set<SynthonSpace.FragId>> already_assembled_sets = new HashSet<>();

        ExecutorService main_pool = Executors.newFixedThreadPool(threads);
        List<Future> tasks_assemble = new ArrayList<>();
        for(SynthonSpaceSimilarityExplorer.SimilarityHit hi : hits) {
            Runnable ri = new Runnable(){
                @Override
                public void run() {
                    finalized_hits.addAll( hi.assembleAllStructures(cdh,csmp,already_assembled_sets,mi) );
                }
            };
            tasks_assemble.add(main_pool.submit(ri));
        }
        // wait for tasks to finish:
        // wait for pool:
        int cnt_fin = 0;

        for(Future ft : tasks_assemble) {
            try {
                ft.get();
                cnt_fin++;
                if(false) {
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

        return finalized_hits;
    }

    public List<SubstructureSearchSpace> getAvailableSubstructureSearchSpaceObjects() {
        return new ArrayList<>(this.mSpaces_SSS.values());
    }

    public Map<String,SubstructureSearchSpace> getAvailableSubstructureSearchSpaces() {
        return this.mSpaces_SSS;
    }

    public List<SimilaritySearchSpace> getAvailableSimilaritySearchSpaceObjects() {
        return new ArrayList<>(this.mSpaces_Similarity.values());
    }

    public Map<String,SimilaritySearchSpace> getAvailableSimilaritySearchSpaces() {
        return this.mSpaces_Similarity;
    }


    public static void printSpaceStats(SynthonSpace space) {
        long num_feasible_products    = space.fragment_map.values().stream().mapToLong( vi -> vi.values().stream().mapToLong( vj -> vj.size() ).reduce( (a,b) -> a*b ).getAsLong()).sum();
        long num_synthons_non_unique  = space.fragment_map.values().stream().mapToLong( vi -> vi.values().stream().mapToLong( vj -> vj.size() ).reduce( (a,b) -> a+b ).getAsLong()).sum();

        System.out.println("Number of Synthon Rxns:              " + space.fragment_map.keySet().size() );
        System.out.println("Number of Synthons (non-unique):     " + num_synthons_non_unique);
        System.out.println("Number of Feasible Structures:       " + num_feasible_products );
    }

    public static void main(String args[]) {

        HyperspaceEngine engine = new HyperspaceEngine();

        try {
            engine.loadSpace_SS("data/enamine_real_space_sss_a.space");
            engine.loadSpace_Similarity("data/enamine_real_space_similarity_a.space");

            IDCodeParser icp     = new IDCodeParser();
            StereoMolecule mi    = new StereoMolecule();
            StereoMolecule mi_2  = new StereoMolecule();
            icp.parse(mi,"fakq`BHJBAGhRclbfbRabTTjTRWAVBwkP@USMP@D@A@");
            icp.parse(mi_2,"ehzzNHDBKOCiYJd\\RQ_PI@z\\^rJJIKPqQQPqIQZJIIJIZyQI[I`SrJZzEfuvBBFZjjfjBjAkJB@hh@@@");

            System.out.println("Substructure Search:");
            List<SynthonSpaceExplorer.FinalizedHit> hits = engine.findSubstructures( engine.getAvailableSubstructureSearchSpaceObjects().get(0).name , mi , 64 , 4 );
            System.out.println("Hits: "+hits.size());

            System.out.println("Similarity Search:");
            List<SynthonSpaceSimilarityExplorer.FinalizedSimilarityHit> hits_sim = engine.findSimilarStructures( engine.getAvailableSimilaritySearchSpaceObjects().get(0).name , mi_2 , 0.7 , 200000, 5000 ,64 , null, 4);
            System.out.println("Hits: "+hits_sim.size());

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
