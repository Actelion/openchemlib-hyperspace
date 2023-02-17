package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import javafx.concurrent.Task;

import java.io.IOException;
import java.io.Serializable;
import java.sql.SQLOutput;
import java.util.*;
import java.util.concurrent.*;
import java.util.function.Consumer;
import java.util.stream.Collectors;


/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace
 */


public class SynthonSimilaritySpace3 extends SynthonSpace {




    public static class TopoConstraints {
        /**
         * From this number of atoms on we check for similarity
         * (atoms include connector atoms)
         */
        int minSizeForSimCheck = 6;

        /**
         *
         * For fragment of size s we consider fragments of size:
         * min size: min( s-sizeLB_Abs , floor( sizeLB_Rel * s ) )
         * max size: max( s+sizeUB_Abs , ceil( sizeUB_Rel * s ) )
         *
         */
        //
        int sizeLB_Abs = 4;
        //int sizeUB_Abs = 4;
        double sizeLB_Rel = 0.75;
        double sizeUB_Rel = 1.33;

        int consideredConnectorDistances[][] =
                new int[][] { { 1,4 } , // dist 0, irrelevant
                              { 1,4 } , // dist 1, (probably irrelevant)
                        { 1,4 } , // dist 2
                        { 2,5 } , // dist 3
                        { 2,6 } , // dist 4
                        { 3,7 } , // dist 5
                        { 4,9 } , // dist 6
                        { 5,10 } , // dist 7
                        { 6,11 } , // dist 8
                        { 7,12 } , // dist 9
                        { 7,13 } , // dist 10
                        { 7,14 } , // dist 11
                        { 8,16 } , // dist 12
                        { 9,17 } , // dist 13
                        { 10,18 } , // dist 14
                        { 11,19 } , // dist 15
                        { 12,20 } , // dist 16
                        { 13,21 } , // dist 17
                        { 14,22 } , // dist 18
                        { 15,23 } , // dist 19
                        { 16,24 } , // dist 21
                        { 17,25 } , // dist 22
                        { 18,26 } , // dist 23
                        { 19,27 } , // dist 24
                        { 20,28 } , // dist 25
                        { 21,29 } , // dist 26
                        { 22,30 } , // dist 27
                        { 23,31 } , // dist 28
                        { 24,32 } , // dist 29
                        { 25,33 } , // dist 30
                        { 26,34 } , // dist 31
                        { 27,35 }   // dist 32 and greater
                        //
        };


        public List<TopoInfo> createCompatibleTopoConstraints(TopoInfo topo) {
            int s = topo.size;
            int min_size = (int) Math.min( s-sizeLB_Abs , Math.floor( sizeLB_Rel * s ) );
            int max_size = (int) Math.max( s+sizeLB_Abs , Math.ceil(  sizeUB_Rel * s ) );

            if( topo.connDistances.size() == 0 ) {
                List<TopoInfo> compat = new ArrayList<>();
                for(int si=min_size;si<=max_size;si++) {compat.add(new TopoInfo(topo.connectors,topo.connDistances,si));}
                return compat;
            }
            else if( topo.connDistances.size() == 1 ) {
                List<TopoInfo> compat = new ArrayList<>();
                int topd_a = Math.min(topo.connDistances.get(0),consideredConnectorDistances.length-1);
                for(int da = consideredConnectorDistances[topd_a][0] ; da <= consideredConnectorDistances[topd_a][1] ; da++ ) {
                    ArrayList<Integer> conn_dists_a = new ArrayList<>();
                    conn_dists_a.add(da);
                    for(int si=min_size;si<=max_size;si++) {compat.add(new TopoInfo(topo.connectors,conn_dists_a,si));}
                }
                return compat;
            }
            else if( topo.connDistances.size() == 2) {
                List<TopoInfo> compat = new ArrayList<>();
                int topd_a = Math.min(topo.connDistances.get(0),consideredConnectorDistances.length-1);
                int topd_b = Math.min(topo.connDistances.get(1),consideredConnectorDistances.length-1);
                for(int da = consideredConnectorDistances[topd_a][0] ; da <= consideredConnectorDistances[topd_a][1] ; da++ ) {
                    ArrayList<Integer> conn_dists_a = new ArrayList<>();
                    conn_dists_a.add(da);
                    for(int db = consideredConnectorDistances[topd_b][0] ; db <= consideredConnectorDistances[topd_b][1] ; db++ ) {
                        ArrayList<Integer> conn_dists_b = new ArrayList<>(conn_dists_a);
                        conn_dists_b.add(db);
                        for(int si=min_size;si<=max_size;si++) {compat.add(new TopoInfo(topo.connectors,conn_dists_b,si));}
                    }
                }
                return compat;
            }
            else if(topo.connDistances.size() == 3) {
                int topd_a = Math.min(topo.connDistances.get(0),consideredConnectorDistances.length-1);
                int topd_b = Math.min(topo.connDistances.get(1),consideredConnectorDistances.length-1);
                int topd_c = Math.min(topo.connDistances.get(2),consideredConnectorDistances.length-1);

                List<TopoInfo> compat = new ArrayList<>();
                for (int da = consideredConnectorDistances[topd_a][0]; da <= consideredConnectorDistances[topd_a][1]; da++) {
                    ArrayList<Integer> conn_dists_a = new ArrayList<>();
                    conn_dists_a.add(da);
                    for (int db = consideredConnectorDistances[topd_b][0]; db <= consideredConnectorDistances[topd_b][1]; db++) {
                        ArrayList<Integer> conn_dists_b = new ArrayList<>(conn_dists_a);
                        conn_dists_b.add(db);
                        for (int dc = consideredConnectorDistances[topd_c][0]; dc <= consideredConnectorDistances[topd_c][1]; dc++)
                        {
                            ArrayList<Integer> conn_dists_c = new ArrayList<>(conn_dists_b);
                            conn_dists_c.add(dc);
                            for (int si = min_size; si <= max_size; si++) {
                                compat.add(new TopoInfo(topo.connectors, conn_dists_c, si));
                            }
                        }
                    }
                }
                return compat;
            }

            // if we encounter case that we do not yet support
            return null;
        }

    }




    /**
     * connectors is a bitset with bit 0 for U, 1 for Np, and so on
     * connDistances is
     */
    public static class TopoInfo implements Comparable<TopoInfo> , Serializable {
        public final BitSet connectors;
        public final List<Integer> connDistances;
        public final int size;
        private final int hash;
        public TopoInfo(BitSet connis, List<Integer> connDistances, int size) {
            this.connectors= (BitSet) connis.clone();
            this.connDistances = new ArrayList<>(connDistances);
            this.size = size;

            hash = connectors.hashCode() ^ connDistances.hashCode() ^ Integer.hashCode(size);
        }
        public boolean equals(Object o) {
            if( ! (o instanceof TopoInfo) ) {return false;}
            TopoInfo to = (TopoInfo) o;
            if (to.size==this.size) {
                if( this.connectors.equals(to.connectors) ) {
                    if( this.connDistances.equals(to.connDistances) ) {
                        return true;
                    }
                }
            }
            return false;
        }
        public int hashCode() {
            return this.hash;
        }

        @Override
        public int compareTo(TopoInfo topoInfo) {
            int connis_e = compareBitSets(this.connectors,topoInfo.connectors);
            if(connis_e!=0) {return connis_e;}
            int nba = 0;
            int nbb = 0;
            for(int za=0;za<this.connectors.cardinality();za++) {
                nba = this.connectors.nextSetBit(nba);
                nbb = topoInfo.connectors.nextSetBit(nbb);
                if(nba != nbb) {
                    return Integer.compare(nba,nbb);
                }
                nba = nba+1;
                nbb = nbb+1;
            }
            return Integer.compare(this.size,topoInfo.size);
        }

        public String toString() {
            List<Integer> ti = this.connectors.stream().boxed().collect(Collectors.toList());
            return ti.toString()+"__"+this.connDistances.toString()+"__"+this.size;
        }
    }

    public static int compareBitSets(BitSet lhs, BitSet rhs) {
            if (lhs.equals(rhs)) return 0;
            BitSet xor = (BitSet)lhs.clone();
            xor.xor(rhs);
            int firstDifferent = xor.length()-1;
            if(firstDifferent==-1)
                return 0;

            return rhs.get(firstDifferent) ? 1 : -1;
    }

    public static TopoInfo computeTopoInfoForFragment(SynthonSpace.FragId fid) {
        StereoMolecule mi = HyperspaceUtils.parseIDCode(fid.idcode);
        return computeTopoInfoForMolecule(mi);
    }

    public static TopoInfo computeTopoInfoForMolecule(StereoMolecule mi) {

        //StereoMolecule mi = HyperspaceUtils.parseIDCode(fid.idcode);
        BitSet connis     = computeConnectorBitSet(mi);
        int    size       = mi.getAtoms();

        // loop over canonical conn distances and compute distances:
        List<Integer> conn_dists = new ArrayList<>();
        if(connis.cardinality()==1) {
            // no distances :)
        }
        else {
            List<Integer> connis_i = connis.stream().boxed().collect(Collectors.toList());
            for(int za=0;za<connis_i.size()-1;za++) {
                for(int zb=za+1;zb<connis_i.size();zb++) {
                    int ca = connis_i.get(za);
                    int cb = connis_i.get(zb);
                    int atom_pos_a = -1;
                    int atom_pos_b = -1;
                    for(int zi=0;zi<mi.getAtoms();zi++) {
                        if(mi.getAtomicNo(zi)== (ca+92) ) {
                            atom_pos_a = zi;
                            break;
                        }
                    }
                    for(int zi=0;zi<mi.getAtoms();zi++) {
                        if(mi.getAtomicNo(zi)== (cb+92) ) {
                            atom_pos_b = zi;
                            break;
                        }
                    }
                    conn_dists.add( mi.getPathLength(atom_pos_a,atom_pos_b) );
                }
            }
        }
        return new TopoInfo(connis,conn_dists,size);
    }

//    private static List<Integer> computeCanonicalConnectorDistancesForFragment() {
//        List<Integer> distances = new ArrayList<>();
//        return null;
//    }



    public static void runFragmentAnalysis(SynthonSpace ss) {
        Map<TopoInfo,List<FragId>> fragments_sorted = new HashMap<>();
        Map<TopoInfo,Set<String>> fragment_structures_sorted = new HashMap<>();


        for(BitSet ki : ss.fragments_by_connectors.keySet()) {
            for( BitSet mbi : ss.fragments_by_connectors.get(ki).keySet() ) {
                for( FragId fid : ss.fragments_by_connectors.get(ki).get(mbi)) {
                    TopoInfo tinf = computeTopoInfoForFragment(fid);
                    if(!fragments_sorted.containsKey(tinf)) {fragments_sorted.put(tinf,new ArrayList<>());}
                    fragments_sorted.get(tinf).add(fid);
                    if(!fragment_structures_sorted.containsKey(tinf)) {fragment_structures_sorted.put(tinf,new HashSet<>());}
                    fragment_structures_sorted.get(tinf).add(fid.idcode);

                }
            }
        }

        Map<String,Set<TopoInfo>> map_reactand_to_compatible_topos = new HashMap<>();
        for(String ri : ss.fragment_map.keySet()) {
            for( Integer ki : ss.fragment_map.get(ri).keySet() ) {
                String rk_id = ri+"_"+ki;
                Set<TopoInfo> tinfis = new HashSet<>();
                for(FragId fid : ss.fragment_map.get(ri).get(ki)) {
                    tinfis.add(computeTopoInfoForFragment(fid));
                }
                System.out.println(rk_id+" --> " + tinfis.size());
            }
        }


        System.out.println("mkay");
        List<TopoInfo> sorted = new ArrayList<>(fragments_sorted.keySet());

        //sorted.sort( (x,y) -> x.compareTo(y) );
        sorted.sort( (x,y) -> - Integer.compare( fragments_sorted.get(x).size() , fragments_sorted.get(y).size() ) );


        for(TopoInfo tinfi : sorted) {
            System.out.println("num frags: "+fragments_sorted.get(tinfi).size()+" num structures: " + fragment_structures_sorted.get(tinfi).size() +" --> "+tinfi.toString());
        }

    }

    public static void main(String args[]) {
        SynthonSpace sspace = null;
        try {
            sspace = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
            //sspace = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/real_space_latest_SMALL_FragFp.data");
        } catch (IOException e) {
            e.printStackTrace();
        }
        sspace.initAfterJavaDeserialization();

        runFragmentAnalysis(sspace);
    }










    private static int logLevel_findCandidates = 1;

    // Lock objects for initialization of hash maps:
    transient Object LOCK_mTopos = new Object();
    transient Object LOCK_mFragmentTopos = new Object();
    transient Object LOCK_mFragmentToposInv = new Object();
    transient Object LOCK_mFastSimilaritySearchers = new Object();

    public Map<FragId,TopoInfo>        mTopos         = new HashMap<>();

    // to store available topos for every reactand
    public Map<FragType,Set<TopoInfo>> mFragmentTopos = new HashMap<>();
    public Map<TopoInfo,Set<FragType>> mFragmentToposInv = new HashMap<>();

    //public Map<FragType,Map<BitSet,FastSimilaritySearcher>> mFastSimilaritySearchers = new HashMap<>();
    public Map<FragType,Map<TopoInfo,FastSimilaritySearcher>> mFastSimilaritySearchers = new HashMap<>();

    //public LSHProvider mLSH_AllFragments_withNonUniqueConnectors = null;
    //public Map<BitSet,Set<FragType>> ffps_from_non_unique_connis_to_fragtypes = null;

    /**
     * !! Guards write access and manipulations of:
     *
     * this.mFastSimlaritySearchers
     *
     */
    private Object lock_WriteFSS = new Serializable(){};

    public SynthonSimilaritySpace3(SynthonSpace ss_space) {
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

        this.rxns = ss_space.rxns;
        this.rxns_by_connector_counts = ss_space.rxns_by_connector_counts;
        this.frags_by_id = ss_space.frags_by_id;
        this.frags_sorted_by_connector_fp = ss_space.frags_sorted_by_connector_fp;
        this.connector_fps_sorted_by_fragtype = ss_space.connector_fps_sorted_by_fragtype;

        this.setFP(ss_space.getDescriptorHandler_PerThread(),ss_space.getBits());


    }

    public void initFastSimilaritySearchers(int threads, Consumer<Double> progress) {
        //mANN = LSHOracle.initDefault(num_lsh_functions, this.BITS, projection_bits, new ArrayList<>(this.fragments.keySet()) );
        int num_total_synthons     = (int) this.fragments_by_connectors.values().stream().flatMap( fi -> fi.values().stream() ).mapToInt(si -> si.size()).sum();
        int num_processed_synthons = 0;

        int num_total_synthon_sets = (int) this.fragment_map.values().stream().mapToInt(vi -> vi.keySet().size()).sum();
        int processed_synthon_sets = 0;

        int num_synthon_sets = this.fragment_map.values().stream().mapToInt( vi -> vi.keySet().size() ).sum();
        int num_synthon_sets_processed = 0;

        ExecutorService executor_A = Executors.newFixedThreadPool(threads);

        List<Callable<Integer>> tasks_init = new ArrayList<>();

        mFastSimilaritySearchers = new HashMap<>();
        for (String rxn : this.fragment_map.keySet()) {
            for (FragType ft : this.fragment_type_map.get(rxn).values()) {

                Callable<Integer> task_init = new Callable<Integer>() {
                    @Override
                    public Integer call() throws Exception {

                        synchronized(LOCK_mFastSimilaritySearchers) {
                            if (!mFastSimilaritySearchers.containsKey(ft)) {
                                mFastSimilaritySearchers.put(ft, new HashMap<>());
                            }
                        }

                        System.out.println("Compute Topologies of Synthon Set " + ft.toString());
                        List<FragId> fragments_i = fragment_map.get(ft.rxn_id).get(ft.frag);

                        Map<TopoInfo, List<FragId>> synthons_sorted_by_topo = new HashMap<>();
                        for (int zi = 0; zi < fragments_i.size(); zi++) {
                            //StereoMolecule fmi = HyperspaceUtils.parseIDCode(fragments_i.get(zi).idcode);
                            TopoInfo ti = computeTopoInfoForFragment(fragments_i.get(zi));

                            synchronized(LOCK_mTopos) {
                                mTopos.put(fragments_i.get(zi), ti);
                            }

                            if (!synthons_sorted_by_topo.containsKey(ti)) {
                                synthons_sorted_by_topo.put(ti, new ArrayList<>());
                            }
                            synthons_sorted_by_topo.get(ti).add(fragments_i.get(zi));

                            // init inverse map:
                            synchronized(LOCK_mFragmentToposInv) {
                                if (!mFragmentToposInv.containsKey(ti)) {
                                    mFragmentToposInv.put(ti, new HashSet<>());
                                }
                                mFragmentToposInv.get(ti).add(ft);
                            }
                        }

                        synchronized(LOCK_mFragmentTopos) {
                            mFragmentTopos.put(ft, new HashSet<>(synthons_sorted_by_topo.keySet()));
                        }

                        // create fast similarity searchers (init without fss stuff..)
                        for (TopoInfo ti : synthons_sorted_by_topo.keySet()) {
                            System.out.println("Init FSS for " + ft.toString() + " -> init FSS , Size= " + synthons_sorted_by_topo.get(ti).size() + " Topo " + ti);

                            FastSimilaritySearcher fss = new FastSimilaritySearcher(getDescriptorHandler().getInfo().shortName);
                            // create the input, with idcode and base64 encoded descriptor:
                            List<String[]> fss_input = new ArrayList<>();
                            for (FragId fid : synthons_sorted_by_topo.get(ti)) {
                                String descriptor_b64 = Base64.getEncoder().encodeToString(fid.fp.toByteArray());
                                String idcode_i = fid.idcode;
                                fss_input.add(new String[]{idcode_i, descriptor_b64});
                            }

                            // omit fss stuff and init with zero projection..
                            if (true) {
                                // "list" init
                                fss.initInMemory(fss_input.iterator(), 1, BITS, 0, threads);
                            }

                            synchronized(LOCK_mFastSimilaritySearchers) {
                                mFastSimilaritySearchers.get(ft).put(ti, fss);
                            }
                        }

                        return 1;
                    }
                };

                tasks_init.add( task_init );
            }
        }

        try {
            List<Future<Integer>> results = executor_A.invokeAll(tasks_init);
            for (Future<Integer> fr : results) {
                fr.get();

                num_synthon_sets_processed++;
                progress.accept((1.0*num_synthon_sets_processed) / num_synthon_sets );
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        executor_A.shutdown();

        System.out.println("\nDone with init of topology data!\n");
        num_synthon_sets_processed = num_synthon_sets;
        progress.accept(1.0);
    }



}
