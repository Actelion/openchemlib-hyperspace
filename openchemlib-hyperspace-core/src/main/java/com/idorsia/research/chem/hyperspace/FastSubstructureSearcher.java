package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.MyDescriptorResolver;
import com.idorsia.research.chem.hyperspace.outofmemory.MergeSortExternal;
import com.idorsia.research.chem.hyperspace.outofmemory.OutOfMemoryStringToStringMap;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


/**
 *
 * Some explanation for the out-of-memory option:
 *
 * It first computes the descriptor fingerprints with idcodes and just writes them into
 * a single large file. This file is then used for two things:
 * 1. the input file for the out-of-memory tree generation is created
 * 2. the rows are sorted according to the fingeprints, and then equally sized chunks
 *    of consecutive (numerically sorted) fingerprints are stored into files, such that
 *    the fingerprint to idcode mapping can be resolved.
 *
 *
 * There is one difference between the in-memory and the out-of-memory version: the out-of-memory
 * version supports multiple ids for the same idcode. With the in-memory version, only
 * a single id is possible.
 *
 */

public class FastSubstructureSearcher implements Serializable {


    private int BITS = 512;

    private boolean mOutOfMemory = false;
    //protected int BITS = 1024;

    private String mDescriptorHandlerShortName = DescriptorHandlerLongFFP512.getDefaultInstance().getInfo().shortName;
    transient private DescriptorHandler<long[], StereoMolecule> mFP = DescriptorHandlerLongFFP512.getDefaultInstance();

    transient Map<Thread, DescriptorHandler<long[], StereoMolecule>> mFPs_forThreads = new ConcurrentHashMap<>();

    //private final int BITTREE_BIN_SIZE = 1024;
    private int BITTREE_BIN_SIZE = 64;
    //private int BITTREE_BIN_SIZE_OUT_OF_MEMORY = 512;
    private int BITTREE_BIN_SIZE_OUT_OF_MEMORY = 8192;


    public void setBitTreeBinSize(int binsize) {
        this.BITTREE_BIN_SIZE = binsize;
    }

    public int getBits() {
        return this.BITS;
    }


    BitSetTree mTree = null;


    /**
     * map from idcode to id
     */
    //Map<String,List<String>> mIDs = new ConcurrentHashMap<>();
    Map<String, String> mIDs = new ConcurrentHashMap<>();

    /**
     * map from bitset to idcodes
     */
    Map<BitSet, List<String>> mStructureMap = new HashMap<>();
    transient Object lock_StructureMap = new Object();

    // in case of out-of-memory, we use the following Map to resolve idcodes from bitsets:
    OutOfMemoryStringToStringMap mExternalMapFPtoIDCode = null;


    Map<BitSet, String> mStructureToMap = new HashMap<>();


    public void initWithIDCodes(List<String> idcodes, Map<String, String> names, int threads_for_init, boolean out_of_memory) {
        try {
            initWithIDCodes(idcodes,names,null,threads_for_init,null,out_of_memory);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * Input File must be line-separated, and each line must contain first the idcode, then the fp encoded via the
     * descriptor specific encoding, then the id.
     *
     *
     *
     * @param input_file
     */
    public void initWithIDCodesAndDescriptorsOutOfMemory(String input_file, String directory) throws IOException {
        this.mOutOfMemory = true;

        // if we want the not-in-memory version, try to create the directory for the tree, and check that it is empty..
        // in addition to the tree folder we also need the xmap folder for the fp to idcode mappings
        String folder_tree = null;
        String folder_xmap_fp2idcode   = null;
        String folder_xmap_idcode2name = null;
        String folder_temp = null;
        if(directory!=null) {
            File output_dir = new File(directory);
            if(!output_dir.isDirectory()) {
                // try to create:
                output_dir.mkdir();
            }
            // ok, check that is is empty:
            if(output_dir.listFiles().length!=0) {
                System.out.println("[WARN] Output Dir is not empty!");
            }
            folder_tree = output_dir.getAbsolutePath()+File.separator+"tree_"+(new Random()).nextInt(100000);
            System.out.println("output_dir not empty: try to create subfolder for tree: ");
            File fnew = new File(folder_tree);
            fnew.mkdir();

            folder_xmap_fp2idcode = output_dir.getAbsolutePath()+File.separator+"xmap_"+(new Random()).nextInt(100000);
            File fnew2 = new File(folder_xmap_fp2idcode);
            fnew2.mkdir();

            folder_temp = output_dir.getAbsolutePath()+File.separator+"temp_"+(new Random()).nextInt(100000);
            File fnew3 = new File(folder_temp);
            fnew3.mkdir();
        }

        // three "input" files for the out-of-memory mode: one for descriptors only (to build the tree), the second
        // one to create the fp to idcode mapping, third is idcode to name
        File f_input_fp = null;
        File f_input_fp_and_idcode = null;
        File f_idcode_and_id = null;
        BufferedWriter input_writer_fp = null;
        BufferedWriter input_writer_fp_and_idcode = null;
        BufferedWriter input_writer_idcode_and_id = null;

        f_input_fp            = new File( folder_temp + File.separator + "input_data_fp.txt" );
        f_input_fp_and_idcode = new File( folder_temp + File.separator + "input_data_fp_and_idcode.txt" );
        f_idcode_and_id       = new File( folder_temp + File.separator + "input_idcode_and_id.txt" );
        try {
            input_writer_fp = new BufferedWriter(new FileWriter(f_input_fp));
            input_writer_fp_and_idcode = new BufferedWriter(new FileWriter(f_input_fp_and_idcode));
            input_writer_idcode_and_id = new BufferedWriter(new FileWriter(f_idcode_and_id));
        } catch (IOException e) {
            e.printStackTrace();
        }

        // parse input and create input files:
        BufferedReader in = new BufferedReader(new FileReader(input_file));
        String line = null;
        while( (line = in.readLine()) != null ) {
            String splits[] = line.split("\\s+");
            String descriptor = null;
            try {
                descriptor = Base64.getEncoder().encodeToString(BitSet.valueOf(mFP.decode(splits[1])).toByteArray());
            }
            catch(Exception ex) {
                System.out.println("hmm..");
            }
            input_writer_fp.write(descriptor+"\n");
            input_writer_fp_and_idcode.write(descriptor+"\t"+splits[0]+"\n");
            input_writer_idcode_and_id.write(splits[0]+"\t"+splits[2]+"\n");
        }
        input_writer_fp.flush(); input_writer_fp.close();
        input_writer_fp_and_idcode.flush(); input_writer_fp_and_idcode.close();
        input_writer_idcode_and_id.flush(); input_writer_idcode_and_id.close();



        System.out.println("\nInitialize BitSetTree!\n");

        if(true) {
            try {
                //mTree = BitSetTree.createTreeOutOfMemory(f_input_fp, BITS, BITTREE_BIN_SIZE_OUT_OF_MEMORY, directory);
                mTree = BitSetTree.createTreeOutOfMemory_ZipFile( f_input_fp, BITS,2048 , directory, "test_tree.zip" );
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // now we have to initialize the fp to idcode mapping.
        // in a first step, we sort the file, in a second step, we split the file into smaller chunks to enable
        // binary search for fingerprints.
        //List<String> filenames = MergeSortExternal.splitAndSort(f_input_fp_and_idcode, 50000);
        // We just initialize the external string2multistring map:
        File f_xmap = new File(folder_xmap_fp2idcode);
        OutOfMemoryStringToStringMap xmap = new OutOfMemoryStringToStringMap(f_xmap);
        if(true) {
            try {
                xmap.initFromFile(f_input_fp_and_idcode, f_xmap, 100000);
                //xmap.initFromFile(f_input_fp_and_idcode,f_xmap,100);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        mExternalMapFPtoIDCode = xmap;
    }

    private static int LOG_LEVEL_INIT = 0;

    /**
     * NOTE! The argument descriptors is not required, it can be null, or it can also can
     * contain missing descriptors, then the missing descriptors are computed.
     *
     * @param idcodes
     * @param names
     * @param descriptors
     * @param threads_for_init
     * @param directory
     * @param out_of_memory
     * @throws IOException
     */
    public void initWithIDCodes(List<String> idcodes, Map<String, String> names, Map<String,BitSet> descriptors , int threads_for_init, String directory, boolean out_of_memory) throws IOException {
        this.mOutOfMemory = out_of_memory;

        // if we want the not-in-memory version, try to create the directory for the tree, and check that it is empty..
        // in addition to the tree folder we also need the xmap folder for the fp to idcode mappings
        String folder_temp = null;
        String folder_tree = null;
        String folder_xmap = null;
        if(directory!=null) {
            File output_dir = new File(directory);
            if(!output_dir.isDirectory()) {
                // try to create:
                output_dir.mkdir();
            }
            // ok, check that is is empty:
            if(output_dir.listFiles().length!=0) {
                System.out.println("[WARN] Output Dir is not empty!");
            }
            folder_tree = output_dir.getAbsolutePath()+File.separator+"tree_"+(new Random()).nextInt(100000);
            System.out.println("output_dir not empty: try to create subfolder for tree: ");
            File fnew = new File(folder_tree);
            fnew.mkdir();
            directory = folder_tree;

            folder_xmap = output_dir.getAbsolutePath()+File.separator+"xmap_"+(new Random()).nextInt(100000);
            File fnew2 = new File(folder_xmap);
            fnew2.mkdir();

            folder_temp = output_dir.getAbsolutePath()+File.separator+"temp_"+(new Random()).nextInt(100000);
            File fnew3 = new File(folder_xmap);
            fnew3.mkdir();
        }

        this.mIDs = names;

        // three "input" files for the out-of-memory mode: one for descriptors only (to build the tree), the second
        // one to create the fp to idcode mapping, third is idcode to name
        File f_input_fp = null;
        File f_input_fp_and_idcode = null;
        File f_idcode_and_id = null;
        BufferedWriter input_writer_fp = null;
        BufferedWriter input_writer_fp_and_idcode = null;
        BufferedWriter input_writer_idcode_and_id = null;
        if(out_of_memory) {
            f_input_fp            = new File( folder_temp + File.separator + "input_data_fp.txt" );
            f_input_fp_and_idcode = new File( folder_temp + File.separator + "input_data_fp_and_idcode.txt" );
            f_idcode_and_id       = new File( folder_temp + File.separator + "input_idcode_and_id.txt" );
            try {
                input_writer_fp = new BufferedWriter(new FileWriter(f_input_fp));
                input_writer_fp_and_idcode = new BufferedWriter(new FileWriter(f_input_fp_and_idcode));
                input_writer_idcode_and_id = new BufferedWriter(new FileWriter(f_idcode_and_id));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


        //ExecutorService pool = Executors.newFixedThreadPool(threads_for_init);
        ExecutorService pool = new ThreadPoolExecutor(threads_for_init,threads_for_init,0L,TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<>(4000), new ThreadPoolExecutor.CallerRunsPolicy() );

        // compute descriptors:
        List<Future> tasks = new ArrayList<>();
        for (String ici : idcodes) {

            BufferedWriter finalInput_writer_fp = input_writer_fp;
            BufferedWriter finalInput_writer_fp_and_idcode = input_writer_fp_and_idcode;
            BufferedWriter finalInput_writer_idcode_and_id = input_writer_idcode_and_id;

            final String id_for_ici = names.get(ici);

            Runnable ri = new Runnable() {
                @Override
                public void run() {

                    // compute descriptor:
                    BitSet bi = null;
                    if(descriptors!=null) {
                        bi = descriptors.get(ici);
                    }
                    if(bi==null) {
                        IDCodeParser icp = new IDCodeParser();
                        StereoMolecule mi = new StereoMolecule();
                        icp.parse(mi, ici);
                        bi = getFP(mi);
                    }

                    if(!out_of_memory) {
                        synchronized (lock_StructureMap) {
                            if (!mStructureMap.containsKey(bi)) {
                                mStructureMap.put(bi, new ArrayList<>());
                            }
                            mStructureMap.get(bi).add(ici);
                        }
                    }
                    else {
                        // hmm... for the moment we just create two files:
                        // File A containing just the fingerprints, File B the fingerprints and the idcode
                        // we will later process the second file, to create the necessary files to implement the fingerprrint
                        // to idcode finding mechanism.
                        synchronized (lock_StructureMap) {
                            try {
                                finalInput_writer_fp.write(Base64.getEncoder().encodeToString(bi.toByteArray())+"\n");
                                finalInput_writer_fp_and_idcode.write(Base64.getEncoder().encodeToString(bi.toByteArray())+"\t"+ici+"\n");
                                finalInput_writer_idcode_and_id.write(ici+"\t"+id_for_ici);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }
            };

            Future fi = pool.submit(ri);
            tasks.add(fi);
        }

        // wait for tasks to finish
        pool.shutdown();

        if(LOG_LEVEL_INIT>0){ System.out.println("Wait for computation to finish:");}

        int cnt = 0;
        for (Future fi : tasks) {
            cnt++;
            try {
                fi.get();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
            if (cnt % 1000 == 0) {
                System.out.print(".");
            }
            if (cnt % 40000 == 0) {
                System.out.print("\n");
            }
        }

        if(out_of_memory) {
            input_writer_fp.flush();
            input_writer_fp.close();
            input_writer_fp_and_idcode.flush();
            input_writer_fp_and_idcode.close();
            input_writer_idcode_and_id.flush();
            input_writer_idcode_and_id.close();
        }

        if(LOG_LEVEL_INIT>0){ System.out.println("Done!"); }


        if(LOG_LEVEL_INIT>0){ System.out.println("Initialize BitSetTree!"); }
        if(!out_of_memory) {
            synchronized (lock_StructureMap) { // should not really be needed?..
                mTree = BitSetTree.createTree(new HashSet<>(mStructureMap.keySet()), BITS, BITTREE_BIN_SIZE, directory, null, null);
            }
        }
        else {

            try {
                mTree = BitSetTree.createTreeOutOfMemory( f_input_fp , BITS, BITTREE_BIN_SIZE_OUT_OF_MEMORY, directory);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if(LOG_LEVEL_INIT>0){ System.out.println("Initialize BitSetTree -> Done!"); }

        if(out_of_memory) {
            // now we have to initialize the fp to idcode mapping.
            // in a first step, we sort the file, in a second step, we split the file into smaller chunks to enable
            // binary search for fingerprints.
            //List<String> filenames = MergeSortExternal.splitAndSort(f_input_fp_and_idcode, 50000);
            // We just initialize the external string2multistring map:
            File f_xmap = new File(folder_xmap);
            OutOfMemoryStringToStringMap xmap = new OutOfMemoryStringToStringMap(f_xmap);
            try {
                xmap.initFromFile(f_input_fp_and_idcode,f_xmap,100000);
            } catch (IOException e) {
                e.printStackTrace();
            }
            mExternalMapFPtoIDCode = xmap;
        }
    }

    public static class SubstructurePerformanceMeasurement {
        long time_filter_ms;
        long time_search_ms;
        long hits_filter;
        long hits_ss;
        int  cardinality_fp;

        public SubstructurePerformanceMeasurement(long time_filter_ms_, long time_search_ms_, long hits_filter_, long hits_ss_, int cardinality_fp_) {
            this.time_filter_ms = time_filter_ms_;
            this.time_search_ms = time_search_ms_;
            this.hits_filter    = hits_filter_;
            this.hits_ss        = hits_ss_;
            this.cardinality_fp = cardinality_fp_;
        }
    }

    public List<String> findSubstructure(CachedDescriptorProvider cdh, StereoMolecule m, int max_hits) {
        return findSubstructure(cdh,m,max_hits,null);
    }

    /**
     * Returns the idcodes of the molecules containing the
     * specified substructure.
     *
     * @param cdh
     * @param m
     * @param max_hits
     * @param out_performance
     * @return
     */
    public List<String> findSubstructure(CachedDescriptorProvider cdh, StereoMolecule m, int max_hits, SubstructurePerformanceMeasurement out_performance[]) {

        List<SynthonSpace.ExpandedHit> hits = new ArrayList<>();

        //BitSet fp = this.getFP( getDescriptorHandler_PerThread(), m );
        BitSet fp = cdh.getFP_cached(m);

        IDCodeParser icp = new IDCodeParser();


        List<SynthonSpace.ExpandedHit> hits_2 = new ArrayList<>();
        IDCodeParser icp_2 = new IDCodeParser();

        //String ri = rxn;
        //for(int fi : space.ffps_sorted_by_rxn_and_frag_BT.get(ri).keySet()) {

        boolean found_hit_for_synthonset = false;
        //FragType ft = new FragType(ri,fi);

        //BitSetTree bsti = this.ffps_sorted_by_rxn_and_frag_BT.get(ri).get(fi);
        BitSetTree bsti = mTree;

        BitSetTree.Node n_out[] = new BitSetTree.Node[1];
        bsti.testSubset(fp, n_out);

        long ts_filter_start = System.currentTimeMillis();

        List<BitSet> results = new ArrayList<>();
        bsti.root.collectSuperSets(fp, results);

        long ts_filter_end = System.currentTimeMillis();

        long ts_ss_start = System.currentTimeMillis();

        SSSearcher sss = new SSSearcher();
        if (results.size() > 0) {
            m.setFragment(true);
            sss.setFragment(m);
        }

        List<String> verified_hits = new ArrayList<>();

        // pre-fetch idcodes in case of out-of-memory:
        Map<BitSet,List<String>> prefetched_idcodes = null;
        if(mOutOfMemory) {
            prefetched_idcodes = new HashMap<>();
            List<String> query_strings = new ArrayList<>();
            for(BitSet bsi : results){ query_strings.add(Base64.getEncoder().encodeToString(bsi.toByteArray())); }
            try {
                Map<String,List<String>> prefetched_a = mExternalMapFPtoIDCode.get(query_strings);
                prefetched_idcodes = new HashMap<>();
                for(String si : prefetched_a.keySet()) {
                    prefetched_idcodes.put( BitSet.valueOf(Base64.getDecoder().decode(si)) , prefetched_a.get(si));
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        for (BitSet bsi : results) {
            if(verified_hits.size()>max_hits){break;}

            //List<PFPSpace.FragId> frags = this.fragments.get(bsi).stream().filter(si -> si.rxn_id.equals(ri) && si.frag==fi ).collect(Collectors.toList());
            List<String> frags = new ArrayList<>();

            if(this.mOutOfMemory) {
                // we already pre-fetched the fragments, so we can just grab them from the pre-fetch:
                frags = prefetched_idcodes.get(bsi);
            }
            else {
                frags = mStructureMap.get(bsi);
            }

            for (String fid : frags) {
                if(Thread.currentThread().isInterrupted()) {return new ArrayList<>();}
                StereoMolecule sm = new StereoMolecule();
                icp.parse(sm, fid);
                // subset test..
                sss.setMolecule(sm);
                boolean test = sss.isFragmentInMolecule();
                if (test) {
                    found_hit_for_synthonset = true;
                    String midc = fid;
                    //System.out.println("Hit! -> for "+  mIDs.get(midc) );
                    verified_hits.add(midc);
                    //String mi_id = mIDs.get(midc);
                    //verified_hits.add(midc);
                    //verified_hits.add(mi_id);
                }
            }
//            if (found_hit_for_synthonset) {
//                break;
//            }
        }

        long ts_ss_end = System.currentTimeMillis();
        //}

        if(out_performance!=null && out_performance.length>0) {
            SubstructurePerformanceMeasurement performance =
                    new SubstructurePerformanceMeasurement((ts_filter_end-ts_filter_start),
                            (ts_ss_end-ts_ss_start), results.size(), verified_hits.size(),
                            fp.cardinality());
            out_performance[0] = performance;
        }

        return verified_hits;
    }

    public Map<String,List<String>> resolveIdcodeToIDs(List<String> idcodes) {

        if(this.mOutOfMemory) {
            try {
                return this.mExternalMapFPtoIDCode.get(idcodes);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else {
            Map<String, List<String>> results = new HashMap<>();
            for (String idi : idcodes) {
                List<String> ids = new ArrayList<>();
                ids.add( this.mIDs.get(idi) );
                results.put(idi,ids);
            }
            return results;
        }
        return null;
    }


    public boolean checkIfAllContainSS(StereoMolecule qi) {
        qi.setFragment(true);
        BitSet fpi = this.getFP(qi);

        BitSetTree bsti = mTree;

        BitSetTree.Node n_out[] = new BitSetTree.Node[1];
        bsti.testSubset(fpi, n_out);

        long ts_filter_start = System.currentTimeMillis();

        List<BitSet> results = new ArrayList<>();
        boolean all_ok = bsti.root.checkAllAreSuperset(fpi);
        if(all_ok) {
            System.out.println("SSTest: "+qi.getIDCode());
            //
            IDCodeParser icp = new IDCodeParser();
            SSSearcher sss = new SSSearcher();
            sss.setFragment(qi);
            boolean still_ok = true;
            for (String fid : this.mStructureMap.values().stream().flatMap( vi -> vi.stream() ).collect(Collectors.toList())) {
                StereoMolecule sm = new StereoMolecule();
                icp.parse(sm, fid);
                // subset test..
                sss.setMolecule(sm);
                boolean test = sss.isFragmentInMolecule();
                if (!test) {
                    System.out.println("M: "+sm.getIDCode() +"  F: "+qi.getIDCode());
                    still_ok = false;
                    break;
                }
            }
            return still_ok;
        }
        else {
            return false;
        }
    }

    //private DescriptorHandler<long[], StereoMolecule> mFP = DescriptorHandlerLongFFP1024_plus.getDefaultInstance();


    private boolean descriptorHandlerOnTheFly = false;

    /**
     * If TRUE, then this object does not keep a descriptor handler instance in memory all the
     * time, and the instance is just created when there is a substructure search.
     *
     * @param on_the_fly
     */
    public void setDescriptorHandlerOnTheFly(boolean on_the_fly) {
        this.descriptorHandlerOnTheFly = on_the_fly;
        if(on_the_fly) {
            this.mFPs_forThreads.clear();
            this.mFP = null;
        }
    }

    public boolean isDescriptorHandlerOnTheFly() {
        return this.descriptorHandlerOnTheFly;
    }


    public synchronized DescriptorHandler<long[], StereoMolecule> getDescriptorHandler_PerThread() {
        if(this.descriptorHandlerOnTheFly) {
            return MyDescriptorResolver.resolveDescriptorHandlerFromName(this.mDescriptorHandlerShortName);
        }
        DescriptorHandler<long[], StereoMolecule> dhi = mFPs_forThreads.get(Thread.currentThread());
        if (dhi == null) {
            dhi = mFP.getThreadSafeCopy();
            mFPs_forThreads.put(Thread.currentThread(), dhi);
        }

        return dhi;
    }

    /**
     * Now Threadsafe!
     *
     * @param mi
     * @return
     */
    public BitSet getFP(StereoMolecule mi) {
        return getFP(getDescriptorHandler_PerThread(), mi);
    }


    public static BitSet getFP(DescriptorHandler<long[], StereoMolecule> dh, StereoMolecule mi) {
        BitSet bs = BitSet.valueOf(dh.createDescriptor(mi));
        return bs;
    }

    public void setFP(DescriptorHandler<long[], StereoMolecule> dh, int bits) {
        setFP(dh,bits,false);
    }

    /**
     *
     * @param dh
     * @param bits
     * @param on_the_fly if true, then this objects will not store instances of the descriptor handler,
     *                   but always create them on the fly.
     */
    public void setFP(DescriptorHandler<long[], StereoMolecule> dh, int bits, boolean on_the_fly) {
        this.mDescriptorHandlerShortName = dh.getInfo().shortName;
        if(!on_the_fly) {
            this.mFP = dh;
            this.BITS = bits;
        }
        else {
            this.mFP = null;
            this.BITS = bits;
        }
    }


    public CachedDescriptorProvider createCachedDescriptorProvider() {
//        PFPSpace space = new PFPSpace();
//        space.setFP( this.getDescriptorHandler_PerThread().getThreadSafeCopy() , this.BITS );
        CachedDescriptorProvider cdp = new CachedDescriptorProvider(this.getDescriptorHandler_PerThread().getThreadSafeCopy().getInfo().shortName);
        //CachedDescriptorProvider cdp = new CachedDescriptorProvider(space);

        return cdp;
    }

    public void clearCaches() {
        this.mFPs_forThreads = new HashMap<>();
    }

    /**
     * Creates the reference descriptor object from the stored short name.
     *
     * @param ois
     * @throws Exception
     */
    private void readObject(ObjectInputStream ois) throws Exception {
        ois.defaultReadObject();
        this.mFP = MyDescriptorResolver.resolveDescriptorHandlerFromName(this.mDescriptorHandlerShortName);
        this.mFPs_forThreads = new HashMap<>();
    }


    public static class DefaultStructureDataProvider implements StructureDataProvider {

        List<String> Structures;
        Map<String,String> Descriptors;
        Map<String, String> Names;
        public DefaultStructureDataProvider(List<String> idcodes, Map<String, String> names) {
            this.Structures = idcodes;
            // compute Descriptors:
        }

        @Override
        public boolean hasNext() {
            return false;
        }

        @Override
        public String[] next() {
            return new String[0];
        }
    }

    public static interface StructureDataProvider {
        public boolean hasNext();
        public String[] next();
    }


    public static void main(String args[]) {
        FastSubstructureSearcher fss = new FastSubstructureSearcher();

        try {
            List<String> idcs = new ArrayList<>();
            Map<String,String>  ids = new HashMap<>();
            Map<String,BitSet>  descriptors = new HashMap<>();
            try(FileReader fi_in = new FileReader("/home/liphath1/dev5/git/hyperspace_data/old_1/combined_output.txt")) {
                BufferedReader in = new BufferedReader(fi_in);
                String line = null;
                while( (line=in.readLine()) != null ) {
                    if(idcs.size()>8000000) {
                        break;
                    }
                    String splits[] = line.split("\\s+");
                    idcs.add(splits[0]);
                    ids.put(splits[0],splits[2]);
                    BitSet bi = BitSet.valueOf( DescriptorHandlerLongFFP512.getDefaultInstance().decode(splits[1]) );
                    descriptors.put(splits[0],bi);
                }
            }

            System.out.println("Done load");

            fss.initWithIDCodes(idcs,ids,descriptors,8,null,false);

            try( FileOutputStream out_fss_f = new FileOutputStream("fss.obj") ) {
                ObjectOutputStream out_fss = new ObjectOutputStream(new BufferedOutputStream( out_fss_f ));
                out_fss.writeObject(fss);

                // try to use fss..
                String query = "fbmaA@KFAHDdDfYwgmeV\\diNXHBjZ@@bDhB";
                StereoMolecule qa = HyperspaceUtils.parseIDCode(query);
                CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
                SubstructurePerformanceMeasurement perf[] = new SubstructurePerformanceMeasurement[1];
                List<String> hits = fss.findSubstructure(cdp,qa,2000,perf);
                System.out.println("mkay.. hits= "+hits.size());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        //fss.initWithIDCodes();
        //fss.initWithIDCodes();
    }

    public static void main_1(String args[]) {
        FastSubstructureSearcher fss = new FastSubstructureSearcher();

        try {
            fss.initWithIDCodesAndDescriptorsOutOfMemory("C:\\dev6\\hyperspace\\combined_output.txt","C:\\temp\\fastss_data");
            try( FileOutputStream out_fss_f = new FileOutputStream("fss.obj") ) {
                ObjectOutputStream out_fss = new ObjectOutputStream(new BufferedOutputStream( out_fss_f ));
                out_fss.writeObject(fss);

                // try to use fss..
                String query = "fbmaA@KFAHDdDfYwgmeV\\diNXHBjZ@@bDhB";
                StereoMolecule qa = HyperspaceUtils.parseIDCode(query);
                CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
                SubstructurePerformanceMeasurement perf[] = new SubstructurePerformanceMeasurement[1];
                List<String> hits = fss.findSubstructure(cdp,qa,2000,perf);
                System.out.println("mkay.. hits= "+hits.size());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        //fss.initWithIDCodes();
        //fss.initWithIDCodes();
    }

}
