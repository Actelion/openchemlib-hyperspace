package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


/**
 * Notes about Out-of-memory mode:
 *
 *
 */

import org.apache.commons.io.FileUtils;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

public class LSHProvider implements Serializable {
//    private static final long serialVersionUID = 5245800001065670825L;  // temporarily used for Synple
    private static final long serialVersionUID = -1958435657855323807L;   // // added 17-May-2024, TLS, evidently used historically

    private int logLevel = 0;


    //Map<BitSet,Long> mIDs = new HashMap<>(); // unique ids to identify samples , IDs are from 0 to n-1
    final List<LSHFunction> mHashFunctions;


    List<BitSet> mData = null;

    Map<String,Set<BitSet>> mHashToData = new HashMap<>();
    List<Map<Integer,List<BitSet>>> mHashBuckets = null;




    public LSHProvider(List<LSHFunction> hash_functions ) {
        this.mHashFunctions = new ArrayList<>(hash_functions);
    }
    public List<BitSet> getData() { return this.mData; }

    public void initData(List<BitSet> data) {
        initData(data,1);
    }

    public void initData(List<BitSet> data, int threads) {
        initData_InMemory( new ResettableBitsetListIterator(data) ,threads);
        //initData_OutOfMemory(data.iterator(),threads);
    }

    private boolean  mOutOfMemory = false;
    private String   mOutOfMemory_Directory = null;
    private List<String>                mOutOfMemory_DataDirectories = null;
    private List<Map<Integer,String>>   mOutOfMemory_DataFiles       = null;

    private boolean mZip = false;
    private List<Map<Integer,String>>   mZip_DataFiles       = null;
    private String  mZipFilePath = null;
    transient private ZipFile mZipFile = null;


    public void setOutOfMemory(boolean outOfMemory, String directory, boolean zip, String zipFile) throws IOException {
        this.mOutOfMemory = outOfMemory;
        if(outOfMemory) {
            this.mOutOfMemory_Directory = directory;
            File fi = new File(directory);
            fi.mkdir();
            if(!fi.isDirectory()) {
                throw new IOException("Cannot create or open directory "+directory);
            }
            // create directories for rp function hash bins:
            mOutOfMemory_DataDirectories = new ArrayList<>();
            mOutOfMemory_DataFiles       = new ArrayList<>();
            Random r = new Random();
            String rand_str = ""+r.nextInt()%1000000;
            for(int zi=0;zi<this.mHashFunctions.size();zi++) {
                String dir_i = fi.getAbsolutePath() + File.separator + "rplsh_"+rand_str+"_"+zi;
                File fdi = new File(dir_i);
                fdi.mkdir();
                if(!fdi.isDirectory()) {
                    throw new IOException("Failed to create directory +"+fdi.getAbsolutePath());
                }
                mOutOfMemory_DataDirectories.add(dir_i);
            }
        }
        if(zip) {
            this.mZip = true;
            this.mZipFilePath = zipFile;
        }
    }

    public void initData_OutOfMemory(ResettableBitSetIterator data) throws IOException {
        int batchsize = 1000000;

        this.mOutOfMemory_DataFiles = new ArrayList<>();

        for (int zi = 0; zi < mHashFunctions.size(); zi++) {
            LSHFunction fi = mHashFunctions.get(zi);
            data.reset(); //!! RESET ITERATOR TO START !!

            Map<Integer,String> datafile_names_i = new HashMap<>();

            Map<Integer,List<BitSet>> batch_data = new HashMap<>();
            int cnt_batch = 0;
            while(data.hasNext()) {
                BitSet bj = data.next();
                int hash_ij = fi.hash(bj);
                if(!batch_data.containsKey(hash_ij)) {batch_data.put(hash_ij,new ArrayList<>());}
                batch_data.get(hash_ij).add(bj);
                cnt_batch++;
                if(cnt_batch>=batchsize) {
                    writeBatch(mOutOfMemory_DataDirectories.get(zi),batch_data,datafile_names_i);
                    batch_data = new HashMap<>();
                    cnt_batch = 0;
                }
            }
            // write final batch
            writeBatch(mOutOfMemory_DataDirectories.get(zi),batch_data,datafile_names_i);
            this.mOutOfMemory_DataFiles.add(datafile_names_i);
        }
    }

    /**
     * Creates and writes the batch file.
     *
     *
     * @param directory
     * @param batch
     * @param out_datafile_names !! will be changed, this function will add absolute file names of bucket data files here.
     *
     * @throws IOException
     */
    private static void writeBatch(String directory, Map<Integer,List<BitSet>> batch, Map<Integer,String> out_datafile_names) throws IOException {
        for(Integer ki : batch.keySet()) {
            String filename_ki = createHashBucketFileName(directory,ki);
            out_datafile_names.put(ki,filename_ki);
            File fi = new File( filename_ki );
            BufferedWriter out = new BufferedWriter(new FileWriter(fi,true));
            for(BitSet bi : batch.get(ki)) {
                out.write( Base64.getEncoder().encodeToString(bi.toByteArray()) + "\n");
            }
            out.flush();
            out.close();
        }
    }


    private List<BitSet> loadHashBucket_OutOfMemory(int hash_function, int hash_bucket) throws IOException {

        BufferedReader in = null;
        if(this.mZip) {
            synchronized(this) {
                if( this.mZipFile == null ) {
                    File fi = new File(this.mZipFilePath);
                    try {
                        this.mZipFile = new ZipFile(fi);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
            String filename = this.mZip_DataFiles.get(hash_function).get(hash_bucket);
            BufferedReader bin = new BufferedReader( new InputStreamReader( this.mZipFile.getInputStream( this.mZipFile.getEntry( filename ) ) ) );
            in = bin;
        }
        else {
            String filename = this.mOutOfMemory_DataFiles.get(hash_function).get(hash_bucket);
            in = new BufferedReader(new FileReader(filename));
        }
        String line = null;
        List<BitSet> loaded = new ArrayList<>();
        while( (line=in.readLine()) != null ) {
            if(line.isEmpty()){continue;}
            BitSet bi = BitSet.valueOf(Base64.getDecoder().decode(line));
            loaded.add(bi);
        }
        in.close();
        return loaded;
    }

    /**
     * This requires that the function initData_outOfMemory was called before.
     *
     *
     */
    public void createZipFile() throws IOException {
       if(! (this.mOutOfMemory && this.mZip) ) {
           System.out.println("Cannot create zip file, not configured for this");
           return;
       }

       //this.mZipFile = new ZipFile(this.mZipFilePath);
       ZipOutputStream zip_out = new ZipOutputStream(new BufferedOutputStream( new FileOutputStream( this.mZipFilePath ) ));

       this.mZip_DataFiles = new ArrayList<Map<Integer,String>>();

       for(int zi=0;zi<this.mOutOfMemory_DataFiles.size();zi++) {
           Map<Integer,String> df_i = this.mOutOfMemory_DataFiles.get(zi);
           String zip_path_i = "/data/hf_" + String.format("%07d",zi);
           Map<Integer,String> zip_file_names = new HashMap<>();
           for(Map.Entry<Integer,String> ei : df_i.entrySet()) {
               int hash_i = ei.getKey();
               File f_in_i = new File(ei.getValue());
               String file_data = FileUtils.readFileToString(f_in_i,"UTF-8");
               String hash_path_zip = createHashBucketFileName(zip_path_i,hash_i);
               zip_file_names.put(hash_i,hash_path_zip);
               zip_out.putNextEntry(new ZipEntry(hash_path_zip));
               BufferedWriter out = new BufferedWriter(new OutputStreamWriter(zip_out));
               out.write(file_data);
               out.flush();
           }
           this.mZip_DataFiles.add(zip_file_names);
       }

       zip_out.flush();
       zip_out.close();
    }


    private String getHashBucketFileName(int hash_function, int hash_bucket) {
        return this.mOutOfMemory_DataFiles.get(hash_function).get(hash_bucket);
    }


    private static String createHashBucketFileName(String directory, int hash) {
        return directory + File.separator + hash +".hsb";
    }


    public void initData_InMemory(ResettableBitSetIterator data, int threads) {

//        // init index:
//        this.mIDs = new HashMap<>();
//        long id = 0;
//        for(BitSet sb : data) {
//            this.mIDs.put(sb,id);
//            id++;
//        }

        //mData = new ArrayList<>(data);
        mData = new ArrayList<>();
        while(data.hasNext()) {
            mData.add(data.next());
        }

        List<Map<Integer,List<BitSet>>> hashBuckets = new ArrayList<>();

        // hash
        //int data_hashs[][] = new int[data.size()][mHashFunctions.size()];
        int data_hashs[][] = new int[mData.size()][mHashFunctions.size()];

        ExecutorService executor = Executors.newFixedThreadPool(threads);
//        ExecutorService executor = new ThreadPoolExecutor(
//                threads, threads,
//                0L, TimeUnit.MILLISECONDS,
//                new LinkedBlockingQueue<Runnable>(16),
//                new ThreadPoolExecutor.CallerRunsPolicy()
//        );

        for (int zi = 0; zi < mHashFunctions.size(); zi++) {
            int fzi = zi;
            Runnable ri = new Runnable(){
                @Override
                public void run() {
                    LSHFunction fi = mHashFunctions.get(fzi);
                    if(logLevel>0){System.out.println("Hash " + fzi + "  -> init");}
                    //for (int zj = 0; zj < data.size(); zj++) {
                    for (int zj = 0; zj < mData.size(); zj++) {
                        BitSet bj = mData.get(zj);
                        int hash_ij = fi.hash(bj);
                        data_hashs[zj][fzi] = hash_ij;
                        //zj++;
                    }
                }
            };
            executor.execute(ri);
            // also init the hash bucket
            hashBuckets.add(new HashMap<>());
        }

        executor.shutdown();
        try {
            executor.awaitTermination(10, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        //System.gc();

        mHashBuckets = hashBuckets;

        // insert hashs:
        List<int[]> unique_hashs = new ArrayList<>();
        mHashToData = new HashMap<>();
        for (int zj = 0; zj < mData.size(); zj++) {
            for (int zi = 0; zi < mHashFunctions.size(); zi++) {
                if (!hashBuckets.get(zi).containsKey(data_hashs[zj][zi])) {
                    hashBuckets.get(zi).put(data_hashs[zj][zi], new ArrayList<>());
                }
                hashBuckets.get(zi).get(data_hashs[zj][zi]).add(mData.get(zj));
            }

            if(false) { // this is actually not needed and was just for debugging
                String hc = getHashString(data_hashs[zj]);
                if (!mHashToData.containsKey(hc)) {
                    mHashToData.put(hc, new HashSet<>());
                    unique_hashs.add(data_hashs[zj]);
                }
                mHashToData.get(hc).add(mData.get(zj));
            }
        }


        if(logLevel>1){ System.out.println("ok"); }

    }





    public void findNearestNeighbors(Random r, BitSet x, int max_dist, int max_results, List<BitSet> results) {
        findNearestNeighbors(r,x,max_dist,max_results,results,false);
    }


    public static class Hit implements Comparable<Hit>{
        public final BitSet bs;
        public final int    hd;

        public Hit(BitSet bs, int hd) {
            this.bs = bs;
            this.hd = hd;
        }

        @Override
        public int compareTo(Hit o) {
            return Integer.compare( hd , o.hd ) ;
        }


    }

    public static class Hit2 implements Comparable<Hit2>{
        public final BitSet bs;
        public final int    hd;
        public final int    hash;
        public final double tanimoto_sim;

        public Hit2(BitSet bs, int hd) {
            this.bs = bs;
            this.hd = hd;
            this.hash = bs.hashCode();
            this.tanimoto_sim = Double.NaN;
        }
        public Hit2(BitSet bs, double tanimoto_sim) {
            this.bs = bs;
            this.hd = -1;
            this.hash = bs.hashCode();
            this.tanimoto_sim = tanimoto_sim;
        }

        @Override
        public int compareTo(Hit2 o) {
            int c1 = Integer.compare( hd , o.hd );
            if(c1!=0){return c1;}
            //lexicographic comparison:
            BitSet xor = (BitSet) this.bs.clone();
            xor.xor(o.bs);
            int first_different = xor.length()-1;
            if(first_different==-1){return 0;}
            return o.bs.get(first_different) ? 1 : -1;
        }

        public boolean equals(Object o) {
            if(!(o instanceof Hit2)) {return false;}
            return this.bs.equals( ((Hit2)o).bs);
        }
        public int hashCode() {
            return this.hash;
        }
    }


    public void findAllNearestNeighbors_OutOfMemory( List<BitSet> queries, List<Integer> max_dist, int max_results, Map<BitSet,List<Hit2>> results) throws IOException {
        // determine all hashs of the query bitsets:
        //int hashs[][] = new int[x.size()][this.mHashFunctions.size()];

        // maps from hash function to map from hash bucket to query bitset ids in this hash bucket
        Map<Integer,Map<Integer,List<Integer>>> required_hash_buckets = new HashMap<>();
        for(int zi=0;zi<queries.size();zi++) {
            //Map<Integer,List<Integer>> hash_buckets_i = new HashMap<>();
            for(int zj=0;zj<this.mHashFunctions.size();zj++) {
                BitSet xi = queries.get(zi);
                int hash_ij = this.mHashFunctions.get(zj).hash(xi);
                //if(!hash_buckets_i.containsKey(hash_ij)){hash_buckets_i.put(hash_ij,new ArrayList<>());}
                //hash_buckets_i.get(hash_ij).add(zi);
                if(!required_hash_buckets.containsKey(zj)) { required_hash_buckets.put(zj,new HashMap<>()); }
                if(!required_hash_buckets.get(zj).containsKey(hash_ij)){ required_hash_buckets.get(zj).put(hash_ij,new ArrayList<>()); }
                required_hash_buckets.get(zj).get(hash_ij).add(zi);
            }
            //required_hash_buckets.put(zi,hash_buckets_i);
        }

        // loop over hash bucket datafiles and fetch all required data:
        // maps from query mol id to candidate bitsets
        //Map<Integer,List<BitSet>> candidates = new HashMap<>();
        Map<Integer,Set<BitSet>> candidates = new HashMap<>();
        for(int fi : required_hash_buckets.keySet()) {
            for(int hi : required_hash_buckets.get(fi).keySet()) {
                // load data:
                List<BitSet> data_ki_hi = loadHashBucket_OutOfMemory(fi,hi);
                for(int qi : required_hash_buckets.get(fi).get(hi)) {
                    if( !(candidates.containsKey(qi) )) {
                        //candidates.put(qi,new ArrayList<>());
                        candidates.put(qi,new HashSet<>());
                    }
                    candidates.get(qi).addAll(data_ki_hi);
                }
            }
        }

        //Process bitsets
        for(int qi = 0; qi < queries.size() ; qi++) {
            List<Hit2> best_hits = new ArrayList<>();
            BitSet xi = queries.get(qi);
            int max_dist_i = max_dist.get(qi);
            for (BitSet bsi : candidates.get(qi)) {
                if (bsi == xi) {
                    continue;
                }
                int md = hamming(bsi, xi);
                if (md <= max_dist_i) {
                    best_hits.add(new Hit2(bsi, md));
                }
                //if(md<=max_dist){ results_set.add(bsi); }
            }
            Collections.sort(best_hits);
            Iterator<Hit2> it = best_hits.iterator();
            List<Hit2> results_i = new ArrayList<>();
            for(int zi=0; zi < max_results; zi++) {
                if(it.hasNext()){
                    results_i.add(it.next());
                }
            }
            results.put(queries.get(qi),results_i);
        }
    }

    /**
     * Returns the best max_result hits, sorted by
     * similarity (best first).
     *
     * @param x
     * @param max_results
     * @param results
     */
    public void findAllNearestNeighbors( BitSet x, int max_results, List<BitSet> results ) {

        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }

        // rand order
        //List<Integer> perm = new ArrayList<>(); for(int zi=0;zi<mHashFunctions.size();zi++){ perm.add(zi); }

        PriorityQueue<Hit> best_hits = new PriorityQueue<>(max_results);

        for(int zi=0;zi<mHashFunctions.size();zi++) {

            // find bucket and search:
            int x_hash = this.mHashFunctions.get(zi).hash(x);
            List<BitSet> bshb = mHashBuckets.get(zi).get(x_hash);

            if(bshb==null){continue;}

            for(BitSet bsi : bshb) {
                if(bsi==x){continue;}
                int md = hamming(bsi,x);
                best_hits.add(new Hit(bsi,md));
                //if(md<=max_dist){ results_set.add(bsi); }
            }
        }

        for(int zi=0; zi < max_results; zi++) {
            if(best_hits.isEmpty()) {
                break;
            }
            if(zi==0) {
                // debug..
                Hit hi = best_hits.poll();
                //System.out.println("best: hd:"+hi.hd);
                results.add(hi.bs);
            }
            else {
                results.add(best_hits.poll().bs);
            }
        }
    }

    /**
     * Returns the best max_result hits, sorted by
     * similarity (best first).
     *
     * @param x
     * @param max_results
     * @param results
     */
    public void findAllNearestNeighbors( BitSet x, int max_dist, int max_results, List<Hit2> results ) {

        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }

        // rand order
        //List<Integer> perm = new ArrayList<>(); for(int zi=0;zi<mHashFunctions.size();zi++){ perm.add(zi); }

        Set<Hit2> best_hits = new HashSet<>();

        for(int zi=0;zi<mHashFunctions.size();zi++) {

            // find bucket and search:
            int x_hash = this.mHashFunctions.get(zi).hash(x);
            List<BitSet> bshb = mHashBuckets.get(zi).get(x_hash);

            if(bshb==null){continue;}

            for(BitSet bsi : bshb) {
                if(bsi==x){continue;}
                int md = hamming(bsi,x);
                if(md<=max_dist) {
                    best_hits.add(new Hit2(bsi, md));
                }
                //if(md<=max_dist){ results_set.add(bsi); }
            }
        }

        List<Hit2> best_hits2 = new ArrayList<>(best_hits);
        Collections.sort(best_hits2);
        Iterator<Hit2> it = best_hits2.iterator();
        for(int zi=0; zi < max_results; zi++) {
            if(it.hasNext()){
                results.add(it.next());
            }
        }
    }


    /**
     * NOTE: If the created hash functions satisfy the MIH (multi-index hashing)
     * criteria (i.e. they all cover disjoint bits, but in total cover the complete
     * bitstring), then this method gives the exact nn results.
     *
     * @param x
     * @param max_dist
     * @return
     */
    public List<Hit2> findAllNearestNeighborsMIH(BitSet x, int max_dist) {

        // mih candidate distance:
        //int mih_dist = Math.floorDiv( max_dist , mHashFunctions.size() );
        int mih_dist = (int) Math.ceil( (1.0*max_dist) / mHashFunctions.size() );

        Set<BitSet> candidates = new HashSet<>();

        for(int zi=0;zi<mHashFunctions.size();zi++) {

            int x_hash = this.mHashFunctions.get(zi).hash(x);

            // consider buckets within mih radius:
            for(int bi = x_hash-mih_dist ; bi < x_hash+mih_dist+1 ; bi++) {
                // add candidate:
                candidates.addAll( mHashBuckets.get(zi).get(bi) );
            }
        }

        // add results:
        List<Hit2> hits = new ArrayList<>();
        for(BitSet bsi : candidates) {
            int md = hamming(x,bsi);
            if(md<=max_dist) {
                hits.add(new Hit2(bsi, md));
            }
        }

        return hits;
    }

    /**
     *
     *
     * @param r
     * @param x
     * @param max_dist
     * @param max_results
     * @param results
     * @param shuffle_buckets with shuffled buckets it is better possible to sample random hits
     */
    public void findNearestNeighbors(Random r, BitSet x, int max_dist, int max_results, List<BitSet> results, boolean shuffle_buckets) {

        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }


        Set<BitSet> results_set = new HashSet<>();

        // rand order
        List<Integer> perm = new ArrayList<>(); for(int zi=0;zi<mHashFunctions.size();zi++){ perm.add(zi); }
        Collections.shuffle(perm,r);

        for(int zi=0;zi<mHashFunctions.size();zi++) {

            // find bucket and search:
            int x_hash = this.mHashFunctions.get(perm.get(zi)).hash(x);
            List<BitSet> bshb = mHashBuckets.get(perm.get(zi)).get(x_hash);

            if(bshb==null){continue;}

            if(shuffle_buckets) { Collections.shuffle(bshb); }

            for(BitSet bsi : bshb) {
                if(bsi==x){continue;}
                int md = hamming(bsi,x);
                if(md<=max_dist){ results_set.add(bsi); }
            }
        }
        results.addAll(results_set);
    }

    public void exactFindNearestNeighbors(BitSet x, int max_dist,int max_results,List<BitSet> results_exact) {
        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }
        for(BitSet bsi : mData) {
            int hdi = LSHProvider.hamming(x,bsi);
            if(x==bsi){continue;}
            if(hdi<=max_dist){results_exact.add(bsi);}
            if(results_exact.size()>max_results){break;}
        }
    }
    public void exactFindNearestNeighbors2(BitSet x, int max_dist,int max_results,List<Hit2> results_exact) {
        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }
        for(BitSet bsi : mData) {
            int hdi = LSHProvider.hamming(x,bsi);
            if(x==bsi){continue;}
            if(hdi<=max_dist){results_exact.add(new Hit2(bsi,hdi));}
            if(results_exact.size()>max_results){break;}
        }
    }
    public List<Hit2> exactFindNearestNeighbors2(BitSet x, int max_dist,int max_results) {
        Set<Hit2> results_exact = new HashSet<>();
        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }
        for(BitSet bsi : mData) {
            int hdi = LSHProvider.hamming(x,bsi);
            if(x==bsi){continue;}
            if(hdi<=max_dist){results_exact.add(new Hit2(bsi,hdi));}
            if(results_exact.size()>max_results){break;}
        }
        List<Hit2> results_exact_list = new ArrayList<>(results_exact);
        return results_exact_list;
    }
    public List<Hit2> exactFindKNearestNeighbors2_Tanimoto(BitSet q, int k) {
        PriorityQueue<Hit2> results_exact = new PriorityQueue<>( (x, y) -> Double.compare( x.tanimoto_sim , y.tanimoto_sim ) );
        if(this.mOutOfMemory) {
            throw new Error("Only available for in-memory LSH");
        }
        for(BitSet bsi : mData) {
            //int hdi = LSHProvider.hamming(x,bsi);
            //if(x==bsi){continue;}
            double tani = tanimoto_similarity(q,bsi);
            results_exact.add( new Hit2(bsi,tani));
            if(results_exact.size()>k) { results_exact.poll(); }
            //if(results_exact.size()>max_results){break;}
        }
        List<Hit2> results_exact_list = new ArrayList<>(results_exact);
        Collections.reverse(results_exact_list);
        return results_exact_list;
    }

    /**
     * Counts the number of buckets with only a single hashed entry.
     *
     * @return
     */
    public long findSingletons() {
        return this.mHashBuckets.stream().flatMap( hi -> hi.values().stream() ).mapToInt( li -> li.size() ).filter( si -> si==1 ).count();
    }


    public static LSHProvider initDefault(List<BitSet> data) {
        //LSHOracle oracle = initDefault(32,data.get(0).size(),32);
        //LSHOracle oracle = initDefault(8,data.get(0).size(),32);
        //LSHOracle oracle = initDefault(16,data.get(0).size(),128);

        //LSHOracle oracle = initDefault(8,data.get(0).size(),32); // this is the real one..
        LSHProvider oracle = initDefault(8,data.get(0).size(),128, data, 1); // this is the real one..
        return oracle;
    }



    /**
     * If "zip" is true, then path indicates the path for the
     * temp directory, and zipFilePath indicates the path of
     * the final zip file.
     *
     */
    public static class LSHProviderConfig {
        public final String path;
        public final boolean outOfMemory;
        public final boolean zip;
        public final String zipFilePath;
        public final int init_threads;
        public LSHProviderConfig(boolean outOfMemory, String path, int init_threads) {
            this(outOfMemory,path,false,null,init_threads);
        }
        public LSHProviderConfig(boolean outOfMemory, String path, boolean zip, String zipFilePath, int init_threads) {
            this.outOfMemory = outOfMemory;
            this.path = path;
            this.init_threads = init_threads;
            this.zip  = zip;
            this.zipFilePath = zipFilePath;
        }
    }

    /**
     *
     * A good value for num_hash_functions is maybe 32. Generally, the more the better (only drawback is longer init time).
     *
     * A good value for projection_bits seems to be 48 (for most bit descriptors with length in the range 512-2048 ).
     * If the value is too large, we will not find any nearest neighbors.
     * If the value is too small, we will not get any performance improvements.
     *
     *
     *
     * @param num_hash_functions the number of hash functions
     * @param data_bits the size of your bitsets
     * @param projection_bits the number of projected bits per hash function
     *
     * @return
     */
    public static LSHProvider initDefault(int num_hash_functions, int data_bits, int projection_bits, List<BitSet> data, int initThreads) {
        LSHProviderConfig conf = new LSHProviderConfig(false,null,initThreads);
        return initDefault(num_hash_functions,data_bits,projection_bits, new ResettableBitsetListIterator(data) ,conf);
    }
//        //int M = 2048; // we could also set this higher (or just use the BitSet default hash function?)
//        int M = 44719;
//
//        // create lsh oracle
//        Random r2 = new Random();
//        List<LSHProvider.LSHFunction> lsh = new ArrayList<>();
//        for(int zi=0;zi<num_hash_functions;zi++) {
//            LSHProvider.RPLSHFunction fi = LSHProvider.RPLSHFunction.createRandom(M,r2,data_bits,projection_bits);
//            lsh.add(fi);
//        }
//
//        LSHProvider oracle = new LSHProvider(lsh);
//        oracle.initData_InMemory(data.iterator(),initThreads);
//        return oracle;
//    }


    /**
     * If this init function is used, then a call to findNearestNeighborsMIH is guaranteed to
     * return all hits.
     *
     * NOTE! It is required that data_bits is dividable by mih_m! Otherwise we might be sightly
     * off with the MIH guarantees.
     *
     * @param mih_m the number of multi-hash indeces
     * @param data_bits
     * @param data
     * @param conf
     * @return
     */
    public static LSHProvider initDefault_MIH( int mih_m, int data_bits, ResettableBitSetIterator data, LSHProviderConfig conf) {

        // probably we should somehow figure out statistics for bits, to create balanced bins..
        // for the moment we shuffle (e.g. w.r.t. FragFp..)
        Random ri = new Random(456);
        List<Integer> perm = new ArrayList<>();
        for(int zi=0;zi<data_bits;zi++) {perm.add(zi);}
        Collections.shuffle(perm,ri);


        List<BitSet> mih_projections = new ArrayList<>();

        // distribute equally:
        int numBuckets = mih_m;
        int numItems = data_bits;
        int itemsPerBucket = (numItems / numBuckets);
        int remainingItems = (numItems % numBuckets);
        List<Integer> bucket_sizes = new ArrayList<>();
        for (int i = 1; i <= numBuckets; i++)
        {
            int extra = (i <= remainingItems) ? 1:0;
            bucket_sizes.add((itemsPerBucket + extra));
            //System.out.println("bucket " + i + " contains " + (itemsPerBucket + extra) + " items.");
        }

        // create random disjoint projections:
        int perm_position = 0;
        for(int zi=0;zi<mih_m;zi++) {
            BitSet bsi = new BitSet(data_bits);
            for(int zb=0;zb<bucket_sizes.get(zi);zb++) {
                bsi.set( perm.get(perm_position) );
                perm_position++;
            }
            mih_projections.add(bsi);
        }

        return initDefault(data_bits,data,conf,mih_projections);
    }

        /**
         *
         * @return
         */
    public static LSHProvider initDefault( int data_bits, ResettableBitSetIterator data, LSHProviderConfig conf , List<BitSet> projections) {

        // create lsh oracle
        Random r2 = new Random();
        List<LSHProvider.LSHFunction> lsh = new ArrayList<>();
        for(int zi=0;zi<projections.size();zi++) {
            LSHProvider.RPLSHFunction2 fi = new LSHProvider.RPLSHFunction2(projections.get(zi));
            lsh.add(fi);
        }
        if(conf.outOfMemory) {


            LSHProvider oracle = new LSHProvider(lsh);
            try {
                if(!conf.zip) {
                    oracle.setOutOfMemory(true, conf.path, false, null);
                }
                else {
                    oracle.setOutOfMemory(true, conf.path, true, conf.zipFilePath);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

            try {
                oracle.initData_OutOfMemory(data);

                if(conf.zip) {
                    // create zip file
                    oracle.createZipFile();
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
            return oracle;
        }
        else {
            LSHProvider oracle = new LSHProvider(lsh);
            oracle.initData_InMemory(data,conf.init_threads);
            return oracle;
        }
    }


    public static LSHProvider initDefault(int num_hash_functions, int data_bits, int projection_bits, ResettableBitSetIterator data, LSHProviderConfig conf) {
        //int M = 2048; // we could also set this higher (or just use the BitSet default hash function?)
        //int M = 44719;
        //int M = 64;
        int M = 32768;


        // create lsh oracle
        Random r2 = new Random();
        List<LSHProvider.LSHFunction> lsh = new ArrayList<>();
        for(int zi=0;zi<num_hash_functions;zi++) {
            LSHProvider.RPLSHFunction fi = LSHProvider.RPLSHFunction.createRandom(M,r2,data_bits,projection_bits);
            lsh.add(fi);
        }

        if(conf.outOfMemory) {


            LSHProvider oracle = new LSHProvider(lsh);
            try {
                if(!conf.zip) {
                    oracle.setOutOfMemory(true, conf.path, false, null);
                }
                else {
                    oracle.setOutOfMemory(true, conf.path, true, conf.zipFilePath);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

            try {
                oracle.initData_OutOfMemory(data);

                if(conf.zip) {
                    // create zip file
                    oracle.createZipFile();
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
            return oracle;
        }
        else {
            LSHProvider oracle = new LSHProvider(lsh);
            oracle.initData_InMemory(data,conf.init_threads);
            return oracle;
        }
    }


    public void test_nn(Random r, int num_tests, int max_dist) {
        for(int zi=0;zi<num_tests;zi++) {
            final BitSet to_find = this.mData.get( r.nextInt(this.mData.size()) );

            List<BitSet> results_ann_list   = new ArrayList<>();
            List<BitSet> results_exact_list = new ArrayList<>();

            long ts_a = System.currentTimeMillis();
            findNearestNeighbors(r,to_find,max_dist,2000,results_ann_list);
            long ts_b = System.currentTimeMillis();
            exactFindNearestNeighbors(to_find,max_dist,2000,results_exact_list);
            long ts_c = System.currentTimeMillis();

            results_ann_list   = new ArrayList<>( results_ann_list.stream().distinct().collect(Collectors.toList()) );
            results_exact_list = new ArrayList<>( results_exact_list.stream().distinct().collect(Collectors.toList()) );

            results_ann_list.sort( (x,y) -> Integer.compare( hamming(to_find,x) , hamming(to_find,y) ) );
            results_exact_list.sort( (x,y) -> Integer.compare( hamming(to_find,x) , hamming(to_find,y) ) );

            System.out.println("Times:  "+ (ts_b-ts_a) + " ms / "+(ts_c-ts_b)+" ms" );
            System.out.println("Results ann:   "+ results_ann_list.stream().map( bi -> ""+(hamming(to_find,bi)) ).reduce( (x,y) -> x+","+y ) );
            System.out.println("Results exact: "+ results_exact_list.stream().map( bi -> ""+(hamming(to_find,bi)) ).reduce( (x,y) -> x+","+y ) );

            int hits_ann   = results_ann_list.size();
            int hits_exact = results_exact_list.size();
            double avg_dist_ann   = results_ann_list.stream().mapToDouble( bi -> (hamming(to_find,bi)) ).average().orElse(Double.NaN);
            double avg_dist_exact = results_exact_list.stream().mapToDouble( bi -> (hamming(to_find,bi)) ).average().orElse(Double.NaN);
            System.out.println("hits: "+hits_ann+" / "+hits_exact + " avg_dist: "+avg_dist_ann+" / "+avg_dist_exact);
            System.out.println("");
        }
    }



//    public List<BitSet> nearestNeighbor(BitSet x) {
//        int[] x_hash = getHash(x);
//        double max_dist[] = new double[this.mHashFunctions.size()];
//        Arrays.fill(max_dist,Double.POSITIVE_INFINITY);
//        int[] nn = mKDT.getRoot().nearest(x_hash,max_dist);
//        return new ArrayList<>(this.mHashToData.get( getHashString(nn) ));
//    }

//    public List<BitSet> proximityQuery(BitSet x, int dist) {
//        int[] x_hash = getHash(x);
//        int[] lb = new int[this.mHashFunctions.size()];
//        int[] ub = new int[this.mHashFunctions.size()];
//
//        for(int zi=0;zi<lb.length;zi++){
//            lb[zi]=x_hash[zi]-dist;
//            ub[zi]=x_hash[zi]+dist;
//        }
//        IntKDTree.Hypercube range = new IntKDTree.Hypercube(lb,ub);
//        IntKDTree.KDSearchResults results = new IntKDTree.KDSearchResults();
//        mKDT.getRoot().range(range,results);
//
//        List<BitSet> bs = new ArrayList<>();
//        for(IntKDTree.N ni : results.nodes) {
//            bs.addAll(this.mHashToData.get( getHashString(ni.point) ));
//        }
//
//        return bs;
//    }



    public int[] getHash(BitSet x) {
        int x_hashs[] = new int[mHashFunctions.size()];
        for (int zi = 0; zi < mHashFunctions.size(); zi++) {
            x_hashs[zi]=mHashFunctions.get(zi).hash(x);
        }
        return x_hashs;
    }

    private String getHashString(int[] hash) {
        StringBuilder sb = new StringBuilder();
        for(int zi=0;zi<hash.length;zi++){
            sb.append(hash[zi]);
            sb.append(',');
        }
        return sb.toString();
    }











    static class BitSetWithLSH implements Comparable<BitSetWithLSH> , Serializable {
        public final BitSet bs;
        public final int x;
        public BitSetWithLSH(BitSet bs, int x) {
            this.bs = bs;
            this.x  = x;
        }
        @Override
        public int compareTo(BitSetWithLSH o) {
            return Integer.compare(this.x,o.x);
        }
    }

    interface LSHFunction {
        int hash(BitSet b);
    }

    public static interface ResettableBitSetIterator extends Iterator<BitSet> {
        /**
         * Resets the iterator into its original state.
         */
        public void reset();
    }

    public static class ResettableBitsetListIterator implements ResettableBitSetIterator {
        ListIterator<BitSet> iterator;
        List<BitSet> list;
        public ResettableBitsetListIterator(List<BitSet> bitset_list) {
            this.list = bitset_list;
            this.reset();
        }
        @Override
        public void reset() {
            this.iterator = list.listIterator();
        }

        @Override
        public boolean hasNext() {
            return this.iterator.hasNext();
        }

        @Override
        public BitSet next() {
            return this.iterator.next();
        }
    }

    static class RPLSHFunction implements LSHFunction , Serializable {
        private static final long serialVersionUID=-5567530831074257650L; // de-facto this one was used before without being declared; TLS 24Nov2024

        public final int M;
        public final BitSet P;
        public final int[] multipliers;
        public RPLSHFunction(int M, BitSet p, int[] multipliers) {
            this.M = M;
            this.P = p;
            this.multipliers = multipliers;
        }
        public static RPLSHFunction createRandom(int M, Random r, int data_bits, int projection_bits) {
            BitSet p = new BitSet(projection_bits);

            if(projection_bits==0) {
                return new RPLSHFunction(M,p,new int[0]);
            }

            List<Integer> perm = new ArrayList<>();
            for(int zi=0;zi<data_bits;zi++) {  perm.add(zi); }
            Collections.shuffle(perm,r);
            int multipliers[] = new int[projection_bits];
            for(int zi=0;zi<projection_bits;zi++) {
                p.set(perm.get(zi));
                multipliers[zi] = r.nextInt(M);
            }
            return new RPLSHFunction(M,p,multipliers);
        }
        @Override
        public int hash(BitSet b) {

            // we support "no hashing"
            if(this.P.cardinality()==0){
                return 0;
            }

            int value = 0;
            int pos = this.P.nextSetBit(0);
            for(int zi=0;zi<this.P.cardinality()-1;zi++) {
                pos = this.P.nextSetBit(pos+1);
                //if(pos<0){break;}
                value+= ((b.get(pos))?1:0) * this.multipliers[zi];
            }
            return value % this.M;
        }
    }

    static class RPLSHFunction2 implements LSHFunction , Serializable {
        public final BitSet P;
        public RPLSHFunction2(BitSet p) {
            this.P = p;
        }

        @Override
        public int hash(BitSet b) {

            // we support "no hashing"
            if(this.P.cardinality()==0){
                return 0;
            }

            BitSet pi = (BitSet) this.P.clone();
            pi.and(b);
            return pi.cardinality();
        }
    }

    public static int hamming(BitSet a, BitSet b) {
        BitSet xored = (BitSet) a.clone();
        xored.xor(b);
        return xored.cardinality();
//        int size = Math.max(a.size(),b.size());
//        int ha = 0;
//        for(int zi=0;zi<size;zi++) {
//            ha += ((a.get(zi)!=b.get(zi))?1:0);
//        }
//        return ha;
    }

    public static double tanimoto_similarity(BitSet a, BitSet b) {
        BitSet all_bits    = (BitSet) a.clone();
        BitSet shared_bits = (BitSet) a.clone();
        all_bits.or(b);
        shared_bits.and(b);
        double a_or_b = all_bits.cardinality();
        double a_and_b = shared_bits.cardinality();
        return a_and_b / a_or_b;
    }
}
