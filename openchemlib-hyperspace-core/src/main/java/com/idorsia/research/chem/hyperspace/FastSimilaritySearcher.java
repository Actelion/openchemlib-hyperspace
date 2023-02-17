package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.idorsia.research.chem.hyperspace.descriptor.MyDescriptorResolver;
import com.idorsia.research.chem.hyperspace.outofmemory.OutOfMemoryStringToStringMap;
import org.apache.commons.io.LineIterator;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;



/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */




public class FastSimilaritySearcher implements Serializable {

    private int BITS = 512;

    private boolean mOutOfMemory = false;
    //protected int BITS = 1024;

    private String mDescriptorHandlerShortName = null;
    transient private DescriptorHandler<long[], StereoMolecule> mFP = null;

    private Map<BitSet,List<String>> mDescriptor2Structures_inMemory = null;
    private OutOfMemoryStringToStringMap mDescriptor2Structures_outOfMemory = null;

    transient Map<Thread, DescriptorHandler<long[], StereoMolecule>> mFPs_forThreads = new ConcurrentHashMap<>();

    private synchronized DescriptorHandler<long[],StereoMolecule> getDH() {
        Thread ti = Thread.currentThread();
        if(mFPs_forThreads.containsKey(ti)) {return mFPs_forThreads.get(ti);}
        mFPs_forThreads.put(ti,mFP.getThreadSafeCopy());
        return mFPs_forThreads.get(ti);
    }

    public int getBits() {
        return this.BITS;
    }

    LSHProvider mANN = null;

    public FastSimilaritySearcher(String descriptorShortName) {
        this.mDescriptorHandlerShortName = descriptorShortName;
        this.mFP = MyDescriptorResolver.resolveDescriptorHandlerFromName(descriptorShortName);
        if(this.mFP==null){throw new Error("Could not resolve DescriptorHandler from Short Name: "+descriptorShortName);}
    }

    public LSHProvider getANN_InMemory() {return this.mANN;}

    /**
     *
     * This function inits two things:
     * 1. an out-of-memory lsh provider using the provided descriptor bitstrings
     * 2. an out-of-memory map that maps from bitstrings to idcodes.
     *
     * @param input_file path to line-separated file, where each line contains
     *                   two tab-separated values: first entry is idcode, second entry is the Base64-encoded descriptor
     *
     */
    public void initOutOfMemory(
                                 String input_file,
                                 int num_hash_functions, int num_bits,
                                 int num_projection_bits,
                                 String path,
                                 int xmap_num_files) throws IOException {
        mOutOfMemory = true;

        BITS = num_bits;
        //mANN = LSHProvider.initDefault(num_hash_functions,num_bits,num_projection_bits);

        // make folders for the ANN and the String2String Map:
        File fi = new File(path);
        if(!fi.exists()) {
            fi.mkdir();
        }
        if(!fi.isDirectory()) {throw new IOException("Path must be directory");}
        Random r = new Random();

        File fi_ann = new File( fi.getAbsolutePath() + File.separator + "ann_"+(r.nextInt(10000)) );
        fi_ann.mkdir();
        File fi_xmap = new File( fi.getAbsolutePath() + File.separator + "xmap_"+(r.nextInt(10000)) );
        fi_xmap.mkdir();

        File fi_tmp = new File( fi.getAbsolutePath() + File.separator + "tmp_"+(r.nextInt(10000)) );
        fi_tmp.mkdir();

        // Open the Input File:
        File f_in = new File(input_file);
        //BufferedReader in_a = new BufferedReader(new FileReader(f_in));
        //BufferedReader in_b = new BufferedReader(new FileReader(f_in));



        // 1. Init ANN
        LSHProvider.LSHProviderConfig conf = new LSHProvider.LSHProviderConfig(true,fi_ann.getAbsolutePath(),4);
        //mANN = LSHProvider.initDefault(num_hash_functions,num_bits,num_projection_bits,new BitsetLineIterator(in_a,1),conf);
        mANN = LSHProvider.initDefault(num_hash_functions,num_bits,num_projection_bits,new BitsetLineIterator(f_in,1),conf);

        // 2. Init XMap:
        OutOfMemoryStringToStringMap stringToStringMap = new OutOfMemoryStringToStringMap(fi_xmap);
        stringToStringMap.initFromFile(f_in,fi_xmap,xmap_num_files);
        mDescriptor2Structures_outOfMemory = stringToStringMap;

    }

    /**
     *
     * @param input arrays first entry is idcode, second entry is the Base64-encoded descriptor
     *
     */
    public void initInMemory(
                             Iterator<String[]> input ,
                             int num_hash_functions, int num_bits,
                             int num_projection_bits,
                             int init_threads) {
        mOutOfMemory = false;
        BITS = num_bits;

        List<BitSet> all_data = new ArrayList<>();
        Map<BitSet,List<String>> structures = new HashMap<>();

        while(input.hasNext()) {
            // 1. add bitsets to all_data
            String input_next[] = input.next();
            //BitSet bsi = BitSet.valueOf(dh.decode(input_next[1]));
            BitSet bsi = BitSet.valueOf(Base64.getDecoder().decode(input_next[1]));
            all_data.add( bsi );
            if(!structures.containsKey(bsi)){structures.put(bsi,new ArrayList<>());}
            structures.get(bsi).add(input_next[0]);
        }
        // 1. init ann
        mANN = LSHProvider.initDefault(num_hash_functions,num_bits,num_projection_bits,all_data,init_threads);
        mDescriptor2Structures_inMemory = structures;
    }

    public Map<StereoMolecule,Map<LSHProvider.Hit2,List<String>>> searchAllSimilar(List<StereoMolecule> q, double min_tanimoto_sim , int max_results ) {

        List<BitSet> q_descriptors = q.parallelStream().map(
                mi -> BitSet.valueOf( this.getDH().createDescriptor(mi) ) ).collect(Collectors.toList());

        Map<BitSet,List<LSHProvider.Hit2>> hits = new ConcurrentHashMap<>();

        if(!mOutOfMemory) {
            q_descriptors.parallelStream().distinct().forEach( bi -> {
                int hamming_dist = (int) Math.ceil( bi.cardinality() * (1.0-min_tanimoto_sim) );
                List<LSHProvider.Hit2> results = new ArrayList<>();
                mANN.findAllNearestNeighbors( bi , hamming_dist , max_results , results );
                hits.put(bi,results);
            } );
        }
        else {
            List<Integer> hamming_dist = q_descriptors.stream().map(bi -> (int) Math.ceil(bi.cardinality() * (1.0 - min_tanimoto_sim))).collect(Collectors.toList());
            Map<BitSet, List<LSHProvider.Hit2>> results = new HashMap<>();
            try {
                mANN.findAllNearestNeighbors_OutOfMemory(q_descriptors, hamming_dist, max_results, results);
                hits.putAll(results);
            }
            catch(IOException ex) {
                ex.printStackTrace();
            }
        }

        // package into result map:
        if(!mOutOfMemory)
        {
            Map<StereoMolecule, Map<LSHProvider.Hit2, List<String>>> results_map = new HashMap<>();
            for (int zi = 0; zi < q.size(); zi++) {
                List<LSHProvider.Hit2> results_i = hits.get(q_descriptors.get(zi));
                Map<LSHProvider.Hit2, List<String>> result_map_i = new HashMap<>();
                for (LSHProvider.Hit2 hi : results_i) {
                    result_map_i.put(hi, mDescriptor2Structures_inMemory.get(hi.bs));
                }
                results_map.put(q.get(zi), result_map_i);
            }
            return results_map;
        }
        else {
            List<String> bitsets_encoded = hits.keySet().parallelStream().
                    map( hi -> Base64.getEncoder().encodeToString(hi.toByteArray()) ).collect(Collectors.toList());

            Map<String,List<String>> fetched_idcodes = null;
            try {
                fetched_idcodes = mDescriptor2Structures_outOfMemory.get(bitsets_encoded);
            } catch (IOException e) {
                e.printStackTrace();
            }

            Map<StereoMolecule, Map<LSHProvider.Hit2, List<String>>> results_map = new HashMap<>();
            for (int zi = 0; zi < q.size(); zi++) {
                List<LSHProvider.Hit2> results_i = hits.get(q_descriptors.get(zi));
                Map<LSHProvider.Hit2, List<String>> result_map_i = new HashMap<>();
                for (LSHProvider.Hit2 hi : results_i) {
                    result_map_i.put(hi, fetched_idcodes.get( Base64.getEncoder().encodeToString(hi.bs.toByteArray())));
                }
                results_map.put(q.get(zi), result_map_i);
            }
            return results_map;
        }
    }


//    public Map<StereoMolecule,Map<LSHProvider.Hit2,List<String>>> searchAllSimilar_Exact(BitSet qbi, double min_tanimoto_sim , int max_results ) {
//        Map<BitSet,List<LSHProvider.Hit2>> hits = new ConcurrentHashMap<>();
//
//        if(!mOutOfMemory) {
//            int hamming_dist = (int) Math.ceil( qbi.cardinality() * (1.0-min_tanimoto_sim) );
//            //List<LSHProvider.Hit2> results = new ArrayList<>();
//            List<LSHProvider.Hit2> results = mANN.exactFindNearestNeighbors2( qbi , hamming_dist , max_results );
//        }
//        else {
//            System.out.println("not supported..");
//        }
//    }


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

    //public static class BitsetLineIterator implements Iterator<BitSet> {

    /**
     * Provides a ResettableBitSetIterator for line-based and
     * tab-separated data files (i.e. tsv or csv files).
     *
     * The input file must be line separated, and the entries on
     * each line must be tab-separated.
     *
     * With the split_pos parameter it is possible to specify a
     * specific entry of the tab-separated line after splitting.
     *
     * The specified split of the line must contain a Base64-encoded
     * bitset.
     *
     */
    public static class BitsetLineIterator implements LSHProvider.ResettableBitSetIterator {

        private File in_file;
        private Reader in = null;
        private LineIterator line_iterator;
        private int split_pos = 0;
        public BitsetLineIterator(File in_file, int split_pos) {
            this.in_file= in_file;
            //this.line_iterator = new LineIterator(in);
            this.split_pos = split_pos;
        }
        @Override
        public void reset() {
            try {
                if(this.in!=null) {
                    this.in.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                this.in = new BufferedReader(new FileReader(in_file));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            this.line_iterator = new LineIterator(this.in);
        }
        @Override
        public boolean hasNext() {
            return this.line_iterator.hasNext();
        }
        @Override
        public BitSet next() {
            String line = this.line_iterator.next();
            String splits[] = line.split("\\t");
            return BitSet.valueOf( Base64.getDecoder().decode( splits[split_pos] ));
        }
    }



    public static void main(String args[]) {

    }

}
