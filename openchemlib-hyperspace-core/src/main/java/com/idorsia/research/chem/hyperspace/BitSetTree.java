package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


/**
 * NOTE: This BitSetTree supports three different modes:
 *
 * 1: Leafs in memory (default)
 * 2: Leafs on hard disk
 * 3: Leafs in zip file
 *
 * This is directly implemented in the (leaf) nodes:
 * in both cases, isLeaf() is equivalent to having the ".bit" integer < 0
 * in case 1, the ".leaf_data" ArrayList is non-null (it can be empty, though), and the datafile is null.
 * in case 2, the ".leaf_data" ArrayList is null, and the datafile is non-null.
 *
 *
 * The computation of either tree is done via function createTree
 *
 *
 */

import java.io.*;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

public class BitSetTree implements Serializable {

    private static final long serialVersionUID = 2219385776457978215L;

    public enum STORAGE_MODE { MEMORY , FILE , ZIP };

    Node root;

    BitSetTree(Node root) {
        this.root = root;
    }

    BitSetTree(String serialized) {
        initFromSerialized(serialized);
    }

    /**
     * NOTE!
     * If the recursion finds the first superset in a leaf node,
     * then not all entries of the first_superset node will contain
     * hits! I.e. it is necessary to call collectHits after this!
     *
     *
     * @param b
     * @param first_superset
     * @return

     */
    public boolean testSubset(BitSet b, Node[] first_superset) {
        return this.root.testSubset(b,first_superset);
    }

    public static final class Node implements Serializable {
//        private static final long serialVersionUID = 6612754111018245823L;  // temporarily used for Synple
        private static final long serialVersionUID = -6432137482703457104L;   // added 17-May-2024, TLS, evidently used historically

        // if -1, then this is a leaf
        public int bit;

        BitSet bits_1;
        BitSet bits_0;

        Node left = null;
        Node right = null;
        private List<BitSet> leaf_data = null;

        private String leaf_data_file = null;

        private STORAGE_MODE storage;

        /**
         * This is required if we use storage mode ZIP.
         */
        transient private ZipFile zip_file = null;
        private String  zip_file_name = null;

        public Node(int bit, BitSet bits_0, BitSet bits_1, Node left, Node right, List<BitSet> data) {
            this.bit = bit;
            this.bits_0 = bits_0;
            this.bits_1 = bits_1;
            this.left = left;
            this.right = right;
            this.leaf_data = data;
            this.leaf_data_file = null;
            this.storage = STORAGE_MODE.MEMORY;
        }

        /**
         * Creates the node object for the hard disk stored nodes.
         *
         * @param bit
         * @param bits_0
         * @param bits_1
         * @param left
         * @param right
         * @param leaf_data_file
         */
        public Node(int bit, BitSet bits_0, BitSet bits_1, Node left, Node right, String leaf_data_file, STORAGE_MODE storage,ZipFile zip_file,String zip_file_name) {
            this.bit = bit;
            this.bits_0 = bits_0;
            this.bits_1 = bits_1;
            this.left = left;
            this.right = right;
            this.leaf_data = null;
            this.leaf_data_file = leaf_data_file;
            this.storage = storage;
            this.zip_file = zip_file;
            this.zip_file_name = zip_file_name;
        }





        public static Node createHardDiskStoredLeafNode(File directory, BitSet bits_0, BitSet bits_1, String tree_pos, List<BitSet> data) {
//            String file_id = "b0_"+ bits_0.stream().mapToObj( bi -> ""+bi ).collect(Collectors.joining("_"))+
//                    "b1_"+ bits_0.stream().mapToObj( bi -> ""+bi ).collect(Collectors.joining("_"))+".leaf";

            String file_id = tree_pos + ".leaf";
            String str_data = listOfBitSetsToString(data);

            // write file
            String file_path = directory.getAbsolutePath()+File.separator+file_id;
            File fi = new File(file_path);
            try {
                BufferedWriter out = new BufferedWriter(new FileWriter(fi));
                out.write(str_data);
                out.flush();
                out.close();
            }
            catch(IOException ex) {
                ex.printStackTrace();
            }
            return new Node(-1,bits_0,bits_1,null,null,file_path, STORAGE_MODE.FILE, null, null);
        }

        public static Node createZipFileStoredLeafNode(String temp_file_path, String zip_file_path, ZipOutputStream z_out, BitSet bits_0, BitSet bits_1, String tree_pos, List<BitSet> data) {
//            String file_id = "b0_"+ bits_0.stream().mapToObj( bi -> ""+bi ).collect(Collectors.joining("_"))+
//                    "b1_"+ bits_0.stream().mapToObj( bi -> ""+bi ).collect(Collectors.joining("_"))+".leaf";

            String file_id = tree_pos + ".leaf";
            String str_data = listOfBitSetsToString(data);

            // write file
            //String file_path = directory.getAbsolutePath()+File.separator+file_id;
            String file_path = "data/"+File.separator+file_id;
            //File fi = new File(file_path);

            try {
                z_out.putNextEntry(new ZipEntry(file_path));
                BufferedWriter out = new BufferedWriter(new OutputStreamWriter(z_out));
                out.write(str_data);
                out.flush();
                //out.close();
            }
            catch(IOException ex) {
                ex.printStackTrace();
            }
            return new Node(-1,bits_0,bits_1,null,null,file_path, STORAGE_MODE.ZIP, null, zip_file_path);
        }



        public boolean isLeaf() {
            return this.bit < 0;
        }

        public void setLeft(Node l) {
            this.left = l;
        }

        public void setRight(Node r) {
            this.right = r;
        }

        public int depth() {
            return this.bits_0.cardinality() + this.bits_1.cardinality();
        }

        public Node(String serialized, Node left, Node right) {
            initFromSerialized(serialized, left, right);
        }

        /**
         * returns the leaf data. For the memory case, this data
         * is loaded from memory, for the non-memory case we load
         * it from the file.
         *
         * @return
         */
        public List<BitSet> getLeafData() {
            if(this.leaf_data!=null) {
                return this.leaf_data;
            }
            else {
                if(this.storage == STORAGE_MODE.FILE) {
                    // load:
                    File fi = new File(this.leaf_data_file);
                    String data = null;
                    try {
                        BufferedReader in = new BufferedReader(new FileReader(fi));
                        data = in.readLine();
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                    List<BitSet> n_leaf_data = stringToListOfBitSets(data);
                    return n_leaf_data;
                }
                else if(this.storage == STORAGE_MODE.ZIP) {
                    synchronized(this) {
                        if(zip_file == null) {
                            File fi = new File(this.zip_file_name);
                            try {
                                this.zip_file = new ZipFile(fi);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                    String data = null;
                    try {
                        BufferedReader bin = new BufferedReader( new InputStreamReader( this.zip_file.getInputStream( this.zip_file.getEntry( this.leaf_data_file ) ) ) );
                        data = bin.readLine();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    List<BitSet> n_leaf_data = stringToListOfBitSets(data);
                    return n_leaf_data;
                }
            }
            return null;
        }

        /**
         * NOTE!
         * If the recursion finds the first superset in a leaf node,
         * then not all entries of the first_superset node will contain
         * hits! I.e. it is necessary to call collectHits after this!
         *
         *
         * @param b
         * @param first_superset
         * @return
         */
        public boolean testSubset(BitSet b , Node[] first_superset) {

            first_superset[0] = null;

            if(this.bit<0) {
                List<BitSet> n_leaf_data = getLeafData();

                //if(this.leaf_data.isEmpty()) {
                if(n_leaf_data.isEmpty()) {
                    return false;
                }
                boolean subset = true;
                //for(BitSet bsi : leaf_data) {
                for(BitSet bsi : n_leaf_data) {
                    BitSet ti = (BitSet) bsi.clone();
                    ti.or(b);
                    if(ti.equals(bsi)) {
                        first_superset[0] = this;
                        return true;
                    }
                }
                return false;
            }

            BitSet test = (BitSet) bits_1.clone();
            test.or(b);
            if( test.equals(bits_1) ) {
                first_superset[0] = this;
                return true;
            }


            if( b.get(this.bit) ) {
                return this.right.testSubset( b , first_superset);
            }
            else {
                return this.right.testSubset( b , first_superset) || this.left.testSubset( b , first_superset);
            }
        }

        /**
         * Collect supersets of the query bs.
         *
         * @param bs
         * @return
         */
        public void collectSuperSets(BitSet bs, List<BitSet> supersets ) {
            if(this.bit<0) {
                //LinkedList<BitSet> supersets = new LinkedList<>();
                List<BitSet> n_leaf_data = this.getLeafData();
                //for(BitSet bsi : this.leaf_data) {
                for(BitSet bsi : n_leaf_data) {
                    BitSet ti = (BitSet) bsi.clone();
                    ti.or(bs);
                    if(ti.equals(bsi)) {
                        supersets.add(bsi);
                    }
                }
            }
            else {
                if(bs.get(this.bit)) {
                    this.right.collectSuperSets(bs, supersets);
                }
                else {
                    this.left.collectSuperSets(bs, supersets);
                    this.right.collectSuperSets(bs, supersets);
                }
            }
        }

        public boolean checkAllAreSuperset(BitSet q) {
            if(this.bit<0) {
                BitSet ti = (BitSet) this.bits_1.clone();
                ti.or(q);
                boolean is_superset = ti.equals(this.bits_1);
                if(is_superset) {return true;} // this is only sufficient condition, not necessary.
                // leaf.. check all
                is_superset = true;
                for(BitSet bi : this.leaf_data) {
                    ti = (BitSet) bi.clone();
                    ti.or(q);
                    if(!ti.equals(bi)) { is_superset = false; break; }
                }
                return is_superset;
            }
            if(q.get(this.bit)) {
                if(this.left==null) { return true; }
                BitSet empty = new BitSet(q.size());
                // we have to check for any elements in this subtree, therefore the empty bitset.
                if(this.left.getSuperSetIterator(empty).hasNext()) {
                    return false;
                }
                return true;
            }
            return this.right.checkAllAreSuperset(q) && this.left.checkAllAreSuperset(q);
        }


        public SuperSetIterator getSuperSetIterator(BitSet q) {
            return new SuperSetIterator(this,q);
        }

        public int countAll() {
            if(this.isLeaf()) {
                List<BitSet> n_leaf_data = this.getLeafData();
                //return this.leaf_data.size();
                return n_leaf_data.size();
            }
            int ca = (this.left!=null)?this.left.countAll():0;
            int cb = (this.right!=null)?this.right.countAll():0;
            return ca+cb;
        }

        /**
         * serializes all fields up to the children
         *
         * @return
         */
        public String serializeWithoutChildren() {

            StringBuilder sb = new StringBuilder();
            sb.append(this.bit);
            sb.append("   ");

            if(isLeaf()) {
                //sb.append("L");
                //sb.append("   ");
                sb.append( bitSetToString(this.bits_0) );
                sb.append("   ");
                sb.append( bitSetToString(this.bits_1) );
                sb.append("   ");
                sb.append( listOfBitSetsToString(this.leaf_data) );
            }
            else {
                //sb.append("N");
                sb.append( bitSetToString(this.bits_0) );
                sb.append("   ");
                sb.append( bitSetToString(this.bits_1) );
            }
            return sb.toString();
        }

        private void initFromSerialized(String serialized, Node left, Node right) {
            this.left  = left;
            this.right = right;

            String splits[] = serialized.split("   ");

            this.bit = Integer.parseInt(splits[0]);

            if(bit < 0) {
                // init a leaf
                this.bits_0    = stringToBitSet(splits[1]);
                this.bits_1    = stringToBitSet(splits[2]);
                this.leaf_data = stringToListOfBitSets(splits[3]);
            }
            else {
                this.bits_0    = stringToBitSet(splits[1]);
                this.bits_1    = stringToBitSet(splits[2]);
                this.leaf_data = null;
            }

            throw new Error("Deprecated..");
        }

    }

    public static BitSetTree createTree( Set<BitSet> bitsets , int num_bits , int binsize) {
        return createTree(bitsets,num_bits,binsize,null,null,null);
    }

    /**
     *
     * @param bitsets
     * @param num_bits
     * @param binsize
     * @param file_directory if non-null, then we try to create data files for the leafs in the specified directory
     * @return
     */
    public static BitSetTree createTree( Set<BitSet> bitsets , int num_bits , int binsize, String file_directory, String zip_file_name, ZipOutputStream zip_out) {
        Node root = split_recursively(bitsets, new BitSet(num_bits) , new BitSet(num_bits) , num_bits,binsize, "r" , file_directory, zip_file_name, zip_out);
        return new BitSetTree(root);
    }


    /**
     *
     * Input file must be line-separated, each line contains a Base64 encoded bitset.
     *
     * @param input
     * @param num_bits
     * @param binsize
     * @param file_directory if non-null, then we try to create data files for the leafs in the specified directory
     * @return
     */
    public static BitSetTree createTreeOutOfMemory( File input , int num_bits , int binsize, String file_directory) throws IOException {
        Node root = split_recursively_out_of_memory(input, new BitSet(num_bits) , new BitSet(num_bits) , num_bits,binsize,"r",file_directory,null, null);
        return new BitSetTree(root);
    }


    private static long num_bytes_switch_to_in_memory_computation = 256000000 ;

    /**
     *
     * Input file must be line-separated, each line contains a Base64 encoded bitset.
     *
     * @param input
     * @param num_bits
     * @param binsize
     * @param
     * @return
     */
    public static BitSetTree createTreeOutOfMemory_ZipFile( File input , int num_bits , int binsize, String temp_directory, String zip_output_file) throws IOException {
        ZipOutputStream z_out = new ZipOutputStream(new FileOutputStream(zip_output_file));
        Node root = split_recursively_out_of_memory(input, new BitSet(num_bits) , new BitSet(num_bits) , num_bits,binsize,"r",temp_directory,zip_output_file,z_out);
        z_out.flush();
        z_out.close();
        return new BitSetTree(root);
    }


    private static Random random = new Random();


    public static Node split_recursively(Collection<BitSet> bi,BitSet bits_0, BitSet bits_1, int num_bits, int binsize, String tree_pos,
                                         String file_directory, String zip_file_name, ZipOutputStream zip_out ) {
        if(bi.size() <= binsize) {
            STORAGE_MODE storage = STORAGE_MODE.MEMORY;
            if(file_directory!=null) {storage = STORAGE_MODE.FILE;}
            if(zip_out!=null){storage = STORAGE_MODE.ZIP;}
            //STORAGE_MODE storage = (zip_out==null) ? STORAGE_MODE.FILE : STORAGE_MODE.ZIP;

            if(storage == STORAGE_MODE.MEMORY) {
                return new Node(-1, bits_0, bits_1, null, null, new ArrayList<>(bi));
            }
            else if(storage == STORAGE_MODE.ZIP) {
                return Node.createZipFileStoredLeafNode(file_directory, zip_file_name , zip_out, bits_0, bits_1, tree_pos, new ArrayList<>(bi));
            }
            else if(storage == STORAGE_MODE.FILE){
                // store data file for the leaf on disk:
                File fi = new File(file_directory);
                return Node.createHardDiskStoredLeafNode(fi,bits_0,bits_1,tree_pos,new ArrayList<>(bi));
            }
        }

        List<Integer> possible_splits = new ArrayList<>();
        for(int zi=0;zi<num_bits;zi++) { if( (!bits_0.get(zi)) && (!bits_1.get(zi)) ) { possible_splits.add(zi); } }
        Collections.shuffle( possible_splits , random );

        double best_split     = 0.0;
        int    best_split_bit = -1;
        for(int zi=0;zi<possible_splits.size();zi++)
        {
            //int split = random.nextInt(num_bits);
            int split = possible_splits.get(zi);

            int sa = 0;
            for(BitSet bsi : bi) {
                sa += (bsi.get(split))?1:0;
            }
            double split_value = 1.0*sa / bi.size();
            double split_score = Math.min( split_value , 1.0-split_value );
            if( split_score > best_split ) {
                best_split_bit = split;
                best_split     = split_score;
            }
            if( best_split > 0.42 ) {
                break;
            }
        }
        //System.out.println("Node size= " + bi.size() + " Split at bit "+best_split_bit+" , score = " + (best_split) );

        if(best_split_bit<0) {
            System.out.println("wut?");
        }
        else {
            if(true) {
                System.out.println("SplitScore: "+best_split+"  (Size="+bi.size()+",level="+(bits_0.cardinality()+bits_1.cardinality())+")");
            }
        }

        List<BitSet> bs_a = new ArrayList<>();
        List<BitSet> bs_b = new ArrayList<>();
        for(BitSet bsi : bi) {
            if(bsi.get(best_split_bit)) {
                bs_b.add(bsi);
            }
            else {
                bs_a.add(bsi);
            }
        }

        BitSet bits_0_left  = (BitSet) bits_0.clone();
        BitSet bits_1_left  = (BitSet) bits_1.clone();
        BitSet bits_0_right = (BitSet) bits_0.clone();
        BitSet bits_1_right = (BitSet) bits_1.clone();

        bits_0_left.set(best_split_bit,true);
        bits_1_right.set(best_split_bit,true);

        Node left  = split_recursively(bs_a,bits_0_left,bits_1_left,num_bits,binsize,tree_pos+"_0",file_directory,zip_file_name,zip_out);
        Node right = split_recursively(bs_b,bits_0_right,bits_1_right,num_bits,binsize,tree_pos+"_1",file_directory,zip_file_name,zip_out);

        return new Node(best_split_bit,bits_0,bits_1,left,right,(List<BitSet>)null);
    }


    /**
     * If zip_out is null, then STORAGE_MODE.FILE will be used, if it is not null,
     * then STOREAGE_MODE.ZIP will be used.
     *
     * @param inputfile
     * @param bits_0
     * @param bits_1
     * @param num_bits
     * @param binsize
     * @param tree_pos
     * @param file_directory for STORAGE_MODE.FILE this is both the temp directory and the directory
     *                       where the final data is stored. For .ZIP this is the path to the temp directory.
     * @param zip_out
     * @return
     * @throws IOException
     */
    public static Node split_recursively_out_of_memory(File inputfile,BitSet bits_0, BitSet bits_1,
                                                       int num_bits, int binsize, String tree_pos,
                                                       String file_directory,
                                                       String zip_file_name, ZipOutputStream zip_out ) throws IOException {

        STORAGE_MODE storage = (zip_out==null) ? STORAGE_MODE.FILE : STORAGE_MODE.ZIP;

        double bits_to_try_ratio = 0.02; // we only consider x of all bits as potential split bits and take the best one
        List<Integer> possible_bits = new ArrayList<>();
        for(int zi=0;zi<num_bits;zi++) {
            if(bits_0.get(zi)||bits_1.get(zi)) {
            }
            else {
                possible_bits.add(zi);
            }
        }
        Collections.shuffle( possible_bits , random );
        possible_bits = possible_bits.subList(0, (int) (bits_to_try_ratio * num_bits ) );

        Base64.Decoder base64_decoder = Base64.getDecoder();

        // read file, and compute bit-counts (and count lines, to see if we are done..)
        long counts_zero[] = new long[num_bits];
        long counts_one[]  = new long[num_bits];

        int line_count = 0;
        try( FileReader in_f = new FileReader(inputfile)) {
            BufferedReader in = new BufferedReader(in_f);

            String line = null;
            int num_characters = 0;
            while ((line = in.readLine()) != null) {
                line_count++;
                if (false) {
                    byte buf[] = base64_decoder.decode(line);
                    BitSet bsi = BitSet.valueOf(buf);
                    for (int zi = 0; zi < num_bits; zi++) {
                        if (bsi.get(zi)) {
                            counts_one[zi]++;
                        } else {
                            counts_zero[zi]++;
                        }
                    }
                } else {
                    byte buf[] = base64_decoder.decode(line);
                    //BitSet bsi = BitSet.valueOf(buf);
                    for (int zzi = 0; zzi < possible_bits.size(); zzi++) {
                        int zi = possible_bits.get(zzi);
                        //boolean bsi_v = bsi.get(zi);
                        boolean bsi_v = isSet(buf, zi);
                        counts_one[zi] += bsi_v ? 1 : 0;
                        counts_zero[zi] += bsi_v ? 0 : 1;
//                    if (bsi.get(zi)) {
//                        counts_one[zi]++;
//                    } else {
//                        counts_zero[zi]++;
//                    }
                    }
                    //num_characters+= line.length() + buf[0];
                }
            }
            //System.out.println("nc="+num_characters);
        }

        if(line_count<binsize) {
            // in this case we read all bitsets into memory and create the leaf:
            ArrayList<BitSet> bitsets = new ArrayList<>();
            try(FileReader in_f = new FileReader(inputfile)) {
                BufferedReader in = new BufferedReader(in_f);
                String line = null;
                while ((line = in.readLine()) != null) {
                    line_count++;
                    byte buf[] = base64_decoder.decode(line);
                    bitsets.add(BitSet.valueOf(buf));
                }
            }

            // delete the temp file (we dont need it anymore..)
            inputfile.delete();

            if(storage == STORAGE_MODE.FILE) {
                // store data file for the leaf on disk:
                File fi = new File(file_directory);
                return Node.createHardDiskStoredLeafNode(fi, bits_0, bits_1, tree_pos, bitsets);
            }
            else if(storage == STORAGE_MODE.ZIP) {
                return Node.createZipFileStoredLeafNode(file_directory, zip_file_name , zip_out, bits_0, bits_1, tree_pos, bitsets);
            }
            else {
                System.out.println("unknown storage mode+ "+storage);
                return null;
            }
        }

        // ok, now we can pick the best bit..
        double best_score = -1;
        int    best_idx   = -1;
        for(int zi=0;zi<num_bits;zi++) {
            double split_value = 1.0*counts_one[zi] / line_count;
            double split_score = Math.min( split_value , 1.0-split_value );
            if(split_score > best_score) {
                best_score = split_score;
                best_idx   = zi;
            }
        }

        // check if we found some split..
        if(best_score<=0) {
            throw new Error("Could not find splitting bit.");
        }
        else {
            System.out.println("SplitScore: "+best_score+"  (Size="+line_count+",level="+(bits_0.cardinality()+bits_1.cardinality())+")");
        }



        String filename_zero = tree_pos+"_0_x"+".tmp";
        String filename_one  = tree_pos+"_1_y"+".tmp";

        File fi_dir = new File(file_directory);
        String file_path_zero = fi_dir.getAbsolutePath()+File.separator+filename_zero;
        String file_path_one  = fi_dir.getAbsolutePath()+File.separator+filename_one;

        File f_a = new File(file_path_zero);
        File f_b = new File(file_path_one);

        //in.close();
        long num_bytes_one  = 0;
        long num_bytes_zero = 0;

        try(FileReader in_f = new FileReader(inputfile)) {
            BufferedReader in = new BufferedReader(in_f);

            try(FileWriter f_out_a = new FileWriter(f_a)) {
                try(FileWriter f_out_b = new FileWriter(f_b)) {
                    BufferedWriter out_zero = new BufferedWriter(f_out_a);
                    BufferedWriter out_one = new BufferedWriter(f_out_b);

                    String line = null;
                    while ((line = in.readLine()) != null) {
                        byte buf[] = base64_decoder.decode(line);
                        BitSet bsi = BitSet.valueOf(buf);
                        if (bsi.get(best_idx)) {
                            out_one.write(line + "\n");
                            num_bytes_one += line.getBytes().length;
                        } else {
                            out_zero.write(line + "\n");
                            num_bytes_zero += line.getBytes().length;
                        }
                    }

                    out_zero.flush();
                    out_one.flush();
                }
            }
        }

        // delete input file:
        boolean deleted_success = inputfile.delete();

        BitSet bits_0_left  = (BitSet) bits_0.clone();
        BitSet bits_1_left  = (BitSet) bits_1.clone();
        BitSet bits_0_right = (BitSet) bits_0.clone();
        BitSet bits_1_right = (BitSet) bits_1.clone();

        bits_0_left.set(best_idx,true);
        bits_1_right.set(best_idx,true);

        Node left  = null;
        Node right = null;
        if(num_bytes_zero>num_bytes_switch_to_in_memory_computation) {
            left = split_recursively_out_of_memory(f_a, bits_0_left, bits_1_left, num_bits, binsize, tree_pos + "_0", file_directory, zip_file_name, zip_out);
        }
        else {
            // in this case we read all bitsets into memory and switch to in memory splitting:
            ArrayList<BitSet> bitsets = new ArrayList<>();
            try(FileReader in_f = new FileReader(f_a)) {
                BufferedReader in = new BufferedReader(in_f);
                String line = null;
                while ((line = in.readLine()) != null) {
                    //line_count++;
                    byte buf[] = base64_decoder.decode(line);
                    bitsets.add(BitSet.valueOf(buf));
                }
            }
            left = split_recursively(bitsets, bits_0_left, bits_1_left, num_bits, binsize, tree_pos + "_0", file_directory, zip_file_name, zip_out) ;
        }
        if(num_bytes_one>num_bytes_switch_to_in_memory_computation) {
            right = split_recursively_out_of_memory(f_b, bits_0_right, bits_1_right, num_bits, binsize, tree_pos + "_1", file_directory, zip_file_name, zip_out);
        }
        else {
            // in this case we read all bitsets into memory and switch to in memory splitting:
            ArrayList<BitSet> bitsets = new ArrayList<>();
            try(FileReader in_f = new FileReader(f_b)) {
                BufferedReader in = new BufferedReader(in_f);
                String line = null;
                while ((line = in.readLine()) != null) {
                    //line_count++;
                    byte buf[] = base64_decoder.decode(line);
                    bitsets.add(BitSet.valueOf(buf));
                }
            }
            right = split_recursively(bitsets, bits_0_right, bits_1_right, num_bits, binsize, tree_pos + "_1", file_directory, zip_file_name, zip_out);;
        }

        return new Node(best_idx,bits_0,bits_1,left,right,(List<BitSet>)null);
    }

    public static boolean isSet(byte[] arr, int bit) {
        int index = bit / 8;  // Get the index of the array for the byte with this bit
        int bitPosition = bit % 8;  // Position of this bit in a byte
        if(index >= arr.length) {
            return false;
        }
        return (arr[index] >> bitPosition & 1) == 1;
    }

    /**
     * Serializes the tree.
     * Note that it is just necessary to serialize the list of nodes.
     * How the nodes are connected can be reconstructed from the
     * BitSets indicating the zero / one splits.
     *
     * @return
     */
    public String serializeToString() {
        List<Node> all_nodes = new ArrayList<>();
        collectNodes_dfs(this.root,all_nodes);

        List<String> node_strings = new ArrayList<>();
        for(Node ni : all_nodes) {
            node_strings.add( ni.serializeWithoutChildren() );
        }

        String data = String.join("    ",node_strings);
        return data;
    }

    private void initFromSerialized(String serialized) {
        // 1. split into nodes and parse all nodes..
        List<Node> nodes = new ArrayList<>();
        String splits[] = serialized.split("    ");

        for(String si : splits) {
            Node ni = new Node(si,null,null);
            nodes.add(ni);
        }

        //System.out.println("ok! assemble tree..");

        Map<BitSet,Map<BitSet,Node>> nodes_sorted = new HashMap<>();
        for(Node ni : nodes) {
            if( !nodes_sorted.containsKey(ni.bits_0) ) { nodes_sorted.put(ni.bits_0,new HashMap<BitSet,Node>()); }
            nodes_sorted.get(ni.bits_0).put(ni.bits_1,ni);
        }

        // now assign children to every node:
        Node root = null;
        for(Node ni : nodes) {
            BitSet b0 = ni.bits_0;
            BitSet b1 = ni.bits_1;

            if(b0.cardinality()==0 && b1.cardinality()==0) {
                root = ni;
            }

            int bit   = ni.bit;
            if(bit<0) {
                continue; // no children..
            }

            BitSet b0_a = (BitSet) b0.clone();
            b0_a.set(bit);
            BitSet b1_a = (BitSet) b1.clone();
            b1_a.set(bit);

            ni.setLeft( nodes_sorted.get(b0_a).get(b1) );
            ni.setRight( nodes_sorted.get(b0).get(b1_a) );
        }
        this.root = root;
    }

    private void collectNodes_dfs(Node n, List<Node> nodes) {
        nodes.add(n);
        if(n.left!=null) {
            collectNodes_dfs(n.left, nodes);
        }
        if(n.right!=null) {
            collectNodes_dfs(n.right, nodes);
        }
    }


    public static void main__(String args[]) {

        Random r = new Random();

        int bits = 64;

        // test:
        List<BitSet> rand_bs = createRandomBitSet(r,5000,bits, 0.7);

        List<BitSet> positive = new ArrayList<>();
        List<BitSet> negative = new ArrayList<>();

        for(int zi=0;zi<128;zi++) {
            BitSet p1 = (BitSet) rand_bs.get(r.nextInt(5000)).clone();
            for(int zj=0;zj<r.nextInt(6);zj++) {
                p1.set(r.nextInt(bits),false);
            }
            positive.add(p1);

            BitSet n1 = (BitSet) rand_bs.get(r.nextInt(5000)).clone();
            n1.set( n1.nextClearBit(r.nextInt(12)) , true );
            negative.add(n1);
        }

        BitSetTree bst = createTree(new HashSet<>(rand_bs),bits,32);

        for(BitSet bsi : positive) {
            Node[] result = new Node[1];
            System.out.println("pos: " + bst.testSubset(bsi,result) );
        }
        for(BitSet bsi : negative) {
            Node[] result = new Node[1];
            System.out.println("neg: " + bst.testSubset(bsi,result) );
        }



    }


    public static class SuperSetIterator implements Iterator<BitSet> {

        Node       n                = null;
        BitSet     q                = null;

        Iterator<BitSet>   current_iterator = null;
        BitSet             current_next  = null;

        // contains the not yet exhausted nodes (fill in reverse order that you want them to be processed!)
        List<Node>       remaining_childs = null;

        public SuperSetIterator(Node n, BitSet q) {
            this.n    = n;
            this.q    = q;

            // fill remaining_childs (we go "1" first, then "0", therefore add in reverse!):
            if(!n.isLeaf()) {
                this.remaining_childs = new ArrayList<>();
                if(n.left  != null) {this.remaining_childs.add(n.left);}
                if(n.right != null) {this.remaining_childs.add(n.right);}
            }
        }

        @Override
        public boolean hasNext() {
            if(this.current_next==null) {
                tryToFindNext();
            }
            return this.current_next!=null;
        }

        private void tryToFindNext() {
            if(this.current_next!=null) {
                // we already have the next..
                return;
            }

            if(this.n.isLeaf()) {
                if(this.current_iterator == null) {
                    this.current_iterator = this.n.leaf_data.iterator();
                }
                // compute superset tests and check for next:
                while( this.current_iterator.hasNext() ) {
                    //this.current_next = this.current_iterator.next();
                    BitSet candidate = this.current_iterator.next();
                    // superset test:
                    BitSet ti = (BitSet) candidate.clone();
                    ti.or(this.q);
                    if(ti.equals(candidate)) {
                        // ok, we have a next and return
                        this.current_next = ti;
                        return;
                    }
                }
                // this leaf is exhausted if we end up here
                this.current_next = null;
                return;
            }
            else {
                // 1. check if we have iterator. if we have, check if we it has next, else set to null
                // 2. if current_iterator == null, check if we have next
                if(this.current_iterator!=null) {
                    if(this.current_iterator.hasNext()) {
                        this.current_next = this.current_iterator.next();
                        return;
                    }
                    else {
                        this.current_iterator = null;
                    }
                }

                while(this.current_next==null) { // && this.remaining_childs.size()>0) {
                    if(this.current_iterator!=null) {
                        if(this.current_iterator.hasNext()) {
                            this.current_next = this.current_iterator.next();
                            return;
                        }
                        else {
                            this.current_iterator = null;
                        }
                    }
                    // we need new iterator?
                    if (this.current_iterator == null) {
                        if (this.remaining_childs.size() > 0) {
                            this.current_iterator = this.remaining_childs.remove(this.remaining_childs.size() - 1).getSuperSetIterator(this.q);
                        } else {
                            this.current_next = null;
                            return;
                        }
                    }
                    // if we end up here, we have next iterator to search through:
                    if(this.current_iterator.hasNext()) {
                        this.current_next = this.current_iterator.next();
                        return;
                    }
                    else {
                        // we continue with this loop..
                        this.current_next = null;
                    }
                }
            }
        }

        @Override
        public BitSet next() {
            BitSet next = this.current_next;
            // Set next to null!!
            this.current_next = null;
            return next;
        }
    }




    public static void main(String args[]) {

        //test_bitsets();
        //test_serialization();
        //test_fetch_supersets();
        //test_bitsets2();
        test_bitsets2_zip();


    }

    public static void test_bitsets2_zip() {
        Random r = new Random();
        List<BitSet> test = createRandomBitSet(r, 2000000,64,0.6);
        // write into file:
        String test_input_file_name = "test_bitsettree_zip_"+r.nextInt(1000)+".data";
        System.out.println("write test file "+test_input_file_name);

        BufferedWriter test_out;
        try {
            test_out = new BufferedWriter(new FileWriter(new File(test_input_file_name)));
            for(BitSet bsi : test) {
                test_out.write( Base64.getEncoder().encodeToString(bsi.toByteArray()) +"\n");
            }
            test_out.flush();
            test_out.close();
            System.out.println("done, now create tree");
        } catch (IOException e) {
            e.printStackTrace();
        }

        File input = new File(test_input_file_name);

        String test_temp_dir = "C:\\temp\\test_zip_bitsettree_temp";
        String test_zip_file = "C:\\temp\\test_zip_bitsettree.zip";

        try {
            File f_temp_dir = new File(test_temp_dir);
            if(!f_temp_dir.isDirectory()) {
                f_temp_dir.mkdir();
            }

            BitSetTree t  = BitSetTree.createTreeOutOfMemory_ZipFile( input, 64,256 , test_temp_dir, test_zip_file );
            System.out.println("done, now test");

            // now test:
            System.out.println("Test:");
            BitSet to_find = new BitSet(64);
            to_find.set(13);
            to_find.set(17);
            to_find.set(19);
            to_find.set(27);

            List<BitSet> supersets = new ArrayList<>();
            for(BitSet ti : test) {
                BitSet tit = (BitSet) ti.clone();
                tit.or(to_find);
                if(tit.equals(ti)) {
                    supersets.add(ti);
                }
            }

            Node supersets_2[] = new Node[1];
            boolean found = t.testSubset(to_find,supersets_2);

            List<BitSet> supersets_3 = new ArrayList<>();
            t.root.collectSuperSets(to_find, supersets_3);

            if(supersets_3.size()==supersets.size()) {
                System.out.println("ok");
            }
            else {
                System.out.println("not ok");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public static void test_bitsets2() {
        Random r = new Random();
        List<BitSet> test = createRandomBitSet(r, 20000,64,0.6);

        BitSetTree t  = BitSetTree.createTree( new HashSet<>(test) , 64,8 );

        BitSet to_find = new BitSet(64);
        to_find.set(13);
        to_find.set(17);
        to_find.set(19);
        to_find.set(27);

        List<BitSet> supersets = new ArrayList<>();
        for(BitSet ti : test) {
            BitSet tit = (BitSet) ti.clone();
            tit.or(to_find);
            if(tit.equals(ti)) {
                supersets.add(ti);
            }
        }

        Node supersets_2[] = new Node[1];
        boolean found = t.testSubset(to_find,supersets_2);

        List<BitSet> supersets_3 = new ArrayList<>();
        t.root.collectSuperSets(to_find, supersets_3);

        if(supersets_3.size()==supersets.size()) {
            System.out.println("ok");
        }
        else {
            System.out.println("not ok");
        }

    }

    public static void test_fetch_supersets() {
        Random r = new Random();
        List<BitSet> test = createRandomBitSet(r,64,16,0.6);

        BitSetTree t = BitSetTree.createTree( new HashSet<>(test) , 16,2 );

        BitSet to_find = test.get(3);
        to_find.set( to_find.nextSetBit(0) , 0 );
        Node supersets[] = new Node[1];
        boolean found = t.testSubset(to_find,supersets);

        //List<BitSet> bs = supersets[0].collectSuperSets(to_find);
        System.out.println("ok");
    }

    public static void test_bitsets() {
        Random r = new Random();
        List<BitSet> test = createRandomBitSet(r,5,13,0.4);
        String test_string = listOfBitSetsToString(test);
        List<BitSet> test2 = stringToListOfBitSets(test_string);
        System.out.println("ok");
    }

    public static void test_serialization() {
        Random r = new Random();
        List<BitSet> test = createRandomBitSet(r,7,64,0.4);

        BitSetTree t = BitSetTree.createTree(new HashSet<>(test),64,2);
        String si = t.serializeToString();

        System.out.println("ok");

        BitSetTree t1_b = new BitSetTree(si);

        System.out.println("ok");

        List<BitSet> test2 = createRandomBitSet(r,1,64,0.4);
        BitSetTree t2 = BitSetTree.createTree(new HashSet<>(test2),64,16);
        String t2_s = t2.serializeToString();
        BitSetTree t2_b = new BitSetTree(t2_s);

        System.out.println("ok");
    }

    public static void test_benchmark_a() {
        Random r = new Random();

        int bits = 512;
        List<BitSet> rand_bs = createRandomBitSet(r,200000,bits, 0.6);

        System.out.println("Create tree:");
        BitSetTree bst = createTree(new HashSet<>(rand_bs),bits,64);
        System.out.println("Create tree: done!");

        List<BitSet> rand_test = createRandomBitSet(r,1000,bits,0.25);

        System.out.println("Start eval:");
        long ts_a = System.currentTimeMillis();

        for(BitSet bsi : rand_test) {
            Node[] result = new Node[1];
            System.out.println(bst.testSubset(bsi,result));
        }

        long ts_b = System.currentTimeMillis();

        System.out.println("Time= "+(ts_b-ts_a));
    }

    public static List<BitSet> createRandomBitSet(Random r, int n, int bits, double density) {
        List<BitSet> rand_bs = new ArrayList<>();
        for(int zi=0;zi<n;zi++) {
            BitSet bsi = new BitSet(bits);
            for(int zj=0;zj<bits;zj++) {
                bsi.set( zj , r.nextFloat() < density );
            }
            rand_bs.add(bsi);
        }
        return rand_bs;
    }

    public static String bitSetToString(BitSet b) {

        StringBuilder sb = new StringBuilder();
        sb.append( b.size() );
        sb.append(" ");
        byte[] ba = b.toByteArray();
        if(ba.length==0){ba = new byte[1];}
        sb.append( Base64.getEncoder().encodeToString(ba));
        return sb.toString();
    }

    public static BitSet stringToBitSet(String s) {
        String splits[] = s.split(" ");
        BitSet bs = BitSet.valueOf( Base64.getDecoder().decode(splits[1]) );
        int size = Integer.parseInt(splits[0]);
        if(bs.size()!=size) {
            BitSet bs2 = new BitSet(size);
            for(int zi=0;zi<size;zi++) {
                bs2.set(zi,bs.get(zi));
            }
            return bs2;
        }
        return bs;
    }

    public static String listOfBitSetsToString(List<BitSet> bitsets ) {
        if(bitsets.isEmpty()) {
            return "_<<EMPTY>>_<<EMPTY>>_";
        }

        List<String> strings = new ArrayList<>();
        for(BitSet bs : bitsets) {
            strings.add(bitSetToString(bs));
        }

        return String.join("  ",strings);
    }

    public static List<BitSet> stringToListOfBitSets(String sb) {
        if(sb.equals("_<<EMPTY>>_<<EMPTY>>_")) {
            return new ArrayList<>();
        }

        String splits[]      = sb.split("  ");
        List<BitSet> bitsets = new ArrayList<>();
        for(String split_i : splits) {
            bitsets.add( stringToBitSet(split_i) );
        }
        return bitsets;
    }
}
