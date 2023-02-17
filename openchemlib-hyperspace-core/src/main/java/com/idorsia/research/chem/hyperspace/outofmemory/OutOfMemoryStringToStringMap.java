package com.idorsia.research.chem.hyperspace.outofmemory;

import java.io.*;
import java.util.*;

public class OutOfMemoryStringToStringMap implements Serializable {

    /**
     * Contains the path to the storage files, INCLUDING the
     * file separator.
     */
    private String m_PathToStorageFile = null;

    /**
     * Number of hash files.
     */
    private int    m_NumHashFiles = 100000;


    public OutOfMemoryStringToStringMap(File f_storage_directory) {
        this.m_PathToStorageFile = f_storage_directory.getAbsolutePath()+File.separator;
    }

    /**
     * Input file must be line-separated, and the
     * key and value must be separated by any \s+
     * NOTE: This implementation supports one-to-
     * many relationships. I.e. if you have the
     * same key pointing to multiple Strings,
     * the data structure will return all value
     * Strings.
     *
     * @param f_input
     */
    public void initFromFile(File f_input, File f_storage_directory, int num_files) throws IOException {
        // 0. check that storage_directory exists and is empty
        if(!f_storage_directory.isDirectory()) {
            throw new IOException("Storage directory does not exist");
        }
        if(f_storage_directory.listFiles().length>0) {
            throw new IOException("Storage directory is not empty");
        }

        this.m_PathToStorageFile = f_storage_directory + File.separator;

        this.m_NumHashFiles = num_files;

        // 1. read through files and process them :)
        BufferedReader in = new BufferedReader( new FileReader( f_input ) );

        Set<String> created_files = new HashSet<>();

        // this is the number of lines that we read, until we start the next set of write operations
        //int config_batchsize = 1000000;
        int config_batchsize = 2000000;

        String line = null;

        Map<String, List<String[]>> data_batch = new HashMap<>();
        int cnt_batch = 0;

        long bytes_to_process = f_input.length();
        long bytes_processed  = 0;

        while( (line=in.readLine()) != null ) {
            bytes_processed += line.getBytes().length;
            line = line.trim();
            if(line.isEmpty()) {continue;}

            // split into two parts:
            String splits[] = line.split("\\s+");
            if(splits.length<2) {
                System.out.println("[WARN] Line splitting did yield less than two splits -> skip line");
                continue;
            }
            if(splits.length>2) {
                System.out.println("[WARN] Line splitting did yield more than two splits -> consider only first part of value");
            }

            String key = splits[0]; String value = splits[1];

            // determine File into which we hash this
            String filename = stringToFilename(key);
            // put into batch
            if(!data_batch.containsKey(filename)){data_batch.put(filename,new ArrayList<>());}
            data_batch.get( filename ).add( new String[]{key,value} );

            // write batch?
            cnt_batch++;
            if(cnt_batch>=config_batchsize) {
                System.out.println("Write batch -> Progress: "+(1.0*bytes_processed / (1.0*bytes_to_process)));
                writeBatch(data_batch);
                cnt_batch = 0;
            }
        }
        // write last batch
        writeBatch(data_batch);
    }

    private static void writeBatch(Map<String,List<String[]>> batch_data) throws IOException {
        for(String ffi : batch_data.keySet()) {
            File fi = new File(ffi);
            fi.createNewFile(); // this also works if the file already exists.
            BufferedWriter out = new BufferedWriter( new FileWriter(fi,true) );
            for(String pi[] : batch_data.get(ffi)) {
                out.write(pi[0]+"\t"+pi[1]+"\n");
            }
            out.flush();
            out.close();
        }
        batch_data.clear();
    }

    public String stringToFilename(String s) {
        int hash = Math.abs(s.hashCode()) % this.m_NumHashFiles;
        return this.m_PathToStorageFile + "smh_"+hash+".xm";
    }

    public Map<String,List<String>> get(List<String> keys) throws IOException {
        Map<String,List<String>> data = new HashMap<>();

        // prepare the fetch, by sorting by files that we need:
        Map<String, List<String>> keys_by_filename = new HashMap<>();
        for(String ki : keys) {
            // determine File into which we hash this
            String filename = stringToFilename(ki);
            // put into batch
            if(!keys_by_filename.containsKey(filename)){keys_by_filename.put(filename,new ArrayList<>());}
            keys_by_filename.get(filename).add(ki);
        }

        // resolve:
        for(String fi_name : keys_by_filename.keySet()) {
            File fi = new File(fi_name);
            Map<String,List<String>> data_i = loadFile(fi);
            for(String ki : keys_by_filename.get(fi_name)) {
                data.put(ki,data_i.get(ki));
            }
        }

        return data;
    }

    public Map<String,List<String>> loadFile(File fi) throws IOException {
        // 1. read file
        Map<String,List<String>> maps = new HashMap<>();
        BufferedReader in = new BufferedReader(new FileReader(fi));
        String line = null;
        while( (line=in.readLine()) != null ) {
            if(line.isEmpty()){continue;}
            String splits[] = line.split("\\t");
            if(!maps.containsKey(splits[0])){maps.put(splits[0],new ArrayList<>());}
            maps.get(splits[0]).add(splits[1]);
        }
        in.close();
        return maps;
    }

    public static void main(String args[]) {
        // 1. create large input file
        File fi_test = new File("C:\\temp\\large_test_input.txt");
        try {
            Random r = new Random();
            BufferedWriter out = new BufferedWriter(new FileWriter(fi_test));
            for(int zi=0;zi<2500000;zi++) {
                out.write( (r.nextInt())+""+(r.nextInt()) +"\t" +(r.nextInt())+""+(r.nextInt())+""+r.nextInt());
                out.write("\n");
            }
            out.flush();
            out.close();


            // 2. init OutOfMemory Map:
            File testdir = new File("C:\\temp\\testdir_"+r.nextInt());
            testdir.mkdir();
            OutOfMemoryStringToStringMap xmap = new OutOfMemoryStringToStringMap(testdir);
            xmap.initFromFile(fi_test,testdir,100000);

            // 3. test it :)
            BufferedReader in = new BufferedReader(new FileReader(fi_test));
            Map<String,List<String>> ground_truth = new HashMap<>();
            for(int zi=0;zi<10000;zi++) {
                String line = in.readLine();
                String splits[] = line.split("\\s+");
                if(!ground_truth.containsKey(splits[0])){ground_truth.put(splits[0],new ArrayList<>());}
                ground_truth.get(splits[0]).add(splits[1]);
            }

            // 4. fetch 10k:
            List<String> keys_to_fetch = new ArrayList<>(ground_truth.keySet());
            long ts_a = System.currentTimeMillis();
            Map<String,List<String>> data_fetched = xmap.get(keys_to_fetch);
            long ts_b = System.currentTimeMillis();
            System.out.println("Fetch took : "+ (ts_b-ts_a) +" ms");

            // check for correctness:
            boolean data_correct = true;
            for(String ki : keys_to_fetch){
                Set<String> sa = new HashSet<>( data_fetched.get(ki) );
                Set<String> sb = new HashSet<>( ground_truth.get(ki) );
                data_correct   &= sa.equals(sb);
            }
            System.out.println("Data: "+(data_correct?"OK":"BAD!!!"));

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
