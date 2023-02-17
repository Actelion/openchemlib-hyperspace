package com.idorsia.research.chem.hyperspace.outofmemory;

import java.io.*;
import java.util.*;

public class MergeSortExternal {


    /**
     * Merges sorted files a and b and writes them (sorted) into file filename_out.
     *
     * @param a
     * @param b
     * @param filename_out
     */
    public static void mergeFiles(File a, File b, String filename_out) {


        try {
            FileReader fin_a = new FileReader(a);
            BufferedReader bin_a = new BufferedReader(fin_a);

            FileReader fin_b = new FileReader(b);
            BufferedReader bin_b = new BufferedReader(fin_b);


            /**
             * Buffer size of the BufferedIterator (i.e. lines read at once when the buffer is empty.)
             */
            int buffer = 4096;

            BufferedIterator bia = new BufferedIterator(a,buffer);
            BufferedIterator bib = new BufferedIterator(b,buffer);


            FileWriter fw = null;
            BufferedWriter bw = null;

            fw = new FileWriter(filename_out);
            bw = new BufferedWriter(fw);

            mergeListsAndWrite(bia,bib,bw);

            bw.flush();
            bw.close();

            bia.close();
            bib.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Implements a buffered iterator with a buffer of buffersize lines.
     */
    static class BufferedIterator implements Iterator {

        File In;

        int Buffersize = 1024;

        FileReader fin = null;
        BufferedReader bin = null;


        LinkedList<String> current_lines = new LinkedList<>();

        // will be set when the last line of the BufferedReader is read (i.e. when it returns null).
        boolean buffered_reader_empty = false;

        public BufferedIterator(File in , int buffersize) {
            In = in;
            Buffersize = buffersize;

            try {
                fin = new FileReader(in);
                bin = new BufferedReader(fin);

                fillBuffer();
            }
            catch(Exception e){
                e.printStackTrace();
            }
        }

        private void fillBuffer() throws Exception {
            while(Buffersize>current_lines.size()){
                String line_next = bin.readLine();
                if(line_next!=null) {
                    current_lines.add(line_next);
                }
                else{
                    buffered_reader_empty = true;
                    break;
                }
            }
        }

        @Override
        public boolean hasNext() {
            return ! ( buffered_reader_empty && current_lines.size()==0);
        }

        @Override
        public Object next() {
            if(current_lines.size()==0){
                try {
                    fillBuffer();
                }
                catch(Exception e){
                    e.printStackTrace();
                }
            }

            String si = current_lines.pollFirst();
            return si;
        }


        public void close() throws IOException {
            this.bin.close();
            this.fin.close();
        }
    }


    /**
     *
     * Performs a mergesort on iterators ia and ib and writes the result in to writer.
     *
     * @param ia
     * @param ib
     * @param writer
     * @param <T>
     * @throws IOException
     */
    public static <T extends Comparable> void mergeListsAndWrite(Iterator<T> ia, Iterator<T> ib , BufferedWriter writer ) throws IOException {

        //List<Comparable> sorted = new ArrayList<>();

        Comparable current_a = ia.next();
        Comparable current_b = ib.next();

        // NEEDED TO ADD THE LAST ELEMENTS..
        boolean last_added_a = false;
        boolean last_added_b = false;

        while( ia.hasNext() || ib.hasNext() ) {
            int comparison = current_a.compareTo(current_b);
            if(comparison==0){
                //sorted.add(current_a);
                //sorted.add(current_b);
                writer.append( current_a +"\n" );
                writer.append( current_b +"\n" );

                current_a = ia.next();
                current_b = ib.next();

                last_added_a = true;
                last_added_b = true;}

            if(comparison<0){
                if(ia.hasNext()) {
                    writer.append( current_a+"\n");
                    current_a = ia.next();
                    last_added_a = true;
                }
                else{
                    // think we cannot end up here..
                    writer.append( current_b+"\n");
                    current_b = ib.next();
                    last_added_b = true;
                }
            }
            else if(comparison>0){
                if(ib.hasNext()) {
                    writer.append( current_b+"\n");
                    current_b = ib.next();
                    last_added_b = true;
                }
                else{
                    writer.append( current_b+"\n");
                    current_a = ia.next();
                    last_added_a = true;
                }
            }
        }


        // add the last element!
        if(last_added_a&&last_added_b){
            writer.append(current_a +"\n"+current_b);
        }
        else if(last_added_a) {
            writer.append("" + current_a);
        }
        else if(last_added_b){
            writer.append(""+current_b);
        }


    }


    /**
     *
     * Splits the file in into sorted files of nlines size.
     *
     * Creates files "sort_data.xxx" and returns the list of filenames-
     *
     * @param in File to split and sort
     * @param nlines number of lines per file
     * @return List of filenames of the created files
     */
    public static List<String> splitAndSort(File in , int nlines ){

        List<String> filenames = new ArrayList<>();

        FileReader fin = null;
        BufferedReader bin = null;
        try {
            fin = new FileReader(in);
            bin = new BufferedReader(fin);


            int split_cnt     = 0;
            int current_lines = 0;

            //StringBuffer sb = new StringBuffer();
            List<String> line_buffer = new ArrayList<>();

            FileWriter     f_out = null;
            BufferedWriter b_out = null;

            boolean first = true;


            String line_i = null;
            while( (line_i = bin.readLine()) != null ){
                if(current_lines==nlines){
                    current_lines = 0;
                    first= false;
                }
                if(current_lines==0){
                    if(!first){

                        System.out.println("Sort batch, n="+line_buffer.size());
                        long timestamp_sort = System.currentTimeMillis();
                        // sort data..
                        //line_buffer.sort(String::compareTo);
                        String buffer_array[] = new String[line_buffer.size()];
                        line_buffer.toArray(buffer_array);
                        Arrays.parallelSort(buffer_array);

                        System.out.println("Took "+(System.currentTimeMillis()-timestamp_sort)+" ms");

                        // write data..
                        long timestamp_write = System.currentTimeMillis();
                        System.out.println("Write file: "+"sort_data."+split_cnt);
                        f_out = new FileWriter("sort_data."+split_cnt);
                        b_out = new BufferedWriter(f_out);
                        filenames.add("sort_data."+split_cnt);

                        //for( String ti : line_buffer){
                        for(int zi=0;zi<buffer_array.length;zi++){
                            String ti = buffer_array[zi];
                            b_out.write( ti );
                            if(ti!=line_buffer.get(line_buffer.size()-1)){
                                b_out.write("\n");
                            }
                        }
                        b_out.flush();
                        b_out.close();

                        System.out.println("done writing.. took "+(System.currentTimeMillis()-timestamp_write)+" ms");

                        // clear line buffer..
                        line_buffer = new ArrayList<>();
                    }
                    split_cnt++;
                }
                line_buffer.add(line_i);
                current_lines++;
            }

            // write the last buffer..
            System.out.println("Last Batch!..");
            System.out.println("Sort batch, n="+line_buffer.size());
            long timestamp_sort = System.currentTimeMillis();
            // sort data..
            //line_buffer.sort(String::compareTo);
            String buffer_array[] = new String[line_buffer.size()];
            line_buffer.toArray(buffer_array);
            Arrays.parallelSort(buffer_array);

            System.out.println("Took "+(System.currentTimeMillis()-timestamp_sort)+" ms");

            // write data..
            long timestamp_write = System.currentTimeMillis();
            System.out.println("Write file: "+"sort_data."+split_cnt);
            f_out = new FileWriter("sort_data."+split_cnt);
            b_out = new BufferedWriter(f_out);
            filenames.add("sort_data."+split_cnt);

            //for( String ti : line_buffer){
            for(int zi=0;zi<buffer_array.length;zi++){
                String ti = buffer_array[zi];
                b_out.write( ti );
                if(ti!=line_buffer.get(line_buffer.size()-1)){
                    b_out.write("\n");
                }
            }
            b_out.flush();
            b_out.close();

            System.out.println("done writing.. took "+(System.currentTimeMillis()-timestamp_write)+" ms");

        }
        catch(Exception e){
            System.out.println(e);
        }

        System.out.println("Done with split and sort!");

        return filenames;
    }



    /**
     * Example program that uses the functions in this class to sort the lines in the file test_01.data, by first
     * splitting and sorting them into chunks (containing n_lines_per_chunk lines).
     *
     * The procedure creates temporary files with filenames file_merged_xxx, the file with the highest _xxx contains
     * the fully sorted list.
     *
     * @param args
     */
    public static void main(String args[]){

        String infile = "test_01.data";

        File f = new File(infile);
        System.out.println("Sorting file of size: "+f.length()/(1024*1024.0) + " MB");


        long timestamp = System.currentTimeMillis();

        int n_lines_per_chunk = 20000000;

        //List<String> filenames = splitAndSort(new File(infile), 10000000);
        List<String> filenames = splitAndSort(new File(infile), n_lines_per_chunk);
        //List<String> filenames = splitAndSort(new File(infile), 100);


        //File a = new File("sort_data.1");
        //File b = new File("sort_data.2");
        //String f_out = "sort_data.stage.1.1";
        //mergeFiles(a,b,f_out);

        int new_file_idx = 1234;


        while( filenames.size()>1 ){

            List<String> new_filenames = new ArrayList<>();

            // merge themz..
            int merges = filenames.size() / 2;

            for(int zi=0;zi<merges;zi++){
                String in_a = filenames.get( zi*2+0 );
                String in_b = filenames.get( zi*2+1 );
                String out_new = "file_merged_"+new_file_idx;
                new_filenames.add(out_new);
                new_file_idx++;

                // merge..
                System.out.println("Merge: "+in_a+" with "+in_b+" into "+out_new);
                File fa = new File(in_a);
                File fb = new File(in_b);
                mergeFiles(fa,fb,out_new);

                System.out.println("done!");
            }

            // ! In case of uneven number we have to add the last one !..
            if( filenames.size()%2>0) {
                new_filenames.add(filenames.get(filenames.size()-1));
            }
            filenames = new_filenames;
        }

        System.out.println("All done!");
        System.out.println("Time needed: "+ ((System.currentTimeMillis()-timestamp)/1000.0) +"s");
        System.out.println("Filename of fully sorted file: "+filenames.get(0));


        if(false)
        {
            try {
                File f_test = new File(filenames.get(0));
                checkIfReallySorted(f_test);
            }
            catch(Exception e){
                e.printStackTrace();
            }
        }


    }



    /*
    public static class MergeTask {
        String f_in_a;
        String f_in_b;
        String f_out;
        public MergeTask( String fia , String fib , String fout ){
            f_in_a = fia;
            f_in_b = fib;
            f_out = fout;
        }
    }
    */

    public static void main_verity(String args[]){
        File f = new File("file_merged_1241");
        try {
            checkIfReallySorted(f);
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public static boolean checkIfReallySorted(File f) throws Exception {

        FileReader     fin = new FileReader(f);
        BufferedReader bin = new BufferedReader(fin);

        String greatest_string = null;

        int line_i = 0;
        String next_line = null;
        while( (next_line=bin.readLine()) !=null ){
            if(greatest_string==null){
                greatest_string=next_line;
            }
            else{
                if(greatest_string.compareTo(next_line)>0){
                    System.out.println("error in line "+line_i+" :\n"+greatest_string+" > "+next_line);
                    System.out.println("should not be..");
                    return false;
                }
            }
            line_i++;

            if(line_i%1000000 == 0){
                System.out.println("tested first "+line_i+" lines, still ok!");
            }
        }

        System.out.println("Tested file is ok! All "+ line_i+" lines sorted alphabetically!");

        bin.close();
        fin.close();

        return true;
    }
}
