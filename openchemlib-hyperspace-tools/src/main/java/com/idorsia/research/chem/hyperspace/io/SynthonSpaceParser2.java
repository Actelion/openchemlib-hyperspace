package com.idorsia.research.chem.hyperspace.io;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongPFP512;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;
import org.apache.commons.lang3.tuple.Pair;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Semaphore;
import java.util.stream.Stream;
import java.util.zip.GZIPOutputStream;

public class SynthonSpaceParser2 {


    private static volatile int linesProcessed = 0;

    /**
     * Mode can be either "mode_pfp" , "mode_skelspheres", "mode_ffp" or "mode_ppcore"
     *
     * @param file
     * @param mode
     * @return
     * @throws IOException
     */
    public static SynthonSpace parseBuildingBlocks(File file, String mode, int threads, int max_number_of_synthon_sets) throws IOException {

        if (!(mode.equals("PathFp") || mode.equals("FragFp") || mode.equals("mode_pfp") || mode.equals("mode_ppcore") || mode.equals("mode_ffp") || mode.equals("mode_skelspheres"))) {
            System.out.println("Mode not supported!!");
            return null;
        }

        SynthonSpace space = new SynthonSpace();

        if(mode.equals("mode_pfp")) {
            space.setFP(new DescriptorHandlerLongFFP1024_plus("pfp"),1024);
        }
        if(mode.equals("mode_ffp")) {
            space.setFP(new DescriptorHandlerLongFFP1024_plus("ffp"),1024);
        }
        if(mode.equals("mode_ppcore")) {
            space.setFP(new DescriptorHandlerPPCore(),2048);
        }
        if(mode.equals("FragFp")) {
            space.setFP(new DescriptorHandlerLongFFP512() , 512);
        }
        if(mode.equals("PathFp")) {
            space.setFP(new DescriptorHandlerLongPFP512() , 512);
        }

        // 1. just parse the input file:
        long lineCount;
        try (Stream<String> stream = Files.lines(file.toPath(), StandardCharsets.UTF_8)) {
            lineCount = stream.count();
        }

        linesProcessed = 0;

        BufferedReader in = new BufferedReader(new FileReader(file));

        Thread t_progress_monitor = new Thread() {
            @Override
            public void run() {
                while(true) {
                    try {
                        Thread.sleep(5000);
                    } catch (InterruptedException e) {
                        System.out.println("monitor done!"); break;
                    }
                    double ratio = (1.0*linesProcessed) / lineCount;
                    System.out.println("read: "+ String.format( "%4.2f" , (ratio*100) ) + " i.e. "+ (linesProcessed)+" / "+lineCount ) ;
                }
            }
        };
        t_progress_monitor.start();


        // map: rxx_id -> ( synthon_set_id -> (idcode,id)  )
        Map<String, Map<Integer, List<Pair<String,String>>>> enamine_data = new HashMap<>();

        // skip first line:
        in.readLine();

        // read data:
        String line = null;
        while( (line=in.readLine()) != null ) {
            linesProcessed++;
            try {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }

                // split line:
                String splits[] = line.split("\\t");

                SmilesParser spi = new SmilesParser();
                StereoMolecule mi = new StereoMolecule();
                spi.parse(mi,splits[0]);
                String idcode = mi.getIDCode();
                Integer synthon_idx = Integer.parseInt(splits[2]);
                String frag_id = splits[1];
                String rxn_id  = splits[3];

                if(! enamine_data.containsKey(rxn_id)) { enamine_data.put(rxn_id,new HashMap<>()); }
                if(! enamine_data.get(rxn_id).containsKey(synthon_idx) ) { enamine_data.get(rxn_id).put( synthon_idx , new ArrayList<>()); }

                enamine_data.get(rxn_id).get(synthon_idx).add( Pair.of( idcode,frag_id ) );
            }
            catch(Exception ex) {
                ex.printStackTrace();
            }
        }

        t_progress_monitor.interrupt();
        System.out.println("Done reading data, now create synthon space!");

        long numStructures = enamine_data.values().stream().mapToLong( ki -> ki.values().stream().mapToLong( hi -> hi.size() ).reduce( (xa,xb) -> xa*xb ).getAsLong()  ).sum();


        List<String> rxn_list = new ArrayList<>(enamine_data.keySet());

        List<Thread> all_space_creation_threads = new ArrayList<>();
        Semaphore semaphore = new Semaphore(threads);

        for( int zi=0;zi<rxn_list.size();zi++ ) {
            int fzi = zi;
            String ri = rxn_list.get(zi);
            if(enamine_data.get(ri).size()>max_number_of_synthon_sets) {
                // continue;
                System.out.println("skip synthon reaction with "+enamine_data.get(ri).size()+" synthon sets");
                continue;
            }

            Thread ti = new Thread() {
                @Override
                public void run() {

                    Map<Integer,List<Object>> synthon_sets = new HashMap<>();
                    Map<String,String> idc_to_identifier = new HashMap<>();

                    for(int ki : enamine_data.get(ri).keySet() ) {
                        synthon_sets.put(ki,new ArrayList<>());
                        for(Pair<String,String> fi : enamine_data.get(ri).get(ki) ) {
                            String idc = fi.getLeft();
                            String identifier = fi.getRight();
                            synthon_sets.get(ki).add( idc );
                            idc_to_identifier.put(idc,identifier);
                        }
                    }

                    space.addReaction(ri,synthon_sets,idc_to_identifier,null);

                    System.out.println("Done with thread "+fzi);
                    semaphore.release();
                }
            };

            try {
                all_space_creation_threads.add(ti);
                System.out.println("Try to start thread: "+zi+" / "+rxn_list.size());
                semaphore.acquire();
                ti.start();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // wait for all threads to finish:
        for(Thread ti : all_space_creation_threads) {
            try {
                ti.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        System.out.println("All space creation threads have finished!");
        System.out.println("reinit bittrees");

        // now reinit bitset trees!
        space.initAfterJavaDeserialization();
        space.reinitHelperMaps();
        space.reinitBitTree();

        return space;

    }


    /**
     * Usage:
     * Parameter:
     * 1. path to file containing the synthon space (must be new enamine real space format)
     * 2. name of the space (just determines output file name)
     * 3. mode ( one out of "mode_pfp" , "mode_skelspheres", "mode_ffp" or "mode_ppcore" )
     * 4. number of threads
     *
     * @param args
     */
    public static void main(String args[]) {
        String path = args[0];
        String space_name = args[1];
        String mode = args[2];
        int threads = -1;
        if(args.length > 3) {
            threads = Integer.parseInt(args[3]);
        }
        initREALSpace(path,space_name,mode,threads);
    }

    public static void initREALSpace(String path, String space_name, String mode, int threads) {
        initREALSpace(path,space_name,mode,threads,null);
    }

     public static void initREALSpace(String path, String space_name, String mode, int threads, String outputDirectory) {

        //process_files(path,space_name,mode,threads,4);
        File file = new File(path);

        SynthonSpace space = null;
        try {
            space = parseBuildingBlocks(file,mode,threads,3);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // then write to file:
        try {
            System.out.println("Start writing the output!");
            File f = null;
            if(outputDirectory == null) {
                f = new File(space_name +"_" +mode +".data");
            }
            else {
                f = new File(outputDirectory + File.separator + space_name + "_" + mode + ".data");
            }
            System.out.println("Into File: "+f.getName());
            if(false) {
                BufferedWriter out = new BufferedWriter(new FileWriter(f));
                out.write(space.serializeToString());
                out.flush();
                out.close();
            }
            if(true) {
                //ObjectOutputStream out = new ObjectOutputStream( new BufferedOutputStream( new FileOutputStream(f) ) );
                ObjectOutputStream out = new ObjectOutputStream( new GZIPOutputStream( new BufferedOutputStream( new FileOutputStream(f) ) ) );
                out.writeObject(space);
                out.flush();
                out.close();
                System.out.println("Done!");
            }
        }
        catch(IOException ex) {
            ex.printStackTrace();
        }


        // then create similarity search space:
         if(true) {
             space.initAfterJavaDeserialization();
             SynthonSimilaritySpace3 space3 = new SynthonSimilaritySpace3(space);
             space3.initFastSimilaritySearchers(threads,(x)->{});

             try {
                 System.out.println("Start writing the output!");
                 //File f = new File(space_name +"_" +mode +"_similarity3"+".data");
                 File f = null;
                 if(outputDirectory == null) {
                     f = new File(space_name +"_" +mode +"_similarity3" + ".data");
                 }
                 else {
                     f = new File(outputDirectory + File.separator + space_name + "_" + mode+ "_similarity3" + ".data");
                 }

                 System.out.println("Into File: "+f.getName());
                 if(false) {
                     BufferedWriter out = new BufferedWriter(new FileWriter(f));
                     out.write(space.serializeToString());
                     out.flush();
                     out.close();
                 }
                 if(true) {
                     ObjectOutputStream out = new ObjectOutputStream( new GZIPOutputStream(new BufferedOutputStream( new FileOutputStream(f) ) ) );
                     out.writeObject(space3);
                     out.flush();
                     out.close();
                     System.out.println("Done!");
                 }
             }
             catch(IOException ex) {
                 ex.printStackTrace();
             }

         }

    }



}
