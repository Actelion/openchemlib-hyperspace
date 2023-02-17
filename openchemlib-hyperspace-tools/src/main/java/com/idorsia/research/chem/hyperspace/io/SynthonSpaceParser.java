package com.idorsia.research.chem.hyperspace.io;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongPFP512;
import com.actelion.research.chem.io.DWARFileParser;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

public class SynthonSpaceParser {



    /**
     * Mode can be either "mode_pfp" , "mode_skelspheres", "mode_ffp" or "mode_ppcore"
     *
     * @param folder
     * @param mode
     * @return
     * @throws IOException
     */
    public static SynthonSpace parseBuildingBlocks(File folder, String mode, int threads, int max_number_of_synthon_sets) throws IOException {

        if(! (  mode.equals("PathFp") || mode.equals("FragFp") || mode.equals("mode_pfp") || mode.equals("mode_ppcore") || mode.equals("mode_ffp") || mode.equals("mode_skelspheres") )) {
            System.out.println("Mode not supported!!");
            return null;
        }


        CachedDescriptorProvider cdp = new CachedDescriptorProvider(mode);

        // Possibility to perform the "addReaction" in parallel
        //ExecutorService executor = Executors.newFixedThreadPool(6);
        ExecutorService executor = new ThreadPoolExecutor(
                threads, threads,
                //3,3,
                //8, 8,
                0L, TimeUnit.MILLISECONDS,
                //new LinkedBlockingQueue<Runnable>(64), // not too high, the queued tasks take a lot of memory!!
                new LinkedBlockingQueue<Runnable>(threads + 1), // not too high, the queued tasks take a lot of memory!!
                new ThreadPoolExecutor.CallerRunsPolicy()
        );

        ExecutorService executor_parse = new ThreadPoolExecutor(
                Math.max(1,threads/2), Math.max(1,threads/2),
                //3,3,
                //8, 8,
                0L, TimeUnit.MILLISECONDS,
                //new LinkedBlockingQueue<Runnable>(64), // not too high, the queued tasks take a lot of memory!!
                new LinkedBlockingQueue<Runnable>(4), // not too high, the queued tasks take a lot of memory!!
                new ThreadPoolExecutor.CallerRunsPolicy()
        );



        //String mode = "mode_pfp";
        //String mode = "mode_ppcore";

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
//        if(mode.equals("mode_skelspheres")) {
//            space.setFP(new DescriptorHandlerBinSkelSpheres2(),4096 , 512);
//        }
//        if(mode.equals("mode_skelspheres")) {
//            //space.setFP(new DescriptorHandlerAdapterByteToBinBuckets<>(new DescriptorHandlerSkeletonSpheres(),512,4),2048);
//            space.setFP(new DescriptorHandlerBinSkelSpheres2(),4096);
//        }


        //SmilesParser sp = new SmilesParser();
        int zfi = 0;

        List<String> rxn_ids = new ArrayList<>();

        for(File fi : folder.listFiles()) {
            //String splits[] = fi.getName().split("_");
            //String ri = splits[0]+"_"+splits[1]+"_";
            // NEW: split at last "_":
            int split_at = fi.getName().lastIndexOf("_");
            String rxn_id = fi.getName().substring(0,split_at);
            rxn_ids.add(rxn_id);
        }

        rxn_ids = rxn_ids.stream().distinct().collect(Collectors.toList());

        System.out.println();

        // rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("s_14") || si.startsWith("m_26") ||  si.startsWith("s_3") || si.startsWith("s_6") ).collect(Collectors.toList());
        //rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("m_24") ||  si.startsWith("s_3") ).collect(Collectors.toList());
        //rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("s_1") || si.startsWith("s_23") ||  si.startsWith("m_274") ).collect(Collectors.toList());
        //rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("m_24") ).collect(Collectors.toList());
        //rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("s_49a_") || si.startsWith( "s_6") ).collect(Collectors.toList());
        //rxn_ids = rxn_ids.stream().filter( si -> si.startsWith("s_1") || si.startsWith("s_3") || si.startsWith("s_4") || si.startsWith("m_") ).collect(Collectors.toList());

        long total_bytes     = Arrays.stream(folder.listFiles()).mapToLong(fi ->  fi.length() ).sum();
        long processed_bytes = 0;
        int total_files     = folder.listFiles().length;
        int processed_files = 0;

        for(String ri : rxn_ids) {
            List<String> possible_ids = new ArrayList<>();
            for(int zi=0;zi<=4;zi++) {
                possible_ids.add( (ri+"_"+zi+".dwar").toLowerCase() );
                possible_ids.add( (ri+"_"+zi+".smi").toLowerCase() );
            }

            File[] r_files = folder.listFiles(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return possible_ids.contains(name.toLowerCase());
                }
            });

            if(r_files.length>max_number_of_synthon_sets) {
                System.out.println("Skip > " + max_number_of_synthon_sets + "synthon sets for the moment: "+ri);
                continue;
            }

            Map<Integer,List<Object>> sorted = new HashMap<>();

            Map<String,String> idcode_to_identifier = new ConcurrentHashMap<>();

            //boolean problem_encountered = false;

            // this is the "name" / id of the synthon set!
            int synthon_set_cnt = 0;

            for (File fi : r_files) {
                System.out.println("File: "+fi);
                System.out.println("Progress: "+ (processed_files) + "/" + total_files + " byte ratio: "  + ( 1.0*processed_bytes/total_bytes ));

                List<Object> fragments_i = Collections.synchronizedList(new ArrayList<>());
                List<String> problems_encountered = Collections.synchronizedList(new ArrayList<>());


                if(fi.getName().endsWith(".smi")) {

                    String name = fi.getName();
                    String splits[] = name.split("_");

                    //ExecutorService executor_parse = Executors.newFixedThreadPool(6);
                    //ExecutorService executor_parse = Executors.newFixedThreadPool(threads);

                    List<Future> parse_futures = new ArrayList<>();

                    BufferedReader in = new BufferedReader(new FileReader(fi.getAbsolutePath()));
                    String line = null;
                    while ((line = in.readLine()) != null) {

                        String final_line = line;



//                        Runnable r_parse = new Runnable() {
//                            @Override
//                            public void run() {
//                                String final_synthon_before_parsing = null;
//                                SmilesParser sp = new SmilesParser();
//                                try {
//                                    StereoMolecule m = new StereoMolecule();
//                                    String l_splits[] = final_line.split("\\t");
//
//                                    // replace the uranium atoms %50,%51,%52,%53 by U,Np,Pu,Am
//                                    String smiles_synthon = l_splits[0];
//
//                                    //smiles_synthon = smiles_synthon.replace(".[U]%50", ".[U]%50");
//
//                                    if (false) {
//                                        smiles_synthon = smiles_synthon.replace(".[U]%51", ".[Np]%51");
//                                        smiles_synthon = smiles_synthon.replace(".[U]%52", ".[Pu]%52");
//                                        smiles_synthon = smiles_synthon.replace(".[U]%53", ".[Am]%53");
//                                    }
//                                    if (true) {
//                                        smiles_synthon = smiles_synthon.replace(".[U]%50", "");
//                                        smiles_synthon = smiles_synthon.replace(".[U]%51", "");
//                                        smiles_synthon = smiles_synthon.replace(".[U]%52", "");
//                                        smiles_synthon = smiles_synthon.replace(".[U]%53", "");
//
//                                        smiles_synthon = smiles_synthon.replace("%50", "([U])");
//                                        smiles_synthon = smiles_synthon.replace("%51", "([Np])");
//                                        smiles_synthon = smiles_synthon.replace("%52", "([Pu])");
//                                        smiles_synthon = smiles_synthon.replace("%53", "([Am])");
//                                    }
//
//                                    final_synthon_before_parsing = smiles_synthon;
//                                    sp.parse(m, smiles_synthon);
//
//                                    fragments_i.add(m);
//                                    String m_idcode = new String(""+m.getIDCode());
//                                    idcode_to_identifier.put("" + m_idcode, l_splits[1]);
//                                } catch (Exception ex) {
//                                    //ex.printStackTrace();
//                                    System.out.println("exception in rxn " + ri + " file=" + fi + " -> synthon " + final_line);
//                                    System.out.println("Tried to parse: " + final_synthon_before_parsing);
//
//                                    //problem_encountered = true;
//                                    problems_encountered.add("exception in rxn " + ri + " file=" + fi + " -> synthon " + final_line);
//                                }
//                            }
//                        };

                        ProcessSmilesLineRunnable runnable_smiles_line = new ProcessSmilesLineRunnable(final_line,ri,fi.getName(),
                                fragments_i,idcode_to_identifier,problems_encountered);

                        //executor_parse.execute(runnable_smiles_line);
                        //parse_futures.add( executor_parse.submit(runnable_smiles_line) );
                        runnable_smiles_line.run();
                    }

                    // wait for all parsing tasks to finish:
                    for(Future fti : parse_futures) {
                        try {
                            fti.get();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        } catch (ExecutionException e) {
                            e.printStackTrace();
                        }
                    }
                    System.out.println("File: "+fi+" --> Parsing complete, parsed "+fragments_i.size()+" molecules");

//                    executor_parse.shutdown();
//                    try {
//                        executor_parse.awaitTermination(10, TimeUnit.MINUTES);
//                    } catch (InterruptedException e) {
//                        e.printStackTrace();
//                    }

                }
                else if(fi.getName().endsWith(".dwar")) {
                    String name = fi.getName();
                    String splits[] = name.split("_");

                    //ExecutorService executor_parse = Executors.newFixedThreadPool(6);

                    //List<StereoMolecule> fragments_i = Collections.synchronizedList(new ArrayList<>());
                    //List<String> problems_encountered = Collections.synchronizedList(new ArrayList<>());

                    //BufferedReader in = new BufferedReader(new FileReader(fi.getAbsolutePath()));
                    DWARFileParser parse = new DWARFileParser(fi.getAbsolutePath());

                    //IDCodeParser icp = new IDCodeParser();

                    while(parse.next()) {
                        StereoMolecule mi = new StereoMolecule();
                        String idcode_i     = parse.getSpecialFieldData( parse.getSpecialFieldIndex("Synthon") );
                        //String identifier   = parse.getSpecialFieldData( parse.getSpecialFieldIndex("Building Block") );
                        String identifier   = parse.getSpecialFieldData( parse.getFieldIndex("SynthonId") );
                        //icp.parse(mi,idcode_i);
                        fragments_i.add(idcode_i);
                        idcode_to_identifier.put(idcode_i,identifier);
                    }

                    System.out.println("Parsed synthon set -> size= "+fragments_i.size());
                    // shrink the large sets a bit ;)
                    if(false) {
                        if (fragments_i.size() > 40000) {
                            System.out.println("Shrink to 40k structures..");
                            Collections.shuffle(fragments_i);
                            fragments_i = fragments_i.subList(0, 10000);
                        }
                    }

                    parse.close();
                }
                else {
                    System.out.println("cannot parse..");
                    continue;
                }



                //if(problem_encountered) {
                if(!problems_encountered.isEmpty()) {
                    System.out.println("\nProblems encountered in synthon set "+ri+"\n" + " successes: "+fragments_i.size()+" fails: "+problems_encountered.size());
                    System.out.println();
                    //System.out.println("Problem with synthon set "+ri+"   -> skip");
                    //continue;
                }


                // To validate the synthon set we sort the parsed synthons by their connector config and sort
                // out synthons that are different
                Map<String,List<Object>> connector_config_sorted_synthons = new ConcurrentHashMap<>();


                //for(StereoMolecule mi : fragments_i) {
                List<Future> jobs_process_synthons = new ArrayList<>();

                final List<Object> f_fragments_i = fragments_i;
                for(int zmi=0;zmi<fragments_i.size();zmi++) {

                    int fzmi = zmi;
                    Runnable r_process_synthon = new Runnable() {
                        @Override
                        public void run() {
                            StereoMolecule mi = new StereoMolecule();
                            IDCodeParser icp = new IDCodeParser();
                            icp.parse(mi,(String) f_fragments_i.get(fzmi));

                            mi.ensureHelperArrays(Molecule.cHelperCIP);
                            int cnt_u=0; int cnt_np = 0; int cnt_pu = 0; int cnt_am = 0;
                            for(int zi=0;zi<mi.getAtoms();zi++){
                                if(mi.getAtomicNo(zi)==92){cnt_u++;}
                                if(mi.getAtomicNo(zi)==93){cnt_np++;}
                                if(mi.getAtomicNo(zi)==94){cnt_pu++;}
                                if(mi.getAtomicNo(zi)==95){cnt_am++;}
                            }
                            String conn_config = "u"+cnt_u+"_np"+cnt_np+"_pu"+cnt_pu+"_am"+cnt_am;
                            if(!connector_config_sorted_synthons.containsKey(conn_config)){ connector_config_sorted_synthons.put(conn_config, Collections.synchronizedList(new ArrayList<>())); }
                            connector_config_sorted_synthons.get(conn_config).add(mi.getIDCode());
                        }
                    };
                    jobs_process_synthons.add( executor_parse.submit(r_process_synthon) );
                }

                // wait for jobs to finish:
                if(true) {
                    System.out.println("Wait for synthon processing jobs to finish:");
                    for (Future ji : jobs_process_synthons) {
                        try {
                            ji.get();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        } catch (ExecutionException e) {
                            e.printStackTrace();
                        }
                    }
                    System.out.println("Wait for synthon processing jobs to finish: done");
                }

                List<Object> fragments_i_validated = null;
                if(connector_config_sorted_synthons.keySet().size()>1) {
                    System.out.println("[WARN] Found multiple connector configs in synthon set "+(ri+synthon_set_cnt));
                    System.out.println(connector_config_sorted_synthons.keySet().toString());
                    int max_synthons = -1; String max_cc = null;
                    for(String cci : connector_config_sorted_synthons.keySet()) {
                        if(connector_config_sorted_synthons.get(cci).size()>max_synthons){ max_cc = cci; max_synthons = connector_config_sorted_synthons.get(cci).size();}
                    }
                    System.out.println("Stick with cc: "+max_cc);
                    fragments_i_validated = new ArrayList<>( connector_config_sorted_synthons.get(max_cc));
                }
                else {
                    fragments_i_validated = fragments_i;
                }



                //String str_int = splits[2].split("\\.")[0];
                //sorted.put(Integer.parseInt(str_int), fragments_i);
                sorted.put(synthon_set_cnt++,fragments_i_validated);


                if (zfi % 20 == 0) {
                    System.out.println();
                }
                zfi++;
                System.out.print(".");

                processed_bytes += fi.length();
                processed_files += 1;
            }

            if(sorted.values().stream().mapToInt( vi -> vi.size() ).anyMatch(svi -> svi==0) ) {
                System.out.println("Problem with rxn "+ri+" -> empty synthon set -> skip..");
                continue;
            }


            boolean merge_synthon_sets = false;

//            if(merge_synthon_sets) {
//                // merged sets will have different integer numbers (starting from 5 on)
//                // merge two smallest sets if we have four sets..
//                if (sorted.keySet().size() == 4) {
//                    sorted = tryMergeSynthonSets(sorted, idcode_to_identifier);
//                }
//            }

            // also important: compactify connectors.. i.e. if we do not have 92,93,94 ascending, make it so..

            boolean fix_connector_types = true;

            if(fix_connector_types) { // means, if the connectors are not ascending from 92, then we replace them by 92, 93 and so on.
                List<Integer> connectors = new ArrayList<>();
                for(int ki : sorted.keySet()) {
                    //StereoMolecule mi = sorted.get(ki).iterator().next();
                    StereoMolecule mi =  parseIDCode( (String) sorted.get(ki).iterator().next() );
                    Map<Integer,Integer> conn_pos = SynthonSpace.computeConnectorPositions(mi);
                    connectors.addAll( conn_pos.keySet() );
                }
                // check if ok:
                List<Integer> conn_unique = connectors.stream().distinct().sorted().collect(Collectors.toList());

                List<Integer> conn_required = new ArrayList<>();
                for(int zi=92;zi<92+conn_unique.size();zi++) { conn_required.add(zi); }
                if(!conn_unique.equals(conn_required)) {
                    System.out.println("Found incorrect connectors. Fix connectors!");
                    // map conn_unique[i] to conn_required[i]
                    Map<Integer,Integer> replacements = new HashMap<>();
                    for(int zi=0;zi<conn_unique.size();zi++) { replacements.put(conn_unique.get(zi),92+zi); }
                    for(int ki : sorted.keySet()) {
                        //for(StereoMolecule mi : sorted.get(ki)) {
                        List<Object> sorted_fixed_i = new ArrayList<>();
                        for(int zmi=0;zmi<sorted.get(ki).size();zmi++) {
                            StereoMolecule mi = parseIDCode( (String) sorted.get(ki).get(zmi) );
                            Map<Integer,Integer> conn_pos = SynthonSpace.computeConnectorPositions(mi);
                            //replace all:
                            for(int cti : conn_pos.keySet()) {
                                mi.setAtomicNo(conn_pos.get(cti) , replacements.get(cti));
                            }
                            //mi.ensureHelperArrays(Molecule.cHelperCIP);
                            sorted_fixed_i.add(mi.getIDCode());
                        }
                        sorted.put(ki,sorted_fixed_i);
                    }
                }
            }



            //space.addReaction(ri,sorted,idcode_to_identifier);

            Map<Integer,List<Object>> final_sorted = sorted;

//            Runnable runnable_add_reaction = new Runnable(){
//                @Override
//                public void run() {
//                    System.out.println("addReaction for Rxn: "+ri + " -> Start!");
//                    long t_start = System.currentTimeMillis();
//                    space.addReaction( ri,final_sorted,idcode_to_identifier );
//                    //System.gc();
//                    long t_end   = System.currentTimeMillis();
//                    System.out.println("addReaction for Rxn: "+ri + " -> Completed! Time= "+(t_end-t_start));
//                }
//            };
            AddReactionRunnable runnable_add_reaction = new AddReactionRunnable(cdp,space,ri,final_sorted,idcode_to_identifier);


            executor.execute( runnable_add_reaction );
            //runnable_add_reaction.run();

            sorted = null;
            final_sorted = null;

            System.gc();
//            if(processed_bytes > 5000000) {
//                break;
//            }
        }

        executor_parse.shutdown();

        try {
            executor_parse.awaitTermination(24, TimeUnit.HOURS);
            System.out.println("Executor 1 Shutdown! All Jobs finished!");
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        executor.shutdown();

        try {
            executor.awaitTermination(24, TimeUnit.HOURS);
            System.out.println("Executor 2 Shutdown! All Jobs finished!");
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // now reinit bst:
        space.reinitBitTree();
        return space;
    }


    public static StereoMolecule parseIDCode(String idc) {
        StereoMolecule mi_new = new StereoMolecule();
        IDCodeParser icp = new IDCodeParser();
        icp.parse(mi_new,idc);
        mi_new.ensureHelperArrays(Molecule.cHelperCIP);
        return mi_new;
    }

    public static class AddReactionRunnable implements Runnable {
        SynthonSpace space;
        String       rxn;
        Map<Integer,List<Object>> synthons_sorted;
        Map<String,String> idcode_to_identifier;

        CachedDescriptorProvider cdp;

        public AddReactionRunnable(CachedDescriptorProvider cdp, SynthonSpace space, String rxn, Map<Integer,List<Object>> synthons_sorted, Map<String,String> idcode_to_identifier) {
            this.space = space;
            this.rxn = rxn;
            this.synthons_sorted = synthons_sorted;
            this.idcode_to_identifier = new HashMap<>(idcode_to_identifier);
            this.cdp = cdp;
        }
        @Override
        public void run() {
            System.out.println("addReaction for Rxn: "+rxn + " -> Start!");
            long t_start = System.currentTimeMillis();
            this.space.addReaction( this.rxn,this.synthons_sorted,idcode_to_identifier , cdp );
            //System.gc();
            long t_end   = System.currentTimeMillis();
            System.out.println("addReaction for Rxn: "+rxn + " -> Completed! Time= "+(t_end-t_start));

            this.synthons_sorted.clear();
            this.synthons_sorted = null;
        }
    }


    public static class ProcessSmilesLineRunnable implements Runnable {

        private String final_line;
        private String ri;
        private String fi;
        private Map<String,String> idcode_to_identifier;
        private List<Object> fragments_i;
        private List<String> problems_encountered;

        public ProcessSmilesLineRunnable(String line, String rxn, String file, List<Object> frags_i, Map<String,String> idcode_to_identifier, List<String> problems_encountered) {
            this.final_line = line;
            this.ri = rxn;
            this.fi = file;
            this.fragments_i = frags_i;
            this.idcode_to_identifier = idcode_to_identifier;
            this.problems_encountered = problems_encountered;
        }

            @Override
            public void run() {
            String final_synthon_before_parsing = null;

            //SmilesParser sp = new SmilesParser();
            SmilesParser sp = getSmilesParser();

            try {
                StereoMolecule m = new StereoMolecule();
                String l_splits[] = final_line.split("\\t");

                // replace the uranium atoms %50,%51,%52,%53 by U,Np,Pu,Am
                String smiles_synthon = l_splits[0];

                //smiles_synthon = smiles_synthon.replace(".[U]%50", ".[U]%50");

                if (false) {
                    smiles_synthon = smiles_synthon.replace(".[U]%51", ".[Np]%51");
                    smiles_synthon = smiles_synthon.replace(".[U]%52", ".[Pu]%52");
                    smiles_synthon = smiles_synthon.replace(".[U]%53", ".[Am]%53");
                }
                if (true) {
                    smiles_synthon = smiles_synthon.replace(".[U]%50", "");
                    smiles_synthon = smiles_synthon.replace(".[U]%51", "");
                    smiles_synthon = smiles_synthon.replace(".[U]%52", "");
                    smiles_synthon = smiles_synthon.replace(".[U]%53", "");

                    smiles_synthon = smiles_synthon.replace("%50", "([U])");
                    smiles_synthon = smiles_synthon.replace("%51", "([Np])");
                    smiles_synthon = smiles_synthon.replace("%52", "([Pu])");
                    smiles_synthon = smiles_synthon.replace("%53", "([Am])");
                }

                final_synthon_before_parsing = smiles_synthon;
                sp.parse(m, smiles_synthon);

                StereoMolecule m2 = new StereoMolecule(m);
                m2.ensureHelperArrays(Molecule.cHelperCIP);

                //fragments_i.add(m);
                fragments_i.add(m2.getIDCode());

                String m_idcode = new String(""+m.getIDCode());
                idcode_to_identifier.put("" + m_idcode, l_splits[1]);
            } catch (Exception ex) {
                //ex.printStackTrace();
                System.out.println("exception in rxn " + ri + " file=" + fi + " -> synthon " + final_line);
                System.out.println("Tried to parse: " + final_synthon_before_parsing);

                //problem_encountered = true;
                problems_encountered.add("exception in rxn " + ri + " file=" + fi + " -> synthon " + final_line);
            }

            // set references to null:
            this.fragments_i = null;
        }
    }

    /**
     * NOTE!! idcode_to_identifier will be modified, i.e. after this it
     * contains additional entries for the merged synthons!
     *
     * @param sorted
     * @param idcode_to_identifier
     * @return
     */
    public static Map<Integer, List<StereoMolecule>> tryMergeSynthonSets(Map<Integer, List<StereoMolecule>> sorted, Map<String,String> idcode_to_identifier) {
        List<Map.Entry<Integer,List<StereoMolecule>>> size_sorted = sorted.entrySet().stream().sorted( (x, y) -> Integer.compare( x.getValue().size() , y.getValue().size() ) ).collect(Collectors.toList());


        List<StereoMolecule> set_a = size_sorted.get(0).getValue();
        List<StereoMolecule> set_b = size_sorted.get(1).getValue();
        List<StereoMolecule> set_c = size_sorted.get(2).getValue();
        Integer key_a = size_sorted.get(0).getKey();
        Integer key_b = size_sorted.get(1).getKey();
        Integer key_c = size_sorted.get(2).getKey();



        // verify that a and b have a common connector:
        Set<Integer> connectors_a = SynthonSpace.computeConnectorSet(set_a.get(0));
        Set<Integer> connectors_b = SynthonSpace.computeConnectorSet(set_b.get(0));
        Set<Integer> connectors_c = SynthonSpace.computeConnectorSet(set_c.get(0));

        connectors_a.retainAll(connectors_b);

        if( !connectors_a.isEmpty() ) {
            System.out.println("Merge synthon sets");
            // we can merge :)
            List<StereoMolecule> mergedSet = mergeSets(set_a,set_b,idcode_to_identifier);

            // remove the two original sets and put the new
            sorted.remove(key_a);
            sorted.remove(key_b);
            sorted.put(5,mergedSet);
        }
        else {
            // then, we try to merge a and c
            System.out.println("Cannot merge two smallest synthon sets.. try to merge a and c");
            connectors_a = SynthonSpace.computeConnectorSet(set_a.get(0));
            connectors_c = SynthonSpace.computeConnectorSet(set_c.get(0));
            connectors_a.retainAll(connectors_c);
            if( !connectors_a.isEmpty() ) {
                System.out.println("Merge synthon sets");
                // we can merge :)
                List<StereoMolecule> mergedSet = mergeSets(set_a,set_c,idcode_to_identifier);

                // remove the two original sets and put the new
                sorted.remove(key_a);
                sorted.remove(key_c);
                sorted.put(5,mergedSet);
            }
            else {
                System.out.println("Cannot merge a and c.. :(");
            }
        }
        return sorted;
    }

    /**
     * Merges sets a and b. It also adds entries for the combined bb ids to the idcode_to_identifier map.
     *
     * @param set_a
     * @param set_b
     * @param idcode_to_identifier
     * @return
     */
    public static List<StereoMolecule> mergeSets(List<StereoMolecule> set_a, List<StereoMolecule> set_b, Map<String,String> idcode_to_identifier) {

        Random r = new Random();
        List<StereoMolecule> merged_set = new ArrayList<>();

        for(StereoMolecule ma : set_a) {
            String a_id = idcode_to_identifier.get(ma.getIDCode());
            if(a_id==null) { a_id = "unknown_a_"+r.nextInt(10000000);}
            for(StereoMolecule mb : set_b) {
                String b_id = idcode_to_identifier.get(mb.getIDCode());
                if(b_id==null) { b_id = "unknown_b_"+r.nextInt(10000000);}
                String merged_id = a_id+":::"+b_id;

                List<StereoMolecule> to_merge = new ArrayList<>();
                to_merge.add(ma);
                to_merge.add(mb);
                System.out.println(HyperspaceUtils.idcodeToSmiles(ma.getIDCode()) + " merged with "+HyperspaceUtils.idcodeToSmiles(mb.getIDCode()));

                try {
                    StereoMolecule merged = SynthonAssembler.assembleSynthons(to_merge);
                    System.out.println(" --> " + HyperspaceUtils.idcodeToSmiles(merged.getIDCode()));
                    merged_set.add(merged);
                    idcode_to_identifier.put(merged.getIDCode(), merged_id);
                }
                catch(Exception ex) {
                    ex.printStackTrace();
                    System.out.println("Failed to merge synthons.. :(");
                }
            }
        }

        return merged_set;
     }


    /**
     * Usage:
     * Parameter:
     * 1. path to folder containing the synthon space (must be either SMILES or DWAR format)
     * 2. name of the space (just determines output file name)
     * 3. mode
     * 4. number of threads
     *
     *
     *
     *
     *
      * @param args
     */
    public static void main(String args[]) {
        String path = args[0];
        String output_path = args[1];
        String space_name = args[2];
        String mode = args[3];
        int threads = Integer.parseInt(args[4]);

        //process_files(path,space_name,mode,threads,4);
        try {
            process_files(path, output_path, space_name, mode, threads, 3);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main_(String args[]) {
        int num_threads = 3;
        try {
            process_files("C:\\datasets\\enamine_real_space_2_input_molecules_baby", "/home/liphath1/hyperspace_out", "enamine_pfp_space_8_a_baby", "mode_ffp", num_threads, 4);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    public static void main_b(String args[]) {
        int num_threads = 3;
        String path_to_idovcs = "C:\\datasets\\IdorsiaVCS";

        try {
            process_files(path_to_idovcs, "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_ppcore", num_threads, 4);
            //process_files( "/home/liphath1/enamine_real_space_2/input/molecules" , "enamine_pfp_space_8_a" ,"mode_ppcore");

            process_files(path_to_idovcs, "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_ffp", num_threads, 4);
            //process_files( "/home/liphath1/enamine_real_space_2/input/molecules" , "enamine_pfp_space_8_a" ,"mode_ffp");

            process_files(path_to_idovcs, "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_pfp", num_threads, 4);
            //process_files( "/home/liphath1/enamine_real_space_2/input/molecules" , "enamine_pfp_space_8_a" ,"mode_pfp");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

     public static void main_a(String args[]) {
         int num_threads = 8;

         try {

             process_files("/mnt/org/Chem/DD_Chem/Chem_Tech/Personal Folders/wahljo1/IdorsiaVCS", "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_ppcore", num_threads, 4);
             process_files("/home/liphath1/enamine_real_space_2/input/molecules", "/home/liphath1/hyperspace_out", "enamine_pfp_space_8_a", "mode_ppcore", num_threads, 4);

             process_files("/mnt/org/Chem/DD_Chem/Chem_Tech/Personal Folders/wahljo1/IdorsiaVCS", "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_ffp", num_threads, 4);
             process_files("/home/liphath1/enamine_real_space_2/input/molecules", "/home/liphath1/hyperspace_out", "enamine_pfp_space_8_a", "mode_ffp", num_threads, 4);

             process_files("/mnt/org/Chem/DD_Chem/Chem_Tech/Personal Folders/wahljo1/IdorsiaVCS", "/home/liphath1/hyperspace_out", "IdoVCS_1_a", "mode_pfp", num_threads, 4);
             process_files("/home/liphath1/enamine_real_space_2/input/molecules", "/home/liphath1/hyperspace_out", "enamine_pfp_space_8_a", "mode_pfp", num_threads, 4);
         } catch (IOException e) {
             e.printStackTrace();
         }

     }


    /**
     *
     * mode currently can be "mode_pfp" , "mode_skelspheres", "mode_ffp" or "mode_ppcore"
     *
     * @param path
     * @param space_name
     * @param mode
     * @param threads
     */
    public static void process_files(String path, String output_path, String space_name, String mode, int threads, int max_number_of_synthon_sets) throws IOException {

        SmilesParser p = new SmilesParser();
        StereoMolecule m = new StereoMolecule();

        // 1. check that output path exists:
        File f_output_path = new File(output_path);
        if( !f_output_path.isDirectory() ) {
            throw new IOException("Output directory does not exist");
        }
        if( !f_output_path.canWrite()) {
            throw new IOException("Cannot write to output directory");
        }

        try {
            p.parse(m,"COCC(C)N=C(S%50)N%51N=%52.[U]%50.[U]%51.[U]%52");
            System.out.println(m.getAtoms());
            for(int zi=0;zi<m.getAtoms();zi++) {
                System.out.println("M="+m.getAtomicNo(zi)+ " Col="+m.getAtomColor(zi));
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

        try {

            //String mode = "mode_pfp";
            //String mode = "mode_ppcore";
            //String mode = "mode_ffp";
            //String mode = "mode_ppcore";
            //String mode = "mode_skelspheres";
            //String mode = "mode_ffp";
            //SynthonSpace space = parseBuildingBlocks(new File("C:\\datasets\\2019-12_REAL_update_Idorsia\\input\\molecules"),mode);
            //SynthonSpace space = parseBuildingBlocks(new File("/home/liphath1/enamine_real_space/input/molecules"),mode);
            System.out.println("Parse building blocks!");

            //SynthonSpace space = parseBuildingBlocks(new File("/home/liphath1/enamine_real_space_2/input/molecules"),mode);
            //String space_name = "enamine_pfp_space_8_a";
            //SynthonSpace space = parseBuildingBlocks(new File("/mnt/org/Chem/DD_Chem/Chem_Tech/Personal Folders/wahljo1/IdorsiaVCS"),mode);
            //String space_name = "IdoVCS_1_a";
            SynthonSpace space = parseBuildingBlocks(new File(path),mode,threads,max_number_of_synthon_sets);

            try {
                System.out.println("Start writing the output!");
                //File f = new File("enamine_pfp_space_4_b.data");
                //File f = new File(space_name +"_" +mode +"."+ (new Random()).nextInt(100000) +".data");
                File f = new File(output_path+ File.separator +space_name +"_" +mode +".data");
                System.out.println("Into File: "+f.getName());
                if(false) {
                    BufferedWriter out = new BufferedWriter(new FileWriter(f));
                    out.write(space.serializeToString());
                    out.flush();
                    out.close();
                }
                if(true) {
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
                    File f = new File(space_name +"_" +mode +"_similarity3"+".data");
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
        catch(IOException ex) {
            ex.printStackTrace();
        }
    }


    private static Map<Thread,SmilesParser> mSmilesParsersPerThread = new HashMap<>();

    public synchronized static SmilesParser getSmilesParser() {
        Thread ti = Thread.currentThread();
        if(!mSmilesParsersPerThread.containsKey(ti)) {
            SmilesParser spi = new SmilesParser();
            mSmilesParsersPerThread.put(ti,spi);
        }
        return mSmilesParsersPerThread.get(ti);
    }

}
