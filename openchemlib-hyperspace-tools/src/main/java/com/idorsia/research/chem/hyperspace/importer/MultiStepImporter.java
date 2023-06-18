package com.idorsia.research.chem.hyperspace.importer;

import com.actelion.research.chem.IsomericSmilesCreator;
import org.apache.commons.lang3.tuple.Pair;
import tech.molecules.leet.chem.ChemUtils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class MultiStepImporter {


    /**
     * This is designed to parse the folder structure created by the demo of the
     * publication "Fully automatic .."
     *
     * Folder structure is:
     * --ReactionFolder
     *   --A
     *     --File with BBs A
     *     --virtual_bbs
     *       --rxnA_synthon_set_1.dwar
     *       --rxnA_synthon_set_2.dwar
     *       --rxnB_synthon_set_1.dwar
     *       --rxnB_synthon_set_2.dwar
     *   --B (same as A)
     *
     * This will be processed as follows:
     * We create the following combinations:
     * (FileWithBBs_A,FileWithBBs_B)
     * (FileWithBBs_A,any_virtualbb_rxn_from_B)
     * (FileWithBBs_B,any_virtualbb_rxn_from_A)
     *
     * @return
     */
    public static List<MultiStepSynthonReaction> parseReactions_01(String pathToReactionFolder, ImporterMode mode) throws Exception {
        // 1. create A
        File sfA = new File(pathToReactionFolder+File.separator+"A");
        File[] dwA_all = sfA.listFiles( new FilenameFilter(){
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".dwar");
            }
        });
        File dwA = dwA_all[0];
        ImporterSynthonReaction rA = null;
        try {
            rA = ImporterTool.importSynthonReaction( Collections.singletonList(dwA.getAbsolutePath()),"Enamine-ID");
        } catch (Exception e) {
            e.printStackTrace();
        }

        // 2. create B
        File sfB = new File(pathToReactionFolder+File.separator+"B");
        File[] dwB_all = sfB.listFiles( new FilenameFilter(){
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".dwar");
            }
        });
        File dwB = dwB_all[0];
        ImporterSynthonReaction rB = null;
        try {
            rB = ImporterTool.importSynthonReaction( Collections.singletonList(dwB.getAbsolutePath()),"Enamine-ID");
        } catch (Exception e) {
            e.printStackTrace();
        }


        ConcreteMultiStepSynthonReaction r_A = new ConcreteMultiStepSynthonReaction(rA,new ConnectorRemap());
        ConcreteMultiStepSynthonReaction r_B = new ConcreteMultiStepSynthonReaction(rB,new ConnectorRemap());

        List<MultiStepSynthonReaction> all_synthon_reactions = new ArrayList<>();

        MultiStepSynthonReaction r_AB = null;
        try {
            r_AB = ImporterTool.combineTwoSynthonReactions(r_A,r_B);
        } catch (Exception e) {
            e.printStackTrace();
        }

        all_synthon_reactions.add(r_AB);

        if(mode == ImporterMode.TwoSplit || mode == ImporterMode.ThreeSplit ) {
            // 3. create virtual_bb rxns A
            List<ConcreteMultiStepSynthonReaction> rxns_a = new ArrayList<>();
            List<ConcreteMultiStepSynthonReaction> rxns_b = new ArrayList<>();

            rxns_a = parseVirtualBBsFolder(pathToReactionFolder + File.separator + "A" + File.separator + "virtual_bbs");
            rxns_b = parseVirtualBBsFolder(pathToReactionFolder + File.separator + "B" + File.separator + "virtual_bbs");

            if( mode == ImporterMode.ThreeSplit ) {
                // 5. create all combinations r_a + r_B
                for (ConcreteMultiStepSynthonReaction ria : rxns_a) {
                    try {
                        for(ConcreteMultiStepSynthonReaction rib : rxns_b) {
                            MultiStepSynthonReaction r_xi = ImporterTool.combineTwoSynthonReactions(ria, rib);
                            all_synthon_reactions.add(r_xi);
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }

            // 4. create all combinations r_A + rxns_b
            for (ConcreteMultiStepSynthonReaction ri : rxns_b) {
                try {
                    MultiStepSynthonReaction r_xi = ImporterTool.combineTwoSynthonReactions(r_A, ri);
                    all_synthon_reactions.add(r_xi);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            // 5. create all combinations r_B + rxns_a

            for (ConcreteMultiStepSynthonReaction ri : rxns_a) {
                try {
                    MultiStepSynthonReaction r_xi = ImporterTool.combineTwoSynthonReactions(r_B, ri);
                    all_synthon_reactions.add(r_xi);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }



        return all_synthon_reactions;
    }

    /**
     * 1. searches for files with similar prefix that end in _1.dwar and _2.dwar and potentially _3.dwar
     *
     */
    public static List<ConcreteMultiStepSynthonReaction> parseVirtualBBsFolder(String folderPath) {
        //String folderPath = "/path/to/folder";

        File folder = new File(folderPath);
        File[] files = folder.listFiles((dir, name) -> name.endsWith(".dwar") && name.matches(".+_[0-9]+\\.dwar"));

        Map<String, List<File>> groups = new HashMap<>();

        for (File file : files) {
            String filename = file.getName();
            String prefix = filename.substring(0, filename.lastIndexOf('_'));

            if (groups.containsKey(prefix)) {
                groups.get(prefix).add(file);
            } else {
                List<File> fileList = new ArrayList<>();
                fileList.add(file);
                groups.put(prefix, fileList);
            }
        }

        // Print the file groups
        for (Map.Entry<String, List<File>> entry : groups.entrySet()) {
            String prefix = entry.getKey();
            List<File> fileList = entry.getValue();

            System.out.println("Group: " + prefix);
            for (File file : fileList) {
                System.out.println(file.getName());
            }
            System.out.println();
        }

        // Now create concrete multistep synthon reactions..
        List<ConcreteMultiStepSynthonReaction> mssr = new ArrayList<>();
        for(Map.Entry<String, List<File>> entry : groups.entrySet()) {
            //if(mssr.size()>=4) {break;}
            try {
                ImporterSynthonReaction ri = ImporterTool.importSynthonReaction(entry.getValue().stream().map(fi -> fi.getAbsolutePath()).collect(Collectors.toList()), "Enamine-ID");
                ConcreteMultiStepSynthonReaction cri = new ConcreteMultiStepSynthonReaction(ri,new ConnectorRemap());
                mssr.add(cri);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return mssr;
    }

    private enum ImporterMode { OneSplit , TwoSplit , ThreeSplit };

    /**
     * Arguments: 1. mode: "1split" or "2split"
     *            2. input-file-path
     *            3. output-file-name
     *
     * @param args
     */
    public static void main(String args[]) {
        //String path_a = "C:\\Temp\\virtual_spaces\\Benzimidazol_with_aldehyde";

        String str_mode   = args[0];
        String str_path_a = args[1];
        String str_output = args[2];

        ImporterMode mode = ImporterMode.OneSplit;
        if(str_mode.equalsIgnoreCase("1split")) { mode = ImporterMode.OneSplit;}
        else if(str_mode.equalsIgnoreCase("2split")) {mode = ImporterMode.TwoSplit;}
        else if(str_mode.equalsIgnoreCase("2split")) {mode = ImporterMode.ThreeSplit;}
        else {
            System.out.println("[INFO] unknown mode: "+str_mode);
            System.out.println("[INFO] unknown mode: -> fallback to 1split");
        }

        List<Pair<String,MultiStepSynthonReaction>> all_reactions = new ArrayList<>();
        //File folderCombLibraries = new File("C:\\Temp\\virtual_spaces\\CombinatorialLibraries");
        File folderCombLibraries = new File(str_path_a);
        for(File fi : folderCombLibraries.listFiles()) {
            List<MultiStepSynthonReaction> rxns_i = new ArrayList<>();
            try {
                rxns_i = parseReactions_01(fi.getAbsolutePath(),mode);
            } catch (Exception e) {
                e.printStackTrace();
            }
            for(int zi=0;zi<rxns_i.size();zi++) {
                all_reactions.add(Pair.of(String.format(fi.getName()+"_%03d",zi),rxns_i.get(zi)));
            }
        }

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(new File("output_a")));
            List<String> parts_header = new ArrayList<>();
            parts_header.add("smiles");
            parts_header.add("synthon");
            parts_header.add("synthonsetid");
            parts_header.add("reaction_id");
            out.write(String.join("\t",parts_header)+"\n");

            int counter_rxn = 0;
            for(Pair<String,MultiStepSynthonReaction> ri : all_reactions) {
                try {
                    exportReaction(ri.getRight(),ri.getLeft(),out);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                counter_rxn++;
            }
            out.flush();
            out.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        System.out.println("mkay..");
    }

    public static void exportReaction(MultiStepSynthonReaction ri, String name, BufferedWriter out) throws IOException {
        List<ImporterSynthonSet> all_synthon_sets = ri.getAllSynthonSets();
        int setcounter = 1;
        for(ImporterSynthonSet ssi : all_synthon_sets) {
            exportSynthonSet(ssi,name,setcounter,name+"_"+setcounter,out);
            setcounter++;
        }
        out.flush();
        return;
    }
    public static void exportSynthonSet(ImporterSynthonSet sset, String rxn_name, int synthonset_id , String synthon_set_name, BufferedWriter out) throws IOException {

        for(ImporterSynthon si : sset.synthons) {
            int synthoncounter = 0;
            try {
                String smiles = new IsomericSmilesCreator(si.getSynthon()).getSmiles();
                List<String> parts = new ArrayList<>();
                parts.add(smiles);
                parts.add(synthon_set_name+"_"+synthoncounter+"_"+si.getId());
                parts.add(""+synthonset_id);
                parts.add(rxn_name);
                out.write( String.join("\t" , parts) + "\n" );
                synthoncounter++;
            }
            catch(Exception ex) {
                ex.printStackTrace();
            }
        }
        out.flush();
    }


}
