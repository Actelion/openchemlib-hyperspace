package com.idorsia.research.chem.hyperspace.reparametrization;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.idorsia.research.chem.hyperspace.*;
import org.apache.commons.io.input.CountingInputStream;
import org.apache.commons.lang3.tuple.Pair;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

public class RealizedRxnReparametrization {

    private RxnReparametrization reparametrization;

    private SynthonSpace space;
    private String rxn_id;

    private Map<Integer, RealizedSynthonSetReparametrization.SynthonSetAnalysis> analysis;

    public static class RxnScaffold {

    }

    public RealizedRxnReparametrization(RxnReparametrization reparametrization, SynthonSpace space, String rxn_id) {
        this.reparametrization = reparametrization;
        this.space = space;
        this.rxn_id = rxn_id;
    }

    public static interface ScaffoldEnumerationCriterion {
        public boolean considerScaffold(List<Pair<String,List<SynthonSpace.FragId>>> scaffoldParts);
    }

    public List<SynthonSpaceScaffold> enumerateCombinatorialScaffolds(ScaffoldEnumerationCriterion scaffoldEnumerationCriterion) {
        List<SynthonSpaceScaffold> scaffolds = new ArrayList<>();
        List<List<Pair<String,List<SynthonSpace.FragId>>>> scaffoldParts = new ArrayList<>();
        for( RealizedSynthonSetReparametrization.SynthonSetAnalysis ai : analysis.values()) {
            scaffoldParts.add( ai.pfrags_sorted_filtered.stream().map( xi -> Pair.of(xi.getKey(),xi.getValue())).collect(Collectors.toList()) );
        }
        List<List<Pair<String,List<SynthonSpace.FragId>>>> tuples = generateTuples( scaffoldParts );

        for( List<Pair<String,List<SynthonSpace.FragId>>> scaffold_i : tuples ) {
            boolean consider = true;
            if( scaffoldEnumerationCriterion!=null ) {
                consider = scaffoldEnumerationCriterion.considerScaffold(scaffold_i);
            }
            if(!consider) {
                System.out.println("Skip scaffold..");
                continue;
            }
            else {
                System.out.println("Construct scaffold.. size="+scaffold_i.stream().mapToLong(xi -> xi.getRight().size()).reduce((x,y)->x*y)+" , "+scaffold_i.stream().map(xi -> ""+xi.getRight().size()).reduce((x,y)->x+","+y));
                List<StereoMolecule> scaffoldParts_i = scaffold_i.stream().map( xi -> HyperspaceUtils.parseIDCode(xi.getKey())).collect(Collectors.toList());
                StereoMolecule reactedScaffold = SynthonReactor.react(scaffoldParts_i);
                reactedScaffold.ensureHelperArrays(Molecule.cHelperCIP);
                if( reactedScaffold.getFragments().length > 1 ) {
                    System.out.println("[WARN] assembly incomplete..");
                }
                String rxn_i = scaffold_i.get(0).getRight().get(0).rxn_id;
                if(scaffold_i.size()==2) {
                    SynthonSpaceScaffold scaffoldData_i = new SynthonSpaceScaffold(rxn_i, reactedScaffold.getIDCode(),
                            scaffoldParts_i.get(0).getIDCode(),scaffoldParts_i.get(1).getIDCode(),"",
                            scaffold_i.get(0).getRight().stream().map(xi -> xi.idcode).collect(Collectors.toList()),
                            scaffold_i.get(1).getRight().stream().map(xi -> xi.idcode).collect(Collectors.toList()),
                            new ArrayList<>());
                    scaffolds.add(scaffoldData_i);
                }
                else if (scaffold_i.size()==3) {
                    SynthonSpaceScaffold scaffoldData_i = new SynthonSpaceScaffold(rxn_i, reactedScaffold.getIDCode(),
                            scaffoldParts_i.get(0).getIDCode(),scaffoldParts_i.get(1).getIDCode(),scaffoldParts_i.get(2).getIDCode(),
                            scaffold_i.get(0).getRight().stream().map(xi -> xi.idcode).collect(Collectors.toList()),
                            scaffold_i.get(1).getRight().stream().map(xi -> xi.idcode).collect(Collectors.toList()),
                            scaffold_i.get(2).getRight().stream().map(xi -> xi.idcode).collect(Collectors.toList()));
                    scaffolds.add(scaffoldData_i);
                }
            }
        }

        return scaffolds;
    }

    private static <T> List<List<T>>  generateTuples(List<List<T>> lists) {
        List<List<T>> result = new ArrayList<>();
        generateTuplesHelper(lists, 0, new ArrayList<>(), result);
        return result;
    }

    private static <T> void generateTuplesHelper(List<List<T>> lists, int currentIndex, List<T> currentTuple, List<List<T>> result) {
        if (currentIndex == lists.size()) {
            // Base case: Add the current tuple to the result
            result.add(new ArrayList<>(currentTuple));
            return;
        }

        List<T> currentList = lists.get(currentIndex);

        for (int i = 0; i < currentList.size(); i++) {
            currentTuple.add(currentList.get(i));
            generateTuplesHelper(lists, currentIndex + 1, currentTuple, result);
            currentTuple.remove(currentTuple.size() - 1); // Backtrack
        }
    }

    /**
     * first stage of analysis, applies the synthon set reparametrization and sorts
     * synthons accordingly.
     *
     * @param threads
     * @throws InterruptedException
     */
    public void processSynthonSets(int threads) throws InterruptedException {
        this.analysis = new HashMap<>();

        FastSubstructureSearcher fss = null;
        CachedDescriptorProvider cdh = new CachedDescriptorProvider("FragFp");


        //List<SynthonSpaceScaffold> sampled_scaffolds     = new ArrayList<>();

        Random r = new Random();
        // cluster synthon sets:
        Map<Pair<String,SynthonSetReparametrization>, String> processed_synthon_cache_a = new ConcurrentHashMap<>();
        //Map<String,StereoMolecule> processed_synthon_cache_b = new HashMap<>();
        //Map<String,StereoMolecule> processed_synthon_cache_c = new HashMap<>();
        //int cnt_rxn = 0;

        //List<String> all_rxn_ids = space.connector_fps_sorted_by_fragtype.keySet().stream().map( fi -> fi.rxn_id ).distinct().collect(Collectors.toList());
        //List<String> all_rxn_ids = new ArrayList<>();
        //all_rxn_ids.add("m_3caa");
        String rxn = rxn_id;

        List<Integer> num_frags = new ArrayList<>();
        List<Integer> num_frags_processed = new ArrayList<>();
        List<List<StereoMolecule>> to_assemble        = new ArrayList<>();
        List<List<List<String>>> to_assemble_var_info = new ArrayList<>();
        //for( int fzi : space.fragment_map.get(rxn).keySet() ) {
        for( int fzi : space.getFragTypes(rxn).keySet() ) {

            List<SynthonSpace.FragId> frags = space.getSynthonSet(rxn,fzi);
            num_frags.add(frags.size());
            RealizedSynthonSetReparametrization sspi = new RealizedSynthonSetReparametrization( reparametrization.getSynthonSetReparametrizatino(fzi) , frags );//new RealizedSynthonSetReparametrization(new StandardSynthonSetReparametrization(),frags);
            sspi.computeAnalysis(processed_synthon_cache_a,threads);
            RealizedSynthonSetReparametrization.SynthonSetAnalysis analysis = sspi.getAnalysisResult();
            //SynthonSetAnalysis analysis_b = process_synthon_set(processed_synthon_cache_b,frags,4,false);
            //SynthonSetAnalysis analysis_c = process_synthon_set(processed_synthon_cache_c,frags,3,true);

            if(true) {
                //writeSynthonSetAnalysis("synset_analysis_NEW_"+rxn+"_"+fzi+"_a.csv", analysis);
            }
            num_frags_processed.add(analysis.pfrags.keySet().size());
            this.analysis.put(fzi,analysis);
            //int limit = 32;
            //if(space.getFragTypes(rxn).keySet().size()>=3) { limit = 24;}
            // take the limit:
            //List<Map.Entry<String,List<String>>> selection = analysis.pfrags_sorted_filtered.subList(0, Math.min(limit,analysis.pfrags_sorted_filtered.size()) );
            //to_assemble.add( new ArrayList<>( selection.stream().map(ei -> analysis.pfrags.get(ei.getKey())  ).collect(Collectors.toList()) ) );
            //to_assemble_var_info.add( new ArrayList<>( selection.stream().map(ei -> ei.getValue() ).collect(Collectors.toList())) );
        }
        System.out.println("RXN "+rxn);
        System.out.println("Size            "+ num_frags.stream().map( vi -> vi.toString() ).collect(Collectors.joining(",")) );
        System.out.println("Size proc.      "+ num_frags_processed.stream().map( vi -> vi.toString() ).collect(Collectors.joining(",")) );
        System.out.println("Products        "+ num_frags.stream().mapToLong( vi -> vi.longValue() ).reduce( (a,b) -> a*b ) );
        System.out.println("Products proc.  "+ num_frags_processed.stream().mapToLong( vi -> vi.longValue() ).reduce( (a,b) -> a*b ) );
    }


    public static void main(String args[]) {
        //String inputfile = "v2_s2_32k_FragFp.data";
        String inputfile = "s2_small_01_FragFp.data";

        ObjectInputStream in = null;
        SynthonSpace space = null;
        try {
            FileInputStream file_in = new FileInputStream(inputfile);
            //counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
            //counting_in_stream = new CountingInputStream( new GZIPInputStream( new BufferedInputStream(file_in)));
            CountingInputStream counting_in_stream = new CountingInputStream( new BufferedInputStream(file_in));
            in = new ObjectInputStream( new GZIPInputStream( counting_in_stream ) );
            SynthonSpace space_a = null;
            space_a = (SynthonSpace) in.readObject();
            space_a.initAfterJavaDeserialization();
            space = space_a;
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        System.out.println("Loaded space: "+space.getSpaceInfoString());
        List<SynthonSpaceScaffold> scaffolds = new ArrayList<>();

        for(String rxn_id : space.getRxnIds()) {
            SynthonSetReparametrization synthonSetReparametrization = new StandardSynthonSetReparametrization(true, 3);
            RxnReparametrization standardReparametrization = new RxnReparametrization(space, rxn_id, synthonSetReparametrization);
            RealizedRxnReparametrization reparametrization = new RealizedRxnReparametrization(standardReparametrization, space, rxn_id);
            try {
                reparametrization.processSynthonSets(12);
                List<SynthonSpaceScaffold> scaffolds_i = reparametrization.enumerateCombinatorialScaffolds(new ScaffoldEnumerationCriterion() {
                    @Override
                    public boolean considerScaffold(List<Pair<String, List<SynthonSpace.FragId>>> scaffoldParts) {
                        return scaffoldParts.stream().filter(xi -> xi.getRight().size() > 8).collect(Collectors.toList()).size() >= 2;
                        //return true;
                    }
                });
                scaffolds.addAll(scaffolds_i);
                System.out.println("mkay");
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
        System.out.println("Analysis: total scaffolds returned: "+scaffolds.size());
        for(SynthonSpaceScaffold sci : scaffolds) {
            System.out.println(sci.toCSV_short());
        }
    }

}
