package com.idorsia.research.chem.hyperspace.reparametrization;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.FastSubstructureSearcher;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class RealizedRxnReparametrization {

    private RxnReparametrization reparametrization;

    private SynthonSpace space;
    private String rxn_id;

    private Map<Integer, RealizedSynthonSetReparametrization.SynthonSetAnalysis> analysis;

    public RealizedRxnReparametrization(RxnReparametrization reparametrization, SynthonSpace space, String rxn_id) {
        this.reparametrization = reparametrization;
        this.space = space;
        this.rxn_id = rxn_id;
    }

    /**
     * first stage of analysis, applies the synthon set
     *
     * @param threads
     * @throws InterruptedException
     */
    public void processSynthonSets(int threads) throws InterruptedException {

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
            RealizedSynthonSetReparametrization sspi = new RealizedSynthonSetReparametrization(new StandardSynthonSetReparametrization(),frags);
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

    }

}
