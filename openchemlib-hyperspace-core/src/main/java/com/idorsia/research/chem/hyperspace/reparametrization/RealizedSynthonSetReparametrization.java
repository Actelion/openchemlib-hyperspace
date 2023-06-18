package com.idorsia.research.chem.hyperspace.reparametrization;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.util.Parallelizer;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Consumer;

public class RealizedSynthonSetReparametrization {

    public static class SynthonSetAnalysis {
        public final Map<String, String> pfrags;
        //public final List<SynthonSpace.FragId> pfrags;
        public final Map<String,List<SynthonSpace.FragId>>  pfrags_variants;
        public final List<Map.Entry<String,List<SynthonSpace.FragId>>> pfrags_sorted_filtered;
        public final Map<String,Integer> pfrags_numexitvectors;

        public SynthonSetAnalysis(Map<String, String> pfrags, Map<String,List<SynthonSpace.FragId>>  pfrags_variants ,
                                  List<Map.Entry<String,List<SynthonSpace.FragId>>> pfrags_sorted_filtered) {
            this.pfrags = pfrags;
            this.pfrags_variants = pfrags_variants;
            this.pfrags_sorted_filtered = pfrags_sorted_filtered;
            pfrags_numexitvectors = new HashMap<>();

            for(String idcode_i : pfrags.values()) {
                StereoMolecule m_processed = HyperspaceUtils.parseIDCode(idcode_i);//pfrags.get(idcode_i);
                m_processed.ensureHelperArrays(Molecule.cHelperNeighbours);
                int num_exitvectors = 0;
                for(int zi=0;zi<m_processed.getAtoms();zi++){ if(m_processed.getAtomicNo(zi)==90){ num_exitvectors++;} }
                if(!pfrags_numexitvectors.containsKey(idcode_i)){ pfrags_numexitvectors.put(idcode_i,num_exitvectors); }
            }
        }
    }

    private SynthonSetAnalysis analysisResult;

    private SynthonSetReparametrization reparametrization;
    private List<SynthonSpace.FragId>   synthonsUnprocessed;

    public RealizedSynthonSetReparametrization(SynthonSetReparametrization reparametrization, List<SynthonSpace.FragId> synthonsUnprocessed) {
        this.reparametrization = reparametrization;
        this.synthonsUnprocessed = synthonsUnprocessed;
    }

    public SynthonSetAnalysis getAnalysisResult() {
        return this.analysisResult;
    }

    /**
     *
     *
     * @param processed_synthon_cache contains map from synthon-idcode+reparametrization to processed synthon idcode
     * @param threads
     */
    public void computeAnalysis(Map<Pair<String,SynthonSetReparametrization>,String> processed_synthon_cache, int threads) throws InterruptedException {
        // frag id to idcode.
        Map<String,String>             pfrags = new ConcurrentHashMap<>();
        //List<SynthonSpace.FragId>      pfrags = Collections.synchronizedList(new ArrayList<>());
        Map<String,List<SynthonSpace.FragId>>       pfrags_variants = new ConcurrentHashMap<>();
        //Map<String,Integer>            pfrags_numexitvectors = new ConcurrentHashMap<>();

        Consumer<SynthonSpace.FragId> processSynthon = new Consumer<SynthonSpace.FragId>() {
            @Override
            public void accept(SynthonSpace.FragId fragId) {
                // visit all nodes up to distance x
                //int max_distance = 4;
                String idcode_i = null;
                Pair<String,SynthonSetReparametrization> rps = Pair.of(fragId.idcode,reparametrization);
                if(processed_synthon_cache.containsKey( rps ) ) {
                    idcode_i = processed_synthon_cache.get( rps );
                }
                else {
                    StereoMolecule fi = HyperspaceUtils.parseIDCode(fragId.idcode);
                    StereoMolecule m_processed = reparametrization.processSynthon(fi);
                    idcode_i = m_processed.getIDCode();
                    //pfrags.add(rps,idcode_i);
                }
                pfrags.put(fragId.fragment_id,idcode_i);
                //idcodes_of_pfrags.add(idcode_i);

                synchronized(pfrags_variants) {
                    if (!pfrags_variants.containsKey(idcode_i)) {
                        ArrayList<SynthonSpace.FragId> fvars = new ArrayList<>();
                        fvars.add(fragId);
                        pfrags_variants.put(idcode_i, fvars);
                    } else {
                        pfrags_variants.get(idcode_i).add(fragId);
                    }
                }
            }
        };
        Parallelizer.computeParallelBlocking(processSynthon,this.synthonsUnprocessed,threads);

        // sort by variants
        List<Map.Entry<String,List<SynthonSpace.FragId>>> pfrags_sorted = new ArrayList<>( pfrags_variants.entrySet() );
        // order descending:
        pfrags_sorted.sort( (x,y) -> -Integer.compare( x.getValue().size() , y.getValue().size() ) );

        //List<Map.Entry<String,Integer>> pfrags_sorted_filtered =  pfrags_sorted.stream().filter( vi -> vi.getValue().intValue() > 0 ).collect(Collectors.toList());
        List<Map.Entry<String,List<SynthonSpace.FragId>>> pfrags_sorted_filtered = pfrags_sorted;
        //System.out.println("Filtered: "+pfrags_sorted_filtered.size()+" of "+pfrags_sorted.size());

        this.analysisResult = new SynthonSetAnalysis(pfrags,pfrags_variants,pfrags_sorted_filtered);
    }

}
