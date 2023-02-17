package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import org.apache.commons.lang3.tuple.Pair;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



public class HyperspaceSubstructureSearchProcess extends AbstractHyperspaceSearchProcess<SynthonSpace.CombinatorialHit> implements AbstractHyperspaceProcess.HasProgress {

    HyperspaceSubstructureSearch.SubstructureSearchConfiguration config;

    SynthonSpace space;
    CachedDescriptorProvider cdp;


    List<SynthonSpace.CombinatorialHit> results = new ArrayList<>();
    //HyperspaceSearchGUI gui;
    //JPanel p_results;

    public HyperspaceSubstructureSearchProcess(HyperspaceSubstructureSearch search_provider, SynthonSpace space, CachedDescriptorProvider cdp) {
        super(search_provider);
        //this.gui = gui;
        //this.p_results = results_panel;

        this.space = space;
        this.cdp    = cdp;
    }

    @Override
    public String getName() {
        return "Hyperspace Substructure Search";
    }

    //@Override
    public List<StereoMolecule> getQueryStructures() {
        return config.getQueryMolecules();
    }


    @Override
    public void setSearchConfiguration(AbstractSearchProvider.SearchConfiguration config) {
        this.config = (HyperspaceSubstructureSearch.SubstructureSearchConfiguration) config;
    }

    public HyperspaceSubstructureSearch.SubstructureSearchConfiguration getSearchConfiguration() {
        return this.config;
    }

    @Override
    public void startSearchAsync() {
        super.startSearchAsync();
        setProcessStatus(ProcessStatus.COMPUTING);

        SubstructureSearchThread search_thread = new SubstructureSearchThread(this.getSearchConfiguration(), this.config.getQueryMolecules().get(0) );
        search_thread.start();
    }




    public class SubstructureSearchThread extends Thread {

        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config;
        StereoMolecule mi;

        public SubstructureSearchThread(HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config, StereoMolecule mi) {
            this.search_config = search_config;
            this.mi = mi;
        }

        @Override
        public void run() {

            System.out.println("[INFO] search: "+this.search_config.getQuery().getIDCode());

            int number_of_threads = search_config.getNumberOfThreads();
            // check if this gets overruled by the max number of threads that we can use
            if( getSearchProvider().getMaxNumberOfThreads() > 0 ) {
                number_of_threads = Math.min(number_of_threads, getSearchProvider().getMaxNumberOfThreads());
            }
            if(number_of_threads<0) { number_of_threads = Runtime.getRuntime().availableProcessors(); }

            int max_combi_hits = this.search_config.getMaxNumberOfCombinatorialHits();



            List<StereoMolecule> allQueryMolecules = new ArrayList<>(); // TODO: use..
            StereoMolecule mq = this.search_config.getQuery();
            mq.ensureHelperArrays(Molecule.cHelperCIP);

            if(mq.getFragments().length>1) {
                System.out.println("search failed, must be single fragment");
                setProcessStatusMessage("Failed, more than 1 fragment..");
                setProcessStatus(ProcessStatus.FAILED);
                setProgressAndFireUpdate(1.0);
                return;
            }

            try {
                allQueryMolecules.addAll(SubstructureSearchHelper.expandBridgedSearches(mq));
            } catch (Exception e) {
                setProcessStatusMessage("Error while expadning bridges");
            }

            System.out.println("To process: "+allQueryMolecules.size()+" query structures");

            for(int zs=0;zs<allQueryMolecules.size();zs++) {
                StereoMolecule qi = allQueryMolecules.get(zs);
                System.out.println("Start processing query structure "+zs+" -> "+qi.getIDCode());

                List<SynthonSpace.CombinatorialHit> hits_no_splits = new ArrayList<>();
                List<SynthonSpace.CombinatorialHit> hits_1s = new ArrayList<>();
                List<SynthonSpace.CombinatorialHit> hits_2s = new ArrayList<>();
                List<SynthonSpace.CombinatorialHit> hits_3s = new ArrayList<>();
                List<SynthonSpace.CombinatorialHit> hits_4s = new ArrayList<>();

                setProcessStatusMessage("Screen BBs");

                //List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space, cdp, mi, number_of_threads);
                List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space, cdp, qi, number_of_threads);
                List<SynthonSpace.CombinatorialHit> hits_bbs_as_combinatorial_hits = new ArrayList<>();

                Set<String> discovered_rxns = new HashSet<>();

                if (!hits_bbs.isEmpty()) {
                    //BitSet fp_molecule = cdp.getFP_cached(mi);
                    BitSet fp_molecule = cdp.getFP_cached(qi);

                    // 1. sort expanded hits by frag type.. (then we report all hits from the same frag type in the same combinatorial hit..)
                    Map<SynthonSpace.FragType, List<SynthonSpace.ExpandedHit>> sorted_bb_hits = new HashMap<>();
                    for (SynthonSpace.ExpandedHit hi : hits_bbs) {
                        Integer hit_frag = hi.hit_fragments.keySet().iterator().next();
                        SynthonSpace.FragId hit_frag_id = hi.hit_fragments.get(hit_frag);

                        SynthonSpace.FragType hit_ft = new SynthonSpace.FragType(hi.rxn, hit_frag_id.frag);
                        if (!sorted_bb_hits.containsKey(hit_ft)) {
                            sorted_bb_hits.put(hit_ft, new ArrayList<>());
                        }
                        sorted_bb_hits.get(hit_ft).add(hi);
                    }

                    for (SynthonSpace.FragType fti_a : sorted_bb_hits.keySet()) {

                        List<SynthonSpace.FragId> all_matching_fragids = new ArrayList<>();
                        for (SynthonSpace.ExpandedHit hi : sorted_bb_hits.get(fti_a)) {
                            all_matching_fragids.add(hi.hit_fragments.get(fti_a.frag));
                        }

                        //SynthonSpace.CombinatorialHit chi = new SynthonSpace.CombinatorialHit();

                        Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> frag_map = new HashMap<>();
                        Map<Integer, Pair<SynthonSpace.FragType, BitSet>> fp_mapping = new HashMap<>();

                        StereoMolecule sri_fragments[] = new StereoMolecule[space.getFragTypes(fti_a.rxn_id).size()];
                        //StereoMolecule mi_a = new StereoMolecule(mi);
                        StereoMolecule mi_a = new StereoMolecule(qi);
                        mi_a.ensureHelperArrays(Molecule.cHelperCIP);
                        for (int szi = 0; szi < sri_fragments.length; szi++) {
                            if (szi == 0) {
                                sri_fragments[szi] = mi_a;
                            } else {
                                sri_fragments[szi] = new StereoMolecule();
                            }
                        }
                        SynthonShredder.SplitResult sri = new SynthonShredder.SplitResult(sri_fragments, new ArrayList<>(), "", new ArrayList<>());

                        // loop over frag types and fill the required maps:
                        int other_splits_cnt = 1;
                        for (Integer ftii : space.getFragTypes(fti_a.rxn_id).keySet()) {
                            SynthonSpace.FragType fti = space.getFragTypes(fti_a.rxn_id).get(ftii);
                            if (ftii == fti_a.frag) {
                                fp_mapping.put(0, Pair.of(fti_a, fp_molecule));
                                frag_map.put(fti_a, all_matching_fragids);
                            }
//                            else {
//                                // put all fragments into combinatorial hit:
//                                frag_map.put(fti, new ArrayList<>(space.getSynthonSet(fti_a.rxn_id, ftii)));
//                                // put empty bitset into fp_mapping
//                                fp_mapping.put(other_splits_cnt, Pair.of(fti, new BitSet(space.getBits())));
//                                other_splits_cnt++;
//                            }
                        }
                        SynthonSpace.CombinatorialHit hit_bb_ci = new SynthonSpace.CombinatorialHit(fti_a.rxn_id, frag_map, sri, fp_mapping);
                        hits_bbs_as_combinatorial_hits.add(hit_bb_ci);
                    }
                }


                if (hits_bbs_as_combinatorial_hits.size() > max_combi_hits) {

                } else {
                    discovered_rxns.addAll(hits_bbs_as_combinatorial_hits.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
                    setProcessStatusMessage("Screen 1 Split");
                    setProgressAndFireUpdate(0.05);
                    List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics = Collections.synchronizedList(new ArrayList<>());
                    //hits_1s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 1, 2, discovered_rxns, 200, number_of_threads, output_statistics);
                    hits_1s = space.findExpandedHits_withConnProximityMatching(space, cdp, qi, 1, 2, discovered_rxns, 200, number_of_threads, output_statistics);
                    discovered_rxns.addAll(hits_1s.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
                }
                if (hits_bbs_as_combinatorial_hits.size() + hits_1s.size() > max_combi_hits) {

                } else {
                    setProcessStatusMessage("Screen 2 Split");
                    setProgressAndFireUpdate(0.1);
                    List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics = Collections.synchronizedList(new ArrayList<>());
                    //hits_2s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 2, 3, discovered_rxns, 200, number_of_threads, output_statistics);
                    hits_2s = space.findExpandedHits_withConnProximityMatching(space, cdp, qi, 2, 3, discovered_rxns, 200, number_of_threads, output_statistics);
                    discovered_rxns.addAll(hits_2s.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
                }

                if (hits_bbs_as_combinatorial_hits.size() + hits_1s.size() + hits_2s.size() > max_combi_hits) {

                } else {
                    setProcessStatusMessage("Screen 3 Split");
                    setProgressAndFireUpdate(0.25);
                    List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics = Collections.synchronizedList(new ArrayList<>());
                    //hits_3s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 3, 4, discovered_rxns, 200, number_of_threads, output_statistics);
                    hits_3s = space.findExpandedHits_withConnProximityMatching(space, cdp, qi, 3, 4, discovered_rxns, 200, number_of_threads, output_statistics);
                    discovered_rxns.addAll(hits_3s.stream().map(hi -> hi.rxn).collect(Collectors.toList()));
                }
                //if(hits_3s.isEmpty()) {
                //    hits_4s = space.findExpandedHits_withConnProximityMatching(space,cdp,mi,4,4,200,1);
                //}

                System.out.println("Hits (bb): " + hits_bbs_as_combinatorial_hits.size());
                System.out.println("Hits (1s): " + hits_1s.size());
                System.out.println("Hits (2s): " + hits_2s.size());
                System.out.println("Hits (3s): " + hits_3s.size());
                System.out.println("Hits (4s): " + hits_4s.size());

                List<SynthonSpace.CombinatorialHit> all_hits = new ArrayList<>();
                all_hits.addAll(hits_bbs_as_combinatorial_hits);
                all_hits.addAll(hits_1s);
                all_hits.addAll(hits_2s);
                all_hits.addAll(hits_3s);
                all_hits.addAll(hits_4s);

                // now fill incomplete mappings:
                for(int zi=0;zi<all_hits.size();zi++) {
                    SynthonSpace.CombinatorialHit ci = all_hits.get(zi);
                    // check if all synthon sets of the synthon reaction are matched. If not, add the fragments to the combinatorial hit
                    List<SynthonSpace.FragType> synthsets_matched = new ArrayList<>(ci.mapping.values().stream().map( xi -> (xi!=null)?xi.getLeft(): null).collect(Collectors.toList()));
                    List<SynthonSpace.FragType> extra_synthsets = space.getFragTypes(ci.rxn).values().stream().collect(Collectors.toList());
                    extra_synthsets.removeAll(synthsets_matched);
                    if(!extra_synthsets.isEmpty()) {
                        for(SynthonSpace.FragType fti : extra_synthsets) {
                            System.out.println("Fill incomplete mapping for set: "+fti.toString());
                            ci.hit_fragments.put(fti,space.getSynthonSet(fti.rxn_id,fti.frag));
                        }
                    }
                }

                //results = all_hits;
                results.addAll(all_hits);

                if (results.size() > search_config.getMaxNumberOfCombinatorialHits()) {
                    results = new ArrayList<>(results.subList(0, search_config.getMaxNumberOfCombinatorialHits()));
                    break;
                }
            }
            setProcessStatusMessage("");
            setProcessStatus(ProcessStatus.DONE);
            setProgressAndFireUpdate(1.0);

        }
    }

    @Override
    public List<SynthonSpace.CombinatorialHit> getSearchResults() {
        return this.results;
    }

    private double progress;

    private void setProgressAndFireUpdate(double progress) {
        this.progress = progress;
        this.fireProcessStatusChanged();
    }

    @Override
    public double getProgress() {
        return this.progress;
    }
}
