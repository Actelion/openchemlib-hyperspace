package com.idorsia.research.chem.hyperspace.services;

import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.cli.InputStructureParser;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import com.idorsia.research.chem.hyperspace.sar.SARQueryFeatureGenerator;
import com.idorsia.research.chem.hyperspace.util.HitExpander;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class HyperspaceSARDirectedHitExpansionService extends HSService {

    public static class SARDHitExpansionConfig extends HyperspaceSubstructureSearchService.HyperspaceSubstructureSearchServiceConfig {

    }

    @Override
    public void runService(Config conf_a) {
        SARDHitExpansionConfig conf = (SARDHitExpansionConfig) conf_a;

        // Part 0.1 check if we have to set threads:
        if( conf.getThreads() <= 0 ) {
            conf.setThreads( Runtime.getRuntime().availableProcessors() );
        }

        // Part 0.2 a parse input structures
        List<StereoMolecule> structures = InputStructureParser.defaultParseInput(conf.getInput());


        // Part 1. configure and run the combi hit search
        HyperspaceSubstructureSearchService.ArmedSubstructureSearch sss = HyperspaceSubstructureSearchService.prepareSubstructureSearch(conf);

        HyperspaceSubstructureSearch search = sss.getHyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = sss.getHyperspaceSubstructureSearchConfig();

        // here we store all hits..
        List<Triple<String, String, HitExpander.SimpleExpandedHit>> all_collected_hits = new ArrayList<>();

        // Part 2. loop over structures:
        for(StereoMolecule mi : structures) {
            // compute sar fragments:
            List<SARQueryFeatureGenerator.SingleSplitResult> splits = SARQueryFeatureGenerator.generateSingleSplitFragments(mi);

            for(SARQueryFeatureGenerator.SingleSplitResult split : splits) {
                List<HitExpander.SimpleExpandedHit> expandedHits = HyperspaceSubstructureSearchService.computeExpandedHits(sss, conf, split.split_query);
                // add files to collected hits..
                String query_idc = mi.getIDCode();
                for (HitExpander.SimpleExpandedHit hix : expandedHits) {
                    all_collected_hits.add(Triple.of(query_idc, split.split_query.getIDCode() , hix));
                }
                System.out.println("[INFO] Added "+expandedHits.size()+" hits to total hits for split query: "+split.split_query);
            }
        }

        System.out.println("[INFO] All results computed, write output..");

        // 2. put into dwar:
        List<Pair<String, HitExpander.SimpleExpandedHit>> all_collected_hits_a = new ArrayList<>();
        List<Triple<String, Boolean, List<String>>> splits_a = new ArrayList<>();
        splits_a.add(Triple.of("SARD QF",true,new ArrayList<>()));
        for(int zi=0;zi<all_collected_hits.size();zi++) {
            all_collected_hits_a.add(Pair.of(all_collected_hits.get(zi).getLeft(),all_collected_hits.get(zi).getRight()));
            splits_a.get(0).getRight().add(all_collected_hits.get(zi).getMiddle());
        }

        HitExpander.exportExpandedHitsWithQueriesToDWAR(conf.getOutput(),all_collected_hits_a,splits_a);

        System.out.println("[INFO] All results computed, write output.. complete!");

        System.out.println("[INFO] All done, thx for using hyperspace!");

    }

    public static void main(String args[]) {

        String json_file = args[0];
        String json_config = "";
        HyperspaceSARDirectedHitExpansionService.SARDHitExpansionConfig conf_a = null;
        try (FileReader fi_in = new FileReader(json_file)) {
            json_config = IOUtils.toString(new BufferedReader(fi_in));
            conf_a =
                    (new ObjectMapper()).reader().readValue(json_config, HyperspaceSARDirectedHitExpansionService.SARDHitExpansionConfig.class);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        HyperspaceSARDirectedHitExpansionService service =
                new HyperspaceSARDirectedHitExpansionService();

        service.runService(conf_a);
    }

}
