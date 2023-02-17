package com.idorsia.research.chem.hyperspace.services;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerSkeletonSpheres;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyDescription;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import com.idorsia.research.chem.hyperspace.cli.InputStructureParser;
import com.idorsia.research.chem.hyperspace.descriptor.SimilarityEvaluationHelper;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import com.idorsia.research.chem.hyperspace.util.ConnectedPairSampler;
import com.idorsia.research.chem.hyperspace.util.HitExpander;
import com.idorsia.research.chem.hyperspace.util.RelevantSubgraphSampler;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class HyperspaceRelevantSubgraphSearchService extends HSService {

    public static class HyperspaceRelevantSubgraphSearchServiceConfig extends HSService.Config {
        @JsonPropertyDescription("config for hyperspace sss")
        @JsonProperty("sssConfig")
        private
        HyperspaceSubstructureSearchService.HyperspaceSubstructureSearchServiceConfig sssConfig;

        @JsonPropertyDescription("minimum size of considered fragments")
        @JsonProperty("minFragSize")
        private
        int minFragSize;

        @JsonPropertyDescription("maximum size of considered fragments")
        @JsonProperty("maxFragSize")
        private
        int maxFragSize;

        @JsonPropertyDescription("minimum ffp bits of considered fragments")
        @JsonProperty("minFFPBits")
        private
        int minFFPBits;

        @JsonPropertyDescription("maximum ffp bits of considered fragments")
        @JsonProperty("maxFFPBits")
        private
        int maxFFPBits;

        @JsonPropertyDescription("number of representative fragments considered")
        @JsonProperty("numClusters")
        private
        int numClusters;


        public HyperspaceRelevantSubgraphSearchServiceConfig() {

        }

        public HyperspaceRelevantSubgraphSearchServiceConfig(HyperspaceSubstructureSearchService.HyperspaceSubstructureSearchServiceConfig sssConfig, int minFragSize, int maxFragSize, int minFFPBits, int maxFFPBits, int numClusters) {
            this.setSssConfig(sssConfig);
            this.setMinFragSize(minFragSize);
            this.setMaxFragSize(maxFragSize);
            this.setMinFFPBits(minFFPBits);
            this.setMaxFFPBits(maxFFPBits);
            this.setNumClusters(numClusters);
        }


        public HyperspaceSubstructureSearchService.HyperspaceSubstructureSearchServiceConfig getSssConfig() {
            return sssConfig;
        }

        public void setSssConfig(HyperspaceSubstructureSearchService.HyperspaceSubstructureSearchServiceConfig sssConfig) {
            this.sssConfig = sssConfig;
        }

        public int getMinFragSize() {
            return minFragSize;
        }

        public void setMinFragSize(int minFragSize) {
            this.minFragSize = minFragSize;
        }

        public int getMaxFragSize() {
            return maxFragSize;
        }

        public void setMaxFragSize(int maxFragSize) {
            this.maxFragSize = maxFragSize;
        }

        public int getMinFFPBits() {
            return minFFPBits;
        }

        public void setMinFFPBits(int minFFPBits) {
            this.minFFPBits = minFFPBits;
        }

        public int getMaxFFPBits() {
            return maxFFPBits;
        }

        public void setMaxFFPBits(int maxFFPBits) {
            this.maxFFPBits = maxFFPBits;
        }

        public int getNumClusters() {
            return numClusters;
        }

        public void setNumClusters(int numClusters) {
            this.numClusters = numClusters;
        }
    }


    @Override
    public void runService(Config conf_a) {
        HyperspaceRelevantSubgraphSearchServiceConfig conf = (HyperspaceRelevantSubgraphSearchServiceConfig) conf_a;


        // Part 1.1 check if we have to set threads:
        if( conf.getSssConfig().getThreads() <= 0 ) {
            conf.getSssConfig().setThreads( Runtime.getRuntime().availableProcessors() );
        }

        // Part 1.2 a parse input structures
        List<StereoMolecule> structures = InputStructureParser.defaultParseInput(conf.getSssConfig().getInput());


        // Part 2. configure and run the combi hit search
        HyperspaceSubstructureSearchService.ArmedSubstructureSearch sss = HyperspaceSubstructureSearchService.prepareSubstructureSearch(conf.getSssConfig());

        HyperspaceSubstructureSearch search = sss.getHyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = sss.getHyperspaceSubstructureSearchConfig();

        // here we store all hits.. (query mol, performed sss query-idc, sxh)
        List<Triple<String, String, HitExpander.SimpleExpandedHit>> all_collected_hits = new ArrayList<>();


        Map<String, SimilarityEvaluationHelper.EvaluatedSimilarities> all_similarity_data = new HashMap<>();

        for(StereoMolecule mq : structures) {

            System.out.println("[INFO] Start with query: " + mq.getIDCode());

            // here we store all hits.. (result idcode, query-idc, sxh)
            List<Triple<String, String, HitExpander.SimpleExpandedHit>> all_collected_hits_i = new ArrayList<>();


            //StereoMolecule mq = HyperspaceUtils.parseIDCode( conf.getSssConfig().getInput() );
            RelevantSubgraphSampler rss = new RelevantSubgraphSampler(mq);
            List<RelevantSubgraphSampler.RelevantSubgraph> queries =  rss.sampleSubgraphs_A( conf.getMinFragSize() , conf.getMaxFragSize() ,
                    conf.getMinFFPBits() , conf.getMaxFFPBits() , conf.getNumClusters() );


            System.out.println("[INFO] connected pair subqueries: "+queries.size());


            for (int zi = 0; zi < queries.size(); zi++) {
                StereoMolecule qi = queries.get(zi).getFrag();
                    List<HitExpander.SimpleExpandedHit> expandedHits = HyperspaceSubstructureSearchService.computeExpandedHits(sss, conf.getSssConfig(), qi );
                    // add files to collected hits..
                    String query_idc = qi.getIDCode();
                    for (HitExpander.SimpleExpandedHit hix : expandedHits) {
                        all_collected_hits_i.add(Triple.of(mq.getIDCode(), query_idc, hix));
                    }
                    System.out.println("[INFO] Added " + expandedHits.size() + " hits to total hits for split query: " + qi.getIDCode());
            }



            System.out.println("[INFO] All results computed, compute descriptors and similarities..");

            List<String> hits_idcodes = new ArrayList<>( all_collected_hits_i.stream().map(ti -> ti.getRight().getIDCodeAssembled()).distinct().collect(Collectors.toList()) );
            System.out.println("[INFO] Hit structures for descriptor calc: "+hits_idcodes.size());


            List<DescriptorHandler> dhs = new ArrayList<>();
            dhs.add(new DescriptorHandlerSkeletonSpheres());
            dhs.add(new DescriptorHandlerLongFFP512());

            SimilarityEvaluationHelper.EvaluatedSimilarities similarity_data = SimilarityEvaluationHelper.evaluateSimilarities(Collections.singletonList(mq),
                    hits_idcodes, dhs,conf.getSssConfig().getThreads(),100,1000);

            all_collected_hits.addAll(all_collected_hits_i);
            all_similarity_data.put(mq.getIDCode(),similarity_data);
        }

        System.out.println("[INFO] All results computed, write output..");

        // 2. put into dwar:
        List<Pair<String, HitExpander.SimpleExpandedHit>> all_collected_hits_a = new ArrayList<>();
        List<Triple<String, Boolean, List<String>>> splits_a = new ArrayList<>();
        splits_a.add(Triple.of("QueryMol",true,new ArrayList<>()));
        splits_a.add(Triple.of("Sim_SkelSpheres",false,new ArrayList<>()));
        splits_a.add(Triple.of("Sim_FFP",false,new ArrayList<>()));
        for(int zi=0;zi<all_collected_hits.size();zi++) {
            SimilarityEvaluationHelper.EvaluatedSimilarities similarity_data = all_similarity_data.get(all_collected_hits.get(zi).getLeft());
            all_collected_hits_a.add(Pair.of(all_collected_hits.get(zi).getLeft(),all_collected_hits.get(zi).getRight()));
            splits_a.get(0).getRight().add(all_collected_hits.get(zi).getMiddle());
            splits_a.get(1).getRight().add( ""+ similarity_data.getSimilarities().get("SkelSpheres").get(all_collected_hits.get(zi).getLeft()).get(all_collected_hits.get(zi).getRight().getIDCodeAssembled()) );
            splits_a.get(2).getRight().add( ""+ similarity_data.getSimilarities().get("FragFp").get(all_collected_hits.get(zi).getLeft()).get(all_collected_hits.get(zi).getRight().getIDCodeAssembled()) );
        }


        HitExpander.exportExpandedHitsWithQueriesToDWAR(conf.getSssConfig().getOutput(),all_collected_hits_a,splits_a);

        System.out.println("[INFO] All results computed, write output.. complete!");
        System.out.println("[INFO] All done, thx for using hyperspace!");


    }


    public static void main(String args[]) {

        String json_file = args[0];
        String json_config = "";
        HyperspaceRelevantSubgraphSearchServiceConfig conf_a = null;
        try (FileReader fi_in = new FileReader(json_file)) {
            json_config = IOUtils.toString(new BufferedReader(fi_in));
            conf_a =
                    (new ObjectMapper()).reader().readValue(json_config, HyperspaceRelevantSubgraphSearchServiceConfig.class );

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        HyperspaceRelevantSubgraphSearchService service =
                new HyperspaceRelevantSubgraphSearchService();

        service.runService(conf_a);
    }


}
