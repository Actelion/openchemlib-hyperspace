package com.idorsia.research.chem.hyperspace.services;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.hyperspace.SimpleCombinatorialHit;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyDescription;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.cli.InputStructureParser;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearchForDW;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearchProcess;
import com.idorsia.research.chem.hyperspace.util.HitExpander;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.tuple.Pair;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Semaphore;
import java.util.stream.Collectors;

public class HyperspaceSubstructureSearchService extends HSService {

    public static class HyperspaceSubstructureSearchServiceConfig extends HSService.Config {

        @JsonPropertyDescription("file containing input structures. ")
        @JsonProperty("input")
        private String input;
        @JsonPropertyDescription("file containing output data. ")
        @JsonProperty("output")
        private String output;

        @JsonPropertyDescription("HS: file containing the space to load")
        @JsonProperty("space")
        private String space;
        @JsonPropertyDescription("HS: max number of combinatorial hits per query structure [default 32]")
        @JsonProperty("maxCombiHits")
        private int maxCombiHits = 32;
        @JsonPropertyDescription("HS: max number of expanded hits enumerated per query structure [default -1]")
        @JsonProperty("maxExpandedHits")
        private int maxExpandedHits = -1;
        @JsonPropertyDescription("HS: max number of expanded hits enumerated per combinatorial hit [default 1000]")
        @JsonProperty("maxExpandedHitsPerCombiHit")
        private int maxExpandedHitsPerCombiHit = 1000;
        @JsonPropertyDescription("HS: min number of ffp bits in query [default 24]")
        @JsonProperty("minFFPQueryBits")
        private int minFFPQueryBits = 24;
        @JsonPropertyDescription("HS: max time in ms for (combinatorial) query [default 60000]")
        @JsonProperty("maxQueryTimeMillis")
        private int maxQueryTimeMillis = 60000;
        @JsonPropertyDescription("HS: number of threads used, -1 to ask system for number of cores [default -1]")
        @JsonProperty("threads")
        private int threads = -1;
        @JsonPropertyDescription("HS: max number of splits applied [default 3]")
        @JsonProperty("maxSplits")
        private int maxSplits = 3;

        public String getInput() {
            return input;
        }

        public void setInput(String input) {
            this.input = input;
        }

        public String getSpace() {
            return space;
        }

        public void setSpace(String space) {
            this.space = space;
        }



        public int getMinFFPQueryBits() {
            return minFFPQueryBits;
        }

        public void setMinFFPQueryBits(int minFFPQueryBits) {
            this.minFFPQueryBits = minFFPQueryBits;
        }

        public int getThreads() {
            return threads;
        }

        public void setThreads(int threads) {
            this.threads = threads;
        }

        public int getMaxSplits() {
            return maxSplits;
        }

        public void setMaxSplits(int maxSplits) {
            this.maxSplits = maxSplits;
        }

        public int getMaxCombiHits() {
            return maxCombiHits;
        }

        public void setMaxCombiHits(int maxCombiHits) {
            this.maxCombiHits = maxCombiHits;
        }

        public int getMaxExpandedHits() {
            return maxExpandedHits;
        }

        public void setMaxExpandedHits(int maxExpandedHits) {
            this.maxExpandedHits = maxExpandedHits;
        }

        public String getOutput() {
            return output;
        }

        public void setOutput(String output) {
            this.output = output;
        }

        public int getMaxQueryTimeMillis() {
            return maxQueryTimeMillis;
        }

        public void setMaxQueryTimeMillis(int maxQueryTimeMillis) {
            this.maxQueryTimeMillis = maxQueryTimeMillis;
        }

        public void setMaxExpandedHitsPerCombiHit(int max) {
            this.maxExpandedHitsPerCombiHit = max;
        }

        public int getMaxExpandedHitsPerCombiHit() {
            return this.maxExpandedHitsPerCombiHit;
        }
    }

    @Override
    public void runService(Config conf_a) {

        HyperspaceSubstructureSearchServiceConfig conf = (HyperspaceSubstructureSearchServiceConfig) conf_a;

        // Part 0.1 check if we have to set threads:
        if( conf.getThreads() <= 0 ) {
            conf.setThreads( Runtime.getRuntime().availableProcessors() );
        }

        // Part 0.2 a parse input structures
        List<StereoMolecule> structures = InputStructureParser.defaultParseInput(conf.getInput());


        // Part 1. configure and run the combi hit search
        ArmedSubstructureSearch sss = prepareSubstructureSearch(conf);

        HyperspaceSubstructureSearch search = sss.getHyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = sss.getHyperspaceSubstructureSearchConfig();

        // here we store all hits..
        List<Pair<String, HitExpander.SimpleExpandedHit>> all_collected_hits = new ArrayList<>();

        // Part 2. loop over structures:
        for(StereoMolecule mi : structures) {
            List<HitExpander.SimpleExpandedHit> expandedHits = computeExpandedHits(sss,conf,mi);
            // add files to collected hits..
            String query_idc = mi.getIDCode();
            for(HitExpander.SimpleExpandedHit hix : expandedHits) {
                all_collected_hits.add(Pair.of(query_idc,hix));
            }

            System.out.println("[INFO] Added "+expandedHits.size()+" hits to total hits for query: "+query_idc);
        }

        System.out.println("[INFO] All results computed, write output..");

        // 2. put into dwar:
        HitExpander.exportExpandedHitsWithQueriesToDWAR(conf.getOutput(),all_collected_hits);

        System.out.println("[INFO] All results computed, write output.. complete!");

        System.out.println("[INFO] All done, thx for using hyperspace!");

    }

    public static ArmedSubstructureSearch prepareSubstructureSearch(HyperspaceSubstructureSearchServiceConfig conf) {
        HyperspaceSubstructureSearch search = new HyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.InitializationConfig config_init = new HyperspaceSubstructureSearch.InitializationConfig("[none]",conf.getSpace(),conf.getThreads());
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = new HyperspaceSubstructureSearch.SubstructureSearchConfiguration();

        search_config.setNumberOfThreads(conf.getThreads());
        search_config.setMaxSplits(conf.getMaxSplits());
        search_config.setMaxNumberOfCombinatorialHits(conf.getMaxCombiHits());

        search.setConfigurationAndGUI(config_init,null);
        AbstractHyperspaceProcess proc_init = search.startInitialization();

        System.out.println("[INFO] Loading space");

        proc_init.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
            private long ts_last_loading_output = -1;
            @Override
            public void processStatusChanged() {
                // nothing
            }
        });

        long loading_info_interval  = 2000;
        TimerTask reportLoadingProgress = new TimerTask(){
            @Override
            public void run() {
                if(proc_init instanceof AbstractHyperspaceProcess.HasProgress) {
                    //if ( ts_last_loading_output < 0 || System.currentTimeMillis() - ts_last_loading_output > loading_info_interval) {
                    double progress_loading = ((AbstractHyperspaceProcess.HasProgress) proc_init).getProgress();
                    System.out.println("[INFO] Loading space: " + String.format( "%.2f" , progress_loading * 100 ) + "%" );
                    //ts_last_loading_output = System.currentTimeMillis();
                    //}
                }
            }
        };
        Timer timer_loadSpace = new Timer();
        timer_loadSpace.schedule(reportLoadingProgress,100,loading_info_interval);
        try {
            proc_init.waitUntilDoneOrFailed(3600000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        timer_loadSpace.cancel();
        System.out.println("[INFO] Loading space complete");
        return new ArmedSubstructureSearch(search,search_config);
    }

    public static class ArmedSubstructureSearch {
        private HyperspaceSubstructureSearch search;
        private HyperspaceSubstructureSearch.SubstructureSearchConfiguration searchConfig;
        public HyperspaceSubstructureSearch getHyperspaceSubstructureSearch() {
            return this.search;
        }
        public HyperspaceSubstructureSearch.SubstructureSearchConfiguration getHyperspaceSubstructureSearchConfig() {
            return this.searchConfig;
        }
        public ArmedSubstructureSearch(HyperspaceSubstructureSearch search, HyperspaceSubstructureSearch.SubstructureSearchConfiguration config) {
            this.search = search;
            this.searchConfig = config;
        }
    }

    public static List<HitExpander.SimpleExpandedHit> computeExpandedHits(ArmedSubstructureSearch armedSearch, HyperspaceSubstructureSearchServiceConfig conf, StereoMolecule mi) {
        System.out.println("[INFO] Start processing structure: "+mi.getIDCode());

        HyperspaceSubstructureSearch search = armedSearch.getHyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = armedSearch.getHyperspaceSubstructureSearchConfig();

        // count ffp bits of query:
        long[] ffpi  = DescriptorHandlerLongFFP512.getDefaultInstance().createDescriptor(mi) ;
        int ffp_bits = BitSet.valueOf(ffpi).cardinality();
        System.out.println("[INFO] FFP bits in query: " + ffp_bits);
        if(ffp_bits < conf.getMinFFPQueryBits() ) {
            System.out.println("[INFO] Too few ffp bits -> Skip!");
            return new ArrayList<>();
        }


        System.out.println("[INFO] Start hyperspace sss for structure: "+mi.getIDCode());
        search_config.setQuery(mi);
        try {

            AbstractHyperspaceSearchProcess<SynthonSpace.CombinatorialHit> p_search = search.runSearch(search_config);
            p_search.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
                @Override
                public void processStatusChanged() {
                }
            });

            // wait until done:
            p_search.waitUntilDoneOrFailed(conf.getMaxQueryTimeMillis());

            if(! p_search.getProcessStatus().equals(AbstractHyperspaceProcess.ProcessStatus.DONE)) {
                System.out.println("[WARN] bad process status after hyperspace search.. Skip..");
                return new ArrayList<>();
            }

            List<SynthonSpace.CombinatorialHit> combi_hits = p_search.getSearchResults();
            System.out.println("[INFO] Combinatorial hits found: "+combi_hits.size());
            for(SynthonSpace.CombinatorialHit chi : combi_hits) {
                System.out.println("[INFO][CombiHit] "+chi.hit_fragments.values().stream().map(li -> ""+li.size()).collect(Collectors.joining("x")));
            }

            if(combi_hits.size()==0) {
                //continue;
                return new ArrayList<>();
            }

            System.out.println("[INFO] Expand combihits, remaining="+conf.getMaxExpandedHits());

            int hits_max = conf.getMaxExpandedHits();
            if(hits_max<0) {hits_max = Integer.MAX_VALUE;}
            int hits_per_hit_max = conf.getMaxExpandedHitsPerCombiHit();
            if(hits_per_hit_max<0) {hits_per_hit_max = Integer.MAX_VALUE;}

            List<HitExpander.SimpleExpandedHit> expandedHits = new ArrayList<>();
            for( SynthonSpace.CombinatorialHit chi : combi_hits) {
                SimpleCombinatorialHit schi = HyperspaceSubstructureSearchForDW.createSimpleCombinatorialHit(chi);
                expandedHits.addAll( HitExpander.expandSimpleHit(schi,Math.min( hits_per_hit_max , hits_max-expandedHits.size())));
                if(hits_max-expandedHits.size()<=0) {
                    System.out.println("[WARN] omitting hits due to limits (maxExpanded)");
                    break;
                }
            }

            System.out.println("[INFO] Computed "+expandedHits.size()+" expanded hits for query: "+mi.getIDCode());
            return expandedHits;

        } catch (AbstractSearchProvider.IncompatibleSearchConfigurationException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // we should not end here..
        System.out.println("[WARN] We should not end here..");
        return new ArrayList<>();
    }

    public static void main(String args[]) {

        String json_file = args[0];
        String json_config = "";
        HyperspaceSubstructureSearchServiceConfig conf_a = null;
        try (FileReader fi_in = new FileReader(json_file)) {
            json_config = IOUtils.toString(new BufferedReader(fi_in));
            conf_a =
                    (new ObjectMapper()).reader().readValue(json_config, HyperspaceSubstructureSearchServiceConfig.class);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        HyperspaceSubstructureSearchService service =
                new HyperspaceSubstructureSearchService();

        service.runService(conf_a);
    }
}
