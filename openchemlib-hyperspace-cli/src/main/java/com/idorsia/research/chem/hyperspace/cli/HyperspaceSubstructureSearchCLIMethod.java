package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.descriptor.DescriptorHandlerFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.hyperspace.SimpleCombinatorialHit;
import com.actelion.research.chem.hyperspace.SimpleSynthon;
import com.actelion.research.chem.io.DWARFileCreator;
import com.actelion.research.chem.io.DWARFileParser;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearchForDW;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSynthonSpaceInitializationProcess;
import com.idorsia.research.chem.hyperspace.util.HitExpander;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.Semaphore;
import java.util.stream.Collectors;

import static com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore.getDefaultInstance;

public class HyperspaceSubstructureSearchCLIMethod implements CLIMethod {



    // Options:

    /**
    public static final String OPT_Space = "space";
    public static final String OPT_DEFAULT_Space = "default.space";
    private int p_Space = -1;

    public static final String OPT_MaxCombi = "max_combi";
    public static final String OPT_DEFAULT_MaxCombi = "32";
    private int p_MaxCombi = -1;

    public static final String OPT_MaxExpanded = "max_expanded";
    public static final String OPT_DEFAULT_MaxExpanded = "4000";
    private int p_MaxExpanded = -1;

    public static final String OPT_MinFFPBits = "min_ffp_bits";
    public static final String OPT_DEFAULT_MinFFPBits = "24";
    private int p_MinFFPBits = -1;

    public static final String OPT_Threads = "threads";
    public static final String OPT_DEFAULT_Threads = "-1";
    private int p_Threads = -1;
     **/

    private static final String P_Space_Msg = "file containing the hyperspace space for search";
    private static final String P_MaxCombi_Msg = "max number of combinatorial hits [default 32]";
    private static final String P_MaxExpanded_Msg = "max number of structures expanded for single query [default 4000]";
    private static final String P_MinFFPBits_Msg = "min number of set ffp bits of query [default 24]";
    private static final String P_Threads_Msg = "number of threads, or -1 to ask system for number of cores [default -1]";
    private static final String P_HSMaxSplits_Msg = "hyperspace sss: max number of splits [default 3]";

    private final CLIParameter P_Space = new CLIParameter("space",true,P_Space_Msg,"default.space",String.class);
    private String p_Space = null;
    private final CLIParameter P_MaxCombi = new CLIParameter("max_combi",true,P_MaxCombi_Msg,"32",Integer.class);
    private int p_MaxCombi = -1;
    private final CLIParameter P_MaxExpanded = new CLIParameter("max_expanded",true,P_MaxExpanded_Msg,"4000",Integer.class);
    private int p_MaxExpanded = -1;
    private final CLIParameter P_MinFFPBits = new CLIParameter("min_ffp",true,P_MinFFPBits_Msg,"24",Integer.class);
    private int p_MinFFPBits = -1;
    private final CLIParameter P_Threads = new CLIParameter("threads",true,P_Threads_Msg,"-1",Integer.class);
    private int p_Threads = -1;
    private final CLIParameter P_HS_MaxSplits = new CLIParameter("hs_maxSplits",true,P_HSMaxSplits_Msg,"3",Integer.class);
    private int p_HS_MaxSplits = -1;

    private List<CLIParameter> parameter;

    public HyperspaceSubstructureSearchCLIMethod() {
        this.parameter = new ArrayList<>();
        parameter.add(P_Space);
        parameter.add(P_MaxCombi);
        parameter.add(P_MaxExpanded);
        parameter.add(P_MinFFPBits);
        parameter.add(P_Threads);
        parameter.add(P_HS_MaxSplits);
    }

    @Override
    public String getMethodName() {
        return "sss";
    }

    @Override
    public List<Option> getCommandLineOptions() {
        List<Option> options = new ArrayList<>();
        options.addAll(this.parameter);
//        options.add(new Option(OPT_Space,true,"file containing the hyperspace space for search"));
//        options.add(new Option(OPT_MaxCombi,true,"max number of combinatorial hits [default 32]"));
//        options.add(new Option(OPT_MaxExpanded,true,"max number of structures expanded for single query [default 4000]"));
//        options.add(new Option(OPT_MinFFPBits,true,"min number of set ffp bits of query [default 24]"));
        return options;
    }

    @Override
    public void parseMethodConfigFromCommandline(CommandLine cli) {

        for(CLIParameter pi : this.parameter) {
            pi.parseConfig(cli);
        }
        this.setP_Space(P_Space.getStringValue());
        this.setP_MaxCombi(P_MaxCombi.getIntValue());
        this.setP_MaxExpanded(P_MaxExpanded.getIntValue());
        this.setP_MinFFPBits(P_MinFFPBits.getIntValue());
        this.setP_Threads(P_Threads.getIntValue());
        this.setP_HS_MaxSplits(P_HS_MaxSplits.getIntValue());

        /**
        String parsed_a = cli.getOptionValue(P_Space, P_Space.getDefaultValue());
        this.setP_Space(Integer.parseInt(parsed_a));
        String parsed_b = cli.getOptionValue(P_MaxCombi, P_MaxCombi.getDefaultValue());
        this.setP_MaxCombi(Integer.parseInt(parsed_b));
        String parsed_c = cli.getOptionValue(P_MaxExpanded, P_MaxExpanded.getDefaultValue());
        this.setP_MaxExpanded(Integer.parseInt(parsed_c));
        String parsed_d = cli.getOptionValue(P_MinFFPBits, P_MinFFPBits.getDefaultValue());
        this.setP_MinFFPBits(Integer.parseInt(parsed_d));
        **/

    }

    /**
     * 1. load space
     * 2. loop over query structures
     *   2.1 check min ffp bits for query
     *   2.2 run search and get combihits for bbs/1/2/3 splits
     *   2.3
     *
     * @param structures
     * @param output
     */
    @Override
    public void runMethod(List<StereoMolecule> structures, String output) {


        // Part 1. configure and run the combi hit search
        HyperspaceSubstructureSearch search = new HyperspaceSubstructureSearch();
        HyperspaceSubstructureSearch.InitializationConfig config_init = new HyperspaceSubstructureSearch.InitializationConfig("[none]",getP_Space(),getP_Threads());
        HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = new HyperspaceSubstructureSearch.SubstructureSearchConfiguration();

        search_config.setMaxSplits(this.getP_HS_MaxSplits());
        search_config.setMaxNumberOfCombinatorialHits(this.getP_MaxCombi());

        search.setConfigurationAndGUI(config_init,null);
        AbstractHyperspaceProcess proc_init = search.startInitialization();

        System.out.println("[INFO] Loading space");

        proc_init.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
            @Override
            public void processStatusChanged() {
                // nothing..
            }
        });

        System.out.println("[INFO] Loading space complete");

        // here we store all hits..
        List<Pair<String, HitExpander.SimpleExpandedHit>> all_collected_hits = new ArrayList<>();

        // Part 2. loop over structures:
        for(StereoMolecule mi : structures) {
            System.out.println("[INFO] Start processing structure: "+mi.getIDCode());

            // count ffp bits of query:
            long[] ffpi  = DescriptorHandlerLongFFP512.getDefaultInstance().createDescriptor(mi) ;
            int ffp_bits = BitSet.valueOf(ffpi).cardinality();
            System.out.println("[INFO] FFP bits in query: " + ffp_bits);
            if(ffp_bits < this.getP_MinFFPBits() ) {
                System.out.println("[INFO] Too few ffp bits -> Skip!");
                continue;
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
                p_search.waitUntilDoneOrFailed(60000);

                if(! p_search.getProcessStatus().equals(AbstractHyperspaceProcess.ProcessStatus.DONE)) {
                    System.out.println("[WARN] bad process status after hyperspace search.. Skip..");
                    continue;
                }

                List<SynthonSpace.CombinatorialHit> combi_hits = p_search.getSearchResults();
                System.out.println("[INFO] Combinatorial hits found: "+combi_hits.size());
                for(SynthonSpace.CombinatorialHit chi : combi_hits) {
                    System.out.println("[INFO][CombiHit] "+chi.hit_fragments.values().stream().map(li -> ""+li.size()).collect(Collectors.joining("x")));
                }

                if(combi_hits.size()==0) {
                    continue;
                }

                System.out.println("[INFO] Expand combihits, remaining="+this.getP_MaxExpanded());

                int hits_max = this.getP_MaxExpanded();
                List<HitExpander.SimpleExpandedHit> expandedHits = new ArrayList<>();
                for( SynthonSpace.CombinatorialHit chi : combi_hits) {
                    SimpleCombinatorialHit schi = HyperspaceSubstructureSearchForDW.createSimpleCombinatorialHit(chi);
                    expandedHits.addAll( HitExpander.expandSimpleHit(schi,hits_max-expandedHits.size()) );
                    if(hits_max-expandedHits.size()<=0) {
                        System.out.println("[WARN] omitting hits due to limits (maxExpanded)");
                        break;
                    }
                }
                // add files to collected hits..
                String query_idc = mi.getIDCode();
                for(HitExpander.SimpleExpandedHit hix : expandedHits) {
                    all_collected_hits.add(Pair.of(query_idc,hix));
                }

                System.out.println("[INFO] Added "+expandedHits.size()+" hits to total hits for query: "+query_idc);

            } catch (AbstractSearchProvider.IncompatibleSearchConfigurationException e) {
                e.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        }

        System.out.println("[INFO] All results computed, write output..");

        // 2. put into dwar:
        HitExpander.exportExpandedHitsWithQueriesToDWAR(output,all_collected_hits);

        System.out.println("[INFO] All results computed, write output.. complete!");

        System.out.println("[INFO] All done, thx for using hyperspace!");
    }




    public int getP_MaxExpanded() {
        return p_MaxExpanded;
    }

    public void setP_MaxExpanded(int p_MaxExpanded) {
        this.p_MaxExpanded = p_MaxExpanded;
    }

    public String getP_Space() {
        return p_Space;
    }

    public void setP_Space(String p_Space) {
        this.p_Space = p_Space;
    }

    public int getP_MaxCombi() {
        return p_MaxCombi;
    }

    public void setP_MaxCombi(int p_MaxCombi) {
        this.p_MaxCombi = p_MaxCombi;
    }

    public int getP_MinFFPBits() {
        return p_MinFFPBits;
    }

    public void setP_MinFFPBits(int p_MinFFPBits) {
        this.p_MinFFPBits = p_MinFFPBits;
    }

    public int getP_Threads() {
        return p_Threads;
    }

    public void setP_Threads(int p_Threads) {
        this.p_Threads = p_Threads;
    }

    public int getP_HS_MaxSplits() {
        return p_HS_MaxSplits;
    }

    public void setP_HS_MaxSplits(int p_HS_MaxSplits) {
        this.p_HS_MaxSplits = p_HS_MaxSplits;
    }
}
