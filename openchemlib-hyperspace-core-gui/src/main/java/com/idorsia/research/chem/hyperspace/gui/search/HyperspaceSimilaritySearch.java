package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerSkeletonSpheres;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import info.clearthought.layout.TableLayout;
import org.apache.commons.lang3.tuple.Pair;
import org.json.JSONObject;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



public class HyperspaceSimilaritySearch extends AbstractSearchProvider<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> {


    public static class SimilaritySearchConfiguration implements SearchConfiguration {

        private StereoMolecule query = null;
        private int maxSplits;
        private int maxNumberOfFragments;
        private int maxHitsPerLevel;
        private int numberOfThreads;

        public SimilaritySearchConfiguration() {}

        public SimilaritySearchConfiguration(StereoMolecule query, int max_splits, int maxNumberOfFragments, int maxHitsPerLevel, int numberOfThreads) {
            this.setQuery(query);
            this.maxSplits = max_splits;
            this.maxNumberOfFragments = maxNumberOfFragments;
            this.maxHitsPerLevel = maxHitsPerLevel;
            this.numberOfThreads = numberOfThreads;
        }

        @Override
        public List<StereoMolecule> getQueryMolecules() {
            ArrayList<StereoMolecule> query_mols = new ArrayList<>();
            query_mols.add(this.getQuery());
            return query_mols;
        }

        public StereoMolecule getQuery() {
            return query;
        }

        public void setMaxHitsPerLevel(int max_hits) {
            this.maxHitsPerLevel = max_hits;
        }

        public void setQuery(StereoMolecule query) {
            this.query = query;
        }

        @Override
        public int getNumberOfThreads() {
            return this.numberOfThreads;
        }

        public void setNumberOfThreads(int t) {
            this.numberOfThreads = t;
        }

        @Override
        public String serializeToJSON() {
            JSONObject jo = new JSONObject();
            jo.put("maxSplits",this.maxSplits);
            jo.put("maxNumberOfFragments",this.maxNumberOfFragments);
            jo.put("maxHitsPerLevel",this.maxHitsPerLevel);
            jo.put("query",this.query.getIDCode());
            jo.put("numberOfThreads",this.numberOfThreads);
            return jo.toString();
        }

        @Override
        public void deserializeFromJSON(String json) throws Exception {
            JSONObject jo = new JSONObject(json);
            this.maxSplits = jo.getInt("maxSplits");
            this.maxNumberOfFragments = jo.getInt("maxNumberOfFragments");
            this.maxHitsPerLevel = jo.getInt("maxHitsPerLevel");
            this.query = HyperspaceUtils.parseIDCode(jo.getString("query"));
            this.numberOfThreads = jo.getInt("numberOfThreads");
        }
    }


    private HyperspaceSubstructureSearch.InitializationConfig config_initialization;

    private SearchProviderStatus status = SearchProviderStatus.NOT_INITIALIZED;

    private HyperspaceSearchGUI gui;


    private SynthonSimilaritySpace3 space3;

    public HyperspaceSimilaritySearch() {

    }

    @Override
    public void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui) {
        this.gui = gui;
        this.config_initialization = (HyperspaceSubstructureSearch.InitializationConfig) config;
        if(config_initialization==null) {
            this.setStatus(SearchProviderStatus.REMOTE);
        }
        if(gui != null) {
            initGUI();
        }
    }

    public void setSpace(SynthonSimilaritySpace3 space3) {
        this.space3 = space3;
    }

    //JPanel gui_jp_Main = new JPanel();

    JPanel gui_search_combined; // will have gui_search_top on top and below gui_search
    JPanel gui_search_top;
    JLabel gui_search_label;
    JButton gui_search_button;
    JToggleButton gui_search_jtb_show_config;
    JPanel gui_search_bottom;
    JPanel gui_search_config;

    JTextField jtf_hits_per_level;



    private void initGUI() {
        //this.gui_jp_Main = new JPanel();
        //JButton jb_search = new JButton("Search");
        //this.gui_jp_Main.add(jb_search);


        gui_search_top = new JPanel();
        TableLayout border_layout = new TableLayout(new double[][]{ { TableLayout.FILL , TableLayout.PREFERRED , TableLayout.PREFERRED } , {TableLayout.PREFERRED} });
        border_layout.setHGap(HiDPIHelper.scale(3));
        gui_search_top.setLayout(border_layout);
        if(this.status != SearchProviderStatus.REMOTE) {
            gui_search_label = new JLabel("Similarity " + this.config_initialization.getSpaceName());
        }
        else {
            gui_search_label = new JLabel("Similarity[Net]"); // the remote search provider is not yet initialized..
            //gui_search_label = new JLabel("SSS[Net] "+ this.remoteSearchProvider.getSpaceName() );
        }
        gui_search_label.setBorder(new EtchedBorder());
        gui_search_top.add(gui_search_label,"0,0");
        gui_search_button = new JButton("Search");
        gui_search_top.add(gui_search_button,"1,0");
        gui_search_jtb_show_config = new JToggleButton("Config");
        gui_search_top.add(gui_search_jtb_show_config,"2,0");

        gui_search_bottom = new JPanel();
        gui_search_bottom.setLayout(new BorderLayout());

        gui_search_config = new JPanel();
        gui_search_config.setLayout(new FlowLayout());
        gui_search_config.add(new JLabel("Configuration: N/A"));

        if(false) {
            JLabel jl_hits_per_level = new JLabel("hits_per_level ");
            jtf_hits_per_level = new JTextField();
            jtf_hits_per_level.setText("2000");
            gui_search_config.add(jl_hits_per_level);
            gui_search_config.add(jtf_hits_per_level);
        }
        //gui_search_bottom.add(gui_search_config,BorderLayout.CENTER);

        gui_search_combined = new JPanel();
        gui_search_combined.setLayout(new BorderLayout());
        gui_search_combined.add(gui_search_top,BorderLayout.NORTH);
        gui_search_combined.add(gui_search_bottom,BorderLayout.CENTER);

        gui_search_jtb_show_config.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(gui_search_jtb_show_config.isSelected()) {
                    gui_search_bottom.remove(gui_search_config);
                    gui_search_bottom.add(gui_search_config,BorderLayout.CENTER);
                }
                else {
                    gui_search_bottom.remove(gui_search_config);
                }
                gui_search_combined.revalidate();
            }
        });

        gui_search_button.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent actionEvent) {

                StereoMolecule query_mol = new StereoMolecule(gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getDrawPanel().getDrawArea().getMolecule());
                query_mol.ensureHelperArrays(Molecule.cHelperCIP);

                SimilaritySearchConfiguration ssc = new SimilaritySearchConfiguration(query_mol, 3, 3, 2000, -1);
                //SimilaritySearchConfiguration ssc = (SimilaritySearchConfiguration) getSearchConfiguration();

                if(getStatus() == SearchProviderStatus.REMOTE) {
                    try {
                        AbstractHyperspaceSearchProcess proc = getRemoteSearchProvider().runSearch(ssc);
                        gui.getProcessListModel().addProcess(proc);
                    } catch (IncompatibleSearchConfigurationException incompatibleSearchConfigurationException) {
                        incompatibleSearchConfigurationException.printStackTrace();
                    }
                }
                else {
                    try {
                        AbstractHyperspaceSearchProcess proc = runSearch(ssc);
                        gui.getProcessListModel().addProcess(proc);
                    } catch (IncompatibleSearchConfigurationException e) {
                        e.printStackTrace();
                    }
                }
            }
        });
    }

    @Override
    public String getSearchName() {
        return "Hyperspace Similarity Search 3";
    }

    @Override
    public String getSpaceName() {
        return config_initialization.getSpaceName();
    }

    @Override
    public JPanel getSearchConfigurationView() {
        if(this.status == SearchProviderStatus.NOT_INITIALIZED) {
            JPanel pa = new JPanel(); pa.setLayout(new FlowLayout(FlowLayout.LEFT));
            pa.add(new JLabel( "Not initialized: "+config_initialization.getSpaceName() ));
            return pa;
        }
        if(this.status == SearchProviderStatus.INITIALIZING) {
            JPanel pa = new JPanel(); pa.setLayout(new FlowLayout(FlowLayout.LEFT));
            pa.add(new JLabel( "Initializting: "+config_initialization.getSpaceName() ));
            return pa;
        }
        if(this.status == SearchProviderStatus.READY || this.status == SearchProviderStatus.REMOTE ) {
//            JPanel pa = new JPanel(); pa.setLayout(new FlowLayout(FlowLayout.LEFT));
//            pa.add(new JLabel( "SSS: "+config_initialization.getSpaceName()));
//            JButton bs = new JButton("Search");
//            pa.add(bs);
            JPanel pa = gui_search_combined;
            return pa;
        }

        JPanel pi = new JPanel(); pi.add(new JLabel("?"));
        return pi;
    }

    @Override
    public SearchConfiguration getSearchConfiguration() {

        int num_threads = -1;
        if(this.config_initialization.maxNumberOfThreads>0) {
            num_threads = this.config_initialization.maxNumberOfThreads;
        }

        if(gui != null) {

            int hits_per_level = 2000;
            try{
                hits_per_level = Integer.parseInt(this.jtf_hits_per_level.getText());
            }
            catch(Exception ex) {
                System.out.println("[WARN] Could not parse hits_per_level");
                //JOptionPane.showMessageDialog( SwingUtilities.getWindowAncestor(this.gui.getHyperspaceMainPanel()) ,"Could not parse hits_per_level, take standard value "+hits_per_level);
            }

            StereoMolecule mi = new StereoMolecule(this.gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getDrawPanel().getDrawArea().getMolecule());
            mi.ensureHelperArrays(Molecule.cHelperCIP);
            SimilaritySearchConfiguration s_config = new SimilaritySearchConfiguration(mi,3,3,hits_per_level,num_threads);

            return s_config;
        }

        return new SimilaritySearchConfiguration();
    }

    @Override
    public int getMaxNumberOfThreads() {
        return this.config_initialization.maxNumberOfThreads;
    }

    @Override
    public Class getSearchConfigurationClass() {
        return SimilaritySearchConfiguration.class;
    }

    HyperspaceSimilaritySearchInitializationProcess init_process = null;

    @Override
    public AbstractHyperspaceProcess startInitialization() {
        //init_process = new InitializeHyperspaceSynthonSpaceProcess(config_initialization.getFile());
        this.init_process = new HyperspaceSimilaritySearchInitializationProcess(this,config_initialization.getFile());
        init_process.startInitializationAsync();

        init_process.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
            @Override
            public void processStatusChanged() {
                if(init_process.getProcessStatus() == AbstractHyperspaceProcess.ProcessStatus.DONE) {
                    space3 = (SynthonSimilaritySpace3) init_process.getSpace();
                    setStatus(SearchProviderStatus.READY);
                }
                else if(init_process.getProcessStatus() == AbstractHyperspaceProcess.ProcessStatus.FAILED ) {
                    setStatus(SearchProviderStatus.ERROR);
                }
            }
        });
        return init_process;
    }

    @Override
    public AbstractHyperspaceSearchProcess<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> runSearch(SearchConfiguration configuration) throws IncompatibleSearchConfigurationException {
        HyperspaceSimilaritySearchProcess process = new HyperspaceSimilaritySearchProcess(this,this.space3);
        process.setSearchConfiguration( (SimilaritySearchConfiguration) configuration);
        process.startSearchAsync();
        return process;
    }


    // Results GUI:
    JPanel p_results = new JPanel();

    @Override
    public void presentSearchResults(List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> results) {

        //List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> all_hits = new ArrayList<>(results);


        Runnable run_gui_update = new Runnable() {
            @Override
            public void run() {
                //JPanel p_tab = new JPanel();
                //p_tab.setLayout(new BorderLayout());
                if(gui != null) {
                    p_results = new JPanel();
                    gui.getHyperspaceMainPanel().addAnalysisTab("Similarity Search", p_results);
                }

                Map<String,String> fragid_to_bb_id = new HashMap<>(); // TODO: implement..
                int max_synthons = 3;


                List<String> descriptor_names = new ArrayList<>();
                descriptor_names.add("SkelSpheres");
                descriptor_names.add("FragFp");
                descriptor_names.add("AssemblyScore");
                JSimilarityHitTable table = new JSimilarityHitTable(descriptor_names);

                // create pairs containing the scores:
                DescriptorHandlerSkeletonSpheres dh_sk  = new DescriptorHandlerSkeletonSpheres();
                DescriptorHandlerLongFFP512      dh_ffp = new DescriptorHandlerLongFFP512();
                StereoMolecule m_q = gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getDrawPanel().getDrawArea().getMolecule();
                BitSet         fp_q_ffp  = BitSet.valueOf((dh_ffp).createDescriptor(m_q));
                byte[]         fp_q_skel = (dh_sk).createDescriptor(m_q);

                List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>> hits_scores = new ArrayList<>();
                for(SynthonSimilaritySpaceExplorer3.FinalSimilarityResult fri : results) {
                    double sim_sk  = dh_sk.getSimilarity(fp_q_skel,fri.getSkelSpheres());
                    double sim_ffp = LSHProvider.tanimoto_similarity(fp_q_ffp,fri.getFFP());
                    hits_scores.add(Pair.of(fri,new double[]{sim_sk,sim_ffp, fri.assembly_score }));
                }

                table.setSimilarityHits(hits_scores);
                //JSubstructureSearchCombinatorialResultsPanel p_results_panel = new JSubstructureSearchCombinatorialResultsPanel(table_model);
                //p_tab.add(p_results_panel,BorderLayout.CENTER);

                //gui.getHyperspaceMainPanel().addAnalysisTab("Substructure Search Results",p_tab);
                p_results.setLayout(new BorderLayout());
                p_results.add(table,BorderLayout.CENTER);

//                SwingUtilities.invokeLater(new Runnable() {
//                    @Override
//                    public void run() {
//                        p_results.setDefaultDivider();
//                    }
//                });
            }
        };

        SwingUtilities.invokeLater(run_gui_update);

    }

    @Override
    public SearchProviderStatus getStatus() {
        return status;
    }

    private void setStatus(SearchProviderStatus status) {
        this.status = status;
        fireSearchStatusChanged();
    }

    @Override
    public double getStatusOfInitialization() {
        if(this.init_process==null) {return 0;}
        return this.init_process.getStatusOfInitialization();
    }

    private void setGUISearchLabel(String label) {
        this.gui_search_label.setText(label);
    }

    private RemoteSearchProviderAdapter<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> remoteSearchProvider;

    /**
     * For Remote Execution
     */
    public RemoteSearchProviderAdapter<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> getRemoteSearchProvider() {
        return remoteSearchProvider;
    }

    /**
     * For Remote Execution
     */
    public void setRemoteSearchProvider(RemoteSearchProviderAdapter<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> remoteSearchProvider) {
        this.remoteSearchProvider = remoteSearchProvider;
        this.setGUISearchLabel("Similarity[Net] "+this.remoteSearchProvider.getSpaceName());
    }




    @Override
    public String encodeResultsToString(List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> results) throws IOException {
        ByteArrayOutputStream out_bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream( new BufferedOutputStream( out_bytes ) );
        out.writeObject(results);
        out.flush();
        out_bytes.flush();
        byte[] data_serialization = out_bytes.toByteArray();
        String b64_data           = Base64.getEncoder().encodeToString(data_serialization);
        return b64_data;
    }

    @Override
    public List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> decodeResultsFromString(String results_string) throws IOException {
        byte[] data_serialization = Base64.getMimeDecoder().decode(results_string);
        ObjectInputStream in      = new ObjectInputStream( new ByteArrayInputStream(data_serialization));
        List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> results = null;
        try {
            results = (List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult>)  in.readObject();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.out.println("[decodeResultsFromString] :: Class not found while deserializing..");
            return null;
        }
        return results;
    }

    @Override
    public boolean hasGUI() {
        return this.gui != null;
    }
}
