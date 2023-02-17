package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import info.clearthought.layout.TableLayout;
import org.json.JSONArray;
import org.json.JSONObject;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.List;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */


public class HyperspaceSubstructureSearch extends AbstractSearchProvider<SynthonSpace.CombinatorialHit> {


    public static class SubstructureSearchConfiguration implements SearchConfiguration {

        private StereoMolecule query = null;
        private boolean screenBuildingBlocks;
        private int maxSplits;
        private int maxNumberOfCombinatorialHits;

        private int numberOfThreads;

        public SubstructureSearchConfiguration() {}

        public SubstructureSearchConfiguration(StereoMolecule query, boolean screenBuildingBlocks, int max_splits, int max_number_of_combinatorial_hits, int number_of_threads) {
            this.setQuery(query);
            this.screenBuildingBlocks = screenBuildingBlocks;
            this.maxSplits = max_splits;
            this.maxNumberOfCombinatorialHits = max_number_of_combinatorial_hits;
            this.numberOfThreads = number_of_threads;
        }

        @Override
        public List<StereoMolecule> getQueryMolecules() {
            ArrayList<StereoMolecule> query_mols = new ArrayList<>();
            query_mols.add(this.getQuery());
            return query_mols;
        }

        public int getMaxNumberOfCombinatorialHits() {
            return maxNumberOfCombinatorialHits;
        }
        public void setMaxNumberOfCombinatorialHits(int maxNumberOfCombinatorialHits) {
            this.maxNumberOfCombinatorialHits = maxNumberOfCombinatorialHits;
        }
        public int getMaxSplits() {
            return maxSplits;
        }
        public void setMaxSplits(int maxSplits) {
            this.maxSplits = maxSplits;
        }

        public int getNumberOfThreads() {return numberOfThreads;}
        public void setNumberOfThreads(int t) {this.numberOfThreads=t;}

        public StereoMolecule getQuery() {
            return query;
        }

        public void setQuery(StereoMolecule query) {
            this.query = query;
        }

        @Override
        public String serializeToJSON() {
            JSONObject jo = new JSONObject();
            jo.put("screenBuildingBlocks",this.screenBuildingBlocks);
            jo.put("maxSplits",this.maxSplits);
            jo.put("maxNumberOfCombinatorialHits",this.maxNumberOfCombinatorialHits);
            jo.put("query",this.query.getIDCode());
            jo.put("numberOfThreads",this.numberOfThreads);
            return jo.toString();
        }

        @Override
        public void deserializeFromJSON(String json) throws Exception {
            JSONObject jo = new JSONObject(json);
            this.screenBuildingBlocks = jo.getBoolean("screenBuildingBlocks");
            this.maxSplits = jo.getInt("maxSplits");
            this.maxNumberOfCombinatorialHits = jo.getInt("maxNumberOfCombinatorialHits");
            this.query = HyperspaceUtils.parseIDCode(jo.getString("query"));
            this.numberOfThreads = jo.getInt("numberOfThreads");
        }
    }

    public static class InitializationConfig implements SearchProviderConfiguration {
        String spaceName = "";
        String file = null;
        int    maxNumberOfThreads;
        public InitializationConfig() {}
        public InitializationConfig(String space_name, String file, int max_number_of_threads) {
            this.file = file;
            this.spaceName = space_name;
            this.maxNumberOfThreads = max_number_of_threads;
        }
        public String getFile() {
            return this.file;
        }
        public String getSpaceName() {
            return this.spaceName;
        }

        public String serializeToJSON() {
            JSONObject jo = new JSONObject();
            jo.put("SpaceName",this.spaceName);
            jo.put("File",this.file);
            jo.put("maxNumberOfThreads",this.maxNumberOfThreads);
            return jo.toString();
        }

        public void deserializeFromJSON(String json) throws Exception {
            JSONObject jo = new JSONObject(json);
            this.spaceName = jo.getString("SpaceName");
            this.file = jo.getString("File");
            if(jo.has("MaxNumberOfThreads")) {
                this.maxNumberOfThreads = jo.getInt("MaxNumberOfThreads");
            }
            else {
                this.maxNumberOfThreads = -1;
            }
        }
    }

    private InitializationConfig config_initialization;

    private SearchProviderStatus status = SearchProviderStatus.NOT_INITIALIZED;

    private HyperspaceSearchGUI gui;

//    public HyperspaceSubstructureSearch(HyperspaceSearchGUI gui, InitializationConfig config) {
//        this.gui = gui;
//        this.config_initialization = config;
//
//        if(config_initialization==null) {
//            this.setStatus(SearchProviderStatus.REMOTE);
//        }
//
//        if(gui != null) {
//            initGUI();
//        }
//    }

    @Override
    public void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui) {
        this.gui = gui;
        this.config_initialization = (InitializationConfig) config;
        if(config_initialization==null) {
            this.setStatus(SearchProviderStatus.REMOTE);
        }
        if(gui != null) {
            initGUI();
        }
    }

    JPanel gui_search_combined; // will have gui_search_top on top and below gui_search
    JPanel gui_search_top;
    JLabel gui_search_label;
    JButton gui_search_button;
    JToggleButton gui_search_jtb_show_config;
    JPanel gui_search_bottom;
    JPanel gui_search_config;

    /**
     * inits the different guis that are required,for search-config
     */
    private void initGUI() {
        // gui for search button and configuration:

        gui_search_top = new JPanel();
        TableLayout border_layout = new TableLayout(new double[][]{ { TableLayout.FILL , TableLayout.PREFERRED , TableLayout.PREFERRED } , {TableLayout.PREFERRED} });
        border_layout.setHGap(HiDPIHelper.scale(3));
        gui_search_top.setLayout(border_layout);
        if(this.status != SearchProviderStatus.REMOTE) {
            gui_search_label = new JLabel("SSS " + this.config_initialization.getSpaceName());
        }
        else {
            gui_search_label = new JLabel("SSS[Net]"); // the remote search provider is not yet initialized..
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
        gui_search_config.add(new JLabel("Configuration: N/A"));


        gui_search_combined = new JPanel();
        gui_search_combined.setLayout(new BorderLayout());
        gui_search_combined.add(gui_search_top,BorderLayout.NORTH);
        gui_search_combined.add(gui_search_bottom,BorderLayout.CENTER);

        //gui_search_top.add(gui_search_button,BorderLayout.EAST);
        //gui_search_top.setBorder(new LineBorder(Color.red.darker(),1));
        //gui_search_bottom.setBorder(new LineBorder(Color.pink.darker(),1));

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

        gui_search_button.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                System.out.println("action..");
                SubstructureSearchConfiguration ssc = (SubstructureSearchConfiguration) getSearchConfiguration();
                StereoMolecule m_q = new StereoMolecule(gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getDrawPanel().getDrawArea().getMolecule());
                m_q.ensureHelperArrays(Molecule.cHelperCIP);
                ssc.setQuery(m_q);

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
                        AbstractHyperspaceSearchProcess search_process = runSearch(ssc);
                        gui.getProcessListModel().addProcess(search_process);
                    } catch (IncompatibleSearchConfigurationException e1) {
                        e1.printStackTrace();
                    }
                }
            }
        });

    }

    private void setGUISearchLabel(String label) {
        this.gui_search_label.setText(label);
    }

    private RemoteSearchProviderAdapter<SynthonSpace.CombinatorialHit> remoteSearchProvider;

    /**
     * For Remote Execution
     */
    public RemoteSearchProviderAdapter<SynthonSpace.CombinatorialHit> getRemoteSearchProvider() {
        return remoteSearchProvider;
    }

    /**
     * For Remote Execution
     */
    public void setRemoteSearchProvider(RemoteSearchProviderAdapter<SynthonSpace.CombinatorialHit> remoteSearchProvider) {
        this.remoteSearchProvider = remoteSearchProvider;
        this.setGUISearchLabel("SSS[Net] "+this.remoteSearchProvider.getSpaceName());
    }

    @Override
    public boolean hasGUI() {
        return this.gui != null;
    }

    @Override
    public SearchConfiguration getSearchConfiguration() {
        int max_threads = -1;
        if(this.config_initialization!=null) {
            if (this.config_initialization.maxNumberOfThreads > 0) {
                max_threads = this.config_initialization.maxNumberOfThreads;
            }
        }
        SubstructureSearchConfiguration ssc = new SubstructureSearchConfiguration(null,true,3,2000,max_threads);
        return ssc;
    }

    @Override
    public int getMaxNumberOfThreads() {
        return this.config_initialization.maxNumberOfThreads;
    }

    @Override
    public Class getSearchConfigurationClass() {
        return SubstructureSearchConfiguration.class;
    }

    //InitializeHyperspaceSynthonSpaceProcess init_process = null;
    HyperspaceSynthonSpaceInitializationProcess init_process = null;

    SynthonSpace space = null;
    CachedDescriptorProvider cdp = null;

    JPanel p_results;

    @Override
    public AbstractHyperspaceSearchProcess<SynthonSpace.CombinatorialHit> runSearch(SearchConfiguration configuration) throws IncompatibleSearchConfigurationException {
        //p_results = null;//new JPanel();
        this.cdp = new CachedDescriptorProvider("FragFp");
//        if(gui != null) {
//            p_results = new JPanel();
//            gui.getHyperspaceMainPanel().addAnalysisTab("Substructure Search", p_results);
//        }
        HyperspaceSubstructureSearchProcess search_process = new HyperspaceSubstructureSearchProcess(this, this.space, this.cdp);
        search_process.setSearchConfiguration(configuration);
        search_process.startSearchAsync();



        return search_process;
    }

    @Override
    public AbstractHyperspaceProcess startInitialization() {
        //init_process = new InitializeHyperspaceSynthonSpaceProcess(config_initialization.getFile());
        init_process = new HyperspaceSynthonSpaceInitializationProcess(this,config_initialization.getFile());
        init_process.startInitializationAsync();

        init_process.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
            @Override
            public void processStatusChanged() {
                if(init_process.getProcessStatus() == AbstractHyperspaceProcess.ProcessStatus.DONE) {
                    space = init_process.getSpace();
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
    public String getSearchName() {
        return null;
    }

    @Override
    public String getSpaceName() {
        return null;
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
    public String encodeResultsToString(List<SynthonSpace.CombinatorialHit> results) throws IOException {
        ByteArrayOutputStream out_bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream( new BufferedOutputStream( out_bytes ) );
        out.writeObject(results);
        out.flush();
        out_bytes.flush();
        byte[] data_serialization = out_bytes.toByteArray();
        //String b64_data           = Base64.getUrlEncoder().encodeToString(data_serialization);
        String encodedData = Base64.getEncoder().encodeToString(data_serialization);
        String urlEncodedData = URLEncoder.encode(encodedData, StandardCharsets.UTF_8.name());
        return urlEncodedData;
    }

    @Override
    public List<SynthonSpace.CombinatorialHit> decodeResultsFromString(String results_string) throws IOException {
        //byte[] data_serialization = Base64.getUrlDecoder().decode(results_string);
        String decodedData = URLDecoder.decode(results_string, StandardCharsets.UTF_8.name());
        byte[] data_serialization = Base64.getDecoder().decode(decodedData);

        ObjectInputStream in      = new ObjectInputStream( new ByteArrayInputStream(data_serialization));
        List<SynthonSpace.CombinatorialHit> results = null;
        try {
            results = (List<SynthonSpace.CombinatorialHit>)  in.readObject();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.out.println("[decodeResultsFromString] :: Class not found while deserializing..");
            return null;
        }
        return results;
    }

    @Override
    public void presentSearchResults(List<SynthonSpace.CombinatorialHit> results) {

        List<SynthonSpace.CombinatorialHit> all_hits = new ArrayList<>(results);

        Runnable run_gui_update = new Runnable() {
            @Override
            public void run() {
                //JPanel p_tab = new JPanel();
                //p_tab.setLayout(new BorderLayout());
                if(gui != null) {
                    p_results = new JPanel();
                    gui.getHyperspaceMainPanel().addAnalysisTab("Substructure Search", p_results);
                }

                Map<String,String> fragid_to_bb_id = new HashMap<>(); // TODO: implement..
                int max_synthons = 3;

                JSubstructureCombinatorialHitTable.JSubstructureCombinatorialHitTableModel combi_table_model = new JSubstructureCombinatorialHitTable.JSubstructureCombinatorialHitTableModel(fragid_to_bb_id,max_synthons);
                JHitExpansionTable.HitExpansionTableModel expansion_table_model = new JHitExpansionTable.HitExpansionTableModel(fragid_to_bb_id,max_synthons,false);

                JSubstructureSearchCombinatorialResultsPanel.SubstructureSearchCombinatorialResultsPanelModel table_model =
                        new JSubstructureSearchCombinatorialResultsPanel.SubstructureSearchCombinatorialResultsPanelModel(
                                combi_table_model, expansion_table_model);

                table_model.getCombiHitsTableModel().setResults(all_hits);
                JSubstructureSearchCombinatorialResultsPanel p_results_panel = new JSubstructureSearchCombinatorialResultsPanel(table_model);
                //p_tab.add(p_results_panel,BorderLayout.CENTER);

                //gui.getHyperspaceMainPanel().addAnalysisTab("Substructure Search Results",p_tab);
                p_results.setLayout(new BorderLayout());
                p_results.add(p_results_panel,BorderLayout.CENTER);

                SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        p_results_panel.setDefaultDivider();
                    }
                });
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
        return init_process.getStatusOfInitialization();
    }


//    public static class InitializeHyperspaceSynthonSpaceProcess extends AbstractHyperspaceProcess implements AbstractHyperspaceProcess.HasProgress {
//
//        private long file_size_in_bytes;
//        private CountingInputStream counting_in_stream = null;
//
//        private String file;
//        private double ratio_complete = 0.0;
//
//        public InitializeHyperspaceSynthonSpaceProcess(String file) {
//            this.file = file;
//            setProcessStatus(ProcessStatus.WAITING);
//        }
//
//        public void startInitializationAsync() {
//            setProcessStatus(ProcessStatus.COMPUTING);
//            Thread t_init = new Thread() {
//                @Override
//                public void run() {
//
//                    try {
//                        File fi_file = new File(file);
//                        file_size_in_bytes = fi_file.length();
//
//                        FileInputStream file_in = new FileInputStream(fi_file);
//                        counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
//                        ObjectInputStream in = new ObjectInputStream(counting_in_stream);
//                        SynthonSpace space_a = null;
//                        space_a = (SynthonSpace) in.readObject();
//                        space_a.initAfterJavaDeserialization();
//                        space = space_a;
//
//                        setStatus(SearchProviderStatus.READY);
//                        setProcessStatus(ProcessStatus.DONE);
//
//                        return;
//
//                    } catch (FileNotFoundException e) {
//                        e.printStackTrace();
//                    } catch (IOException e) {
//                        e.printStackTrace();
//                    } catch (ClassNotFoundException e) {
//                        e.printStackTrace();
//                    }
//                    setStatus(SearchProviderStatus.ERROR);
//                }
//            };
//            setProcessStatus(ProcessStatus.COMPUTING);
//            t_init.start();
//        }
//
//        public double getStatusOfInitialization() {
//            if(status==SearchProviderStatus.READY) {return 1.0;}
//            if(counting_in_stream==null) {
//                return 0.0;
//            }
//            // all read means 80 percent done:
//            double ratio_done = 0.8 * counting_in_stream.getByteCount() / file_size_in_bytes;
//            return ratio_done;
//        }
//
//        @Override
//        public String getName() {
//            return "Load Synthon Space";
//        }
//
//
//        @Override
//        public List<StereoMolecule> getQueryStructures() {
//            return new ArrayList<>();
//        }
//
//        @Override
//        public double getProgress() {
//            return this.getStatusOfInitialization();
//        }
//    }

}
