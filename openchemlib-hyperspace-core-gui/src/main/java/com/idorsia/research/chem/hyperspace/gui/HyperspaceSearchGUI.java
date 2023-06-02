package com.idorsia.research.chem.hyperspace.gui;

import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
//import com.idorsia.research.chem.hyperspace.gui.action.AddSubstructureSearchProviderAction;
//import com.idorsia.research.chem.hyperspace.gui.action.LoadHyperspaceConfigFileAction;
//import com.idorsia.research.chem.hyperspace.gui.action.SaveHyperspaceConfigFileAction;
import com.idorsia.research.chem.hyperspace.gui.action.AddSubstructureSearchProviderAction;
import com.idorsia.research.chem.hyperspace.gui.action.LoadHyperspaceConfigFileAction;
import com.idorsia.research.chem.hyperspace.gui.action.SaveHyperspaceConfigFileAction;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.JProcessList;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import org.pushingpixels.substance.api.SubstanceCortex;
import org.pushingpixels.substance.api.SubstanceSlices;
import org.pushingpixels.substance.api.skin.SubstanceBusinessLookAndFeel;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.stream.Collectors;


/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



public class HyperspaceSearchGUI {

    public static class JMainPanel extends JPanel {

        // MODEL

        JProcessList.ProcessListModel process_list_model;



        // LAYOUT STUFF

        JSplitPane jsp_Main; // splits all into left A and right B

        JSplitPane jsp_a; // splits A into top and bottom
        //JSplitPane jsp_a_a; // splits top into left and right

        JPanel A_top;
        JPanel A_bottom;

        // this is the compute right hand side of the main split pane
        JPanel B_main;

        JTabbedPane B_main_tabbed_pane;

        // LAYOUT STUFF END


        // AREA A Stuff
        JHyperspaceSearchPanel jh_search;

        JProcessList jpl_processes;

        // AREA A Stuff END



        public JMainPanel() {
            // init model:
            this.process_list_model = new JProcessList.ProcessListModel();


            initGUI();
        }


        private void initGUI() {
            //initIconAndTitle();

            this.removeAll();
            this.setLayout(new BorderLayout());

            // init a/b split pane:
            jsp_Main = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
            this.add(jsp_Main,BorderLayout.CENTER);

            // init split panes:
            jsp_a      = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

            B_main = new JPanel(); B_main.setLayout(new BorderLayout());

            A_top      = new JPanel(); A_top.setLayout(new BorderLayout());
            A_bottom   = new JPanel(); A_bottom.setLayout(new BorderLayout());

            jsp_a.setTopComponent(A_top);
            jsp_a.setBottomComponent(A_bottom);

            // init JHyperspaceSearchPanel
            jh_search = new JHyperspaceSearchPanel();
            A_top.add(jh_search);

            jsp_Main.setLeftComponent(jsp_a);
            jsp_Main.setRightComponent(B_main);

            // init B:
            B_main.setBorder(BorderFactory.createTitledBorder("Results & Analysis"));
            B_main_tabbed_pane = new JTabbedPane();
            B_main.add(B_main_tabbed_pane,BorderLayout.CENTER);

            this.jpl_processes = new JProcessList(this.process_list_model);
            this.A_bottom.add(this.jpl_processes,BorderLayout.CENTER);

        }

        /**
         * Add a tab into the results and analysis tabbed pane
         *
         * @param tab_title
         * @param analysis
         * @return
         */
        public void addAnalysisTab(String tab_title, JPanel analysis) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    Set<String> existing_titles = new HashSet<>();
                    for(int zi=0;zi<B_main_tabbed_pane.getTabCount();zi++) {
                        B_main_tabbed_pane.getTitleAt(zi);
                    }

                    String title_i = tab_title;
                    if(existing_titles.contains(title_i)) {
                        boolean title_ok = false;
                        int postfix = 1;
                        while(!title_ok && postfix < 200) {
                            title_i = tab_title+"_"+postfix;
                            title_ok = !existing_titles.contains(title_i);
                            postfix++;
                        }
                    }

                    B_main_tabbed_pane.removeAll();
                    B_main_tabbed_pane.addTab(title_i,analysis);
                }
            });
        }

        public JHyperspaceSearchPanel getHyperspaceSearchPanel() {
            return this.jh_search;
        }

        public void restoreDefaultLayout() {
            this.jsp_a.setDividerLocation(0.7);
            this.jsp_Main.setDividerLocation(0.7);
            this.getHyperspaceSearchPanel().restoreLayout();
        }

        public JProcessList.ProcessListModel getProcessListModel() {
            return this.process_list_model;
        }

    }


    public HyperspaceSearchGUI() {
        //init();
    }

    JFrame f_main;
    JMainPanel p_main;


    public void addSearchProvider(AbstractSearchProvider searchProvider) {
        AbstractHyperspaceProcess init_p = searchProvider.startInitialization();
        this.p_main.process_list_model.addProcess(init_p);
        this.p_main.getHyperspaceSearchPanel().getSearchProviderListPanel().addSearchProvider(searchProvider);
    }


    private void initIconAndTitle() {
        f_main.setTitle("Hyperspace 2.0.6");
        java.net.URL url = this.getClass().getClassLoader().getSystemResource("com/idorsia/research/hyperspace_logo_02.png");
        Image img = Toolkit.getDefaultToolkit().createImage(url);
        this.f_main.setIconImage(img);
        //this.setTitle("Hyperspace Client 1.1.2");
    }


    private void init(String config_file) {

        f_main = new JFrame();
        f_main.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        f_main.getContentPane().setLayout(new BorderLayout());

        initIconAndTitle();

        GraphicsDevice gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
        int screen_width = gd.getDisplayMode().getWidth();
        int screen_height = gd.getDisplayMode().getHeight();

        f_main.setSize(Math.min(1200,screen_width),Math.min(800,screen_height));
        f_main.setVisible(true);

        initMenu();

        p_main = new JMainPanel();
        f_main.getContentPane().add(p_main,BorderLayout.CENTER);

        // init:
        String json_config = HyperspaceUtils.readFileIntoString(config_file);
        if(json_config==null) {
            File fi_a = new File(config_file);
            if(!fi_a.exists()) {
                System.out.println("[INFO] Try to create config file..");
                File fi = new File(config_file);
                json_config = config_file;
            }
        }

        //JSONObject joc = new JSONObject(json_config);
        try {
            Map<String,AbstractSearchProvider> providers = HyperspaceInit.loadSearchProviderInitFile(this,json_config);
            for(String service_name : providers.keySet().stream().sorted().collect(Collectors.toList()) ) {
                AbstractSearchProvider spi = providers.get(service_name);
                this.addSearchProvider(spi);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }


        if(false) {
            // perform some initialization:
            //HyperspaceSubstructureSearch.InitializationConfig conf_a = new HyperspaceSubstructureSearch.InitializationConfig("REAL space","C:\\dev4\\git4\\hyperspace_core\\real_space_latest_FragFp.data");
            //HyperspaceSubstructureSearch.InitializationConfig conf_b = new HyperspaceSubstructureSearch.InitializationConfig("REAL space","C:\\dev4\\git4\\hyperspace_core\\real_space_small_latest_FragFp.data");
            HyperspaceSubstructureSearch.InitializationConfig conf_b = new HyperspaceSubstructureSearch.InitializationConfig("REAL space","/home/liphath1/dev5/git/hyperspace/REAL_Space_latest_FragFp.data",-1);




            String hs_base_linux = "/home/liphath1/hyperspace_base_2/data_hyperspace";

//        HyperspaceSubstructureSearch.InitializationConfig conf_ivcs_2s = new HyperspaceSubstructureSearch.InitializationConfig("IVCS 3.2 (2s)"      ,"C:\\dev4\\git4\\hyperspace_server\\ivcs_3_2_FULL_only_2s_FragFp.data");
//        HyperspaceSubstructureSearch.InitializationConfig conf_ivcs_3s = new HyperspaceSubstructureSearch.InitializationConfig("IVCS 3.2 (3s) Mini" ,"C:\\dev4\\git4\\hyperspace_server\\ivcs_3_2_FULL_only_2s_FragFp.data");

            HyperspaceSubstructureSearch.InitializationConfig conf_ivcs_2s = new HyperspaceSubstructureSearch.InitializationConfig("IVCS 3.2 (2s)", "/home/liphath1/dev5/git/hyperspace/ivcs3_2s_FragFp.data",-1);
            HyperspaceSubstructureSearch.InitializationConfig conf_ivcs_3s_a = new HyperspaceSubstructureSearch.InitializationConfig("IVCS 3.2 (3s) Mini", "/home/liphath1/dev5/git/hyperspace/ivcs3_3s_FragFp.data",-1);
            HyperspaceSubstructureSearch.InitializationConfig conf_ivcs_3s_b = new HyperspaceSubstructureSearch.InitializationConfig("IVCS 3.2 (3s) Full", "/home/liphath1/hyperspace_base_2/data_hyperspace/ivcs3_3s_FragFp.data",-1);
            //HyperspaceSubstructureSearch sss_a = new HyperspaceSubstructureSearch(conf_a);

//            HyperspaceSubstructureSearch sss_b = new HyperspaceSubstructureSearch(this, conf_b);
//            HyperspaceSubstructureSearch sss_ivcs_2s = new HyperspaceSubstructureSearch(this, conf_ivcs_2s);
//            HyperspaceSubstructureSearch sss_ivcs_3s_a = new HyperspaceSubstructureSearch(this, conf_ivcs_3s_a);

            HyperspaceSubstructureSearch sss_b = new HyperspaceSubstructureSearch();
            HyperspaceSubstructureSearch sss_ivcs_2s = new HyperspaceSubstructureSearch();
            HyperspaceSubstructureSearch sss_ivcs_3s_a = new HyperspaceSubstructureSearch();

            sss_b.setConfigurationAndGUI(conf_b,this);
            sss_ivcs_2s.setConfigurationAndGUI(conf_ivcs_2s,this);
            sss_ivcs_3s_a.setConfigurationAndGUI(conf_ivcs_3s_a,this);

//        HyperspaceSubstructureSearch sss_ivcs_3s_b = new HyperspaceSubstructureSearch(this,conf_ivcs_3s_b);

            this.addSearchProvider(sss_b);
            this.addSearchProvider(sss_ivcs_2s);
            //this.addSearchProvider(sss_ivcs_3s_a);

            //AbstractHyperspaceProcess process_a = sss_a.startInitialization();

            //AbstractHyperspaceProcess process_b       = sss_b.startInitialization();
            //AbstractHyperspaceProcess process_ivcs_2s = sss_ivcs_2s.startInitialization();
            //AbstractHyperspaceProcess process_ivcs_3s = sss_ivcs_3s_a.startInitialization();
            //AbstractHyperspaceProcess process_ivcs_3s = sss_ivcs_3s_b.startInitialization();

            //main_panel.process_list_model.addProcess(process_a);

            //p_main.process_list_model.addProcess(process_b);
            //p_main.process_list_model.addProcess(process_ivcs_2s);
            //p_main.process_list_model.addProcess(process_ivcs_3s);
            //main_panel.getHyperspaceSearchPanel().getSearchProviderListPanel().addSearchProvider(sss_a);

            //p_main.getHyperspaceSearchPanel().getSearchProviderListPanel().addSearchProvider(sss_b);
            //p_main.getHyperspaceSearchPanel().getSearchProviderListPanel().addSearchProvider(sss_ivcs_2s);
            //p_main.getHyperspaceSearchPanel().getSearchProviderListPanel().addSearchProvider(sss_ivcs_3s_b);
        }

        SwingUtilities.invokeLater(new Runnable() {
                                       @Override
                                       public void run() {
                                           p_main.restoreDefaultLayout();
                                       }
                                   }
        );


        // test..
        //main

    }

    public JMainPanel getHyperspaceMainPanel() {
        return this.p_main;
    }

    public JProcessList.ProcessListModel getProcessListModel() {
        return this.p_main.getProcessListModel();
    }


    private static HyperspaceSearchGUI gui;

    public static void main(String args[]) {

        System.out.println("Hyperspace Search 2.0");
        System.out.println("Idorsia SciComp-Chemistry 2021");
        System.out.println();

        String config_file = "hyperspace_config.json";
        if(args.length==0) {
        }
        else {
            config_file = args[0];
        }
        String final_config_file = config_file;
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {

                try {
                    UIManager.setLookAndFeel(new SubstanceBusinessLookAndFeel());
                    SubstanceCortex.GlobalScope.disallowAnimations(SubstanceSlices.AnimationFacet.ROLLOVER);
                    SubstanceCortex.GlobalScope.disallowAnimations(SubstanceSlices.AnimationFacet.SELECTION);

                } catch (UnsupportedLookAndFeelException e) {
                    e.printStackTrace();
                }

                gui = new HyperspaceSearchGUI();
                gui.init(final_config_file);
            }
        });
    }

    public void initMenu() {
        JMenuBar menubar = new JMenuBar();
        JMenu mfile = new JMenu("File");
        mfile.add(new LoadHyperspaceConfigFileAction(gui));
        mfile.add(new SaveHyperspaceConfigFileAction(gui));
        mfile.addSeparator();
        mfile.add(new AddSubstructureSearchProviderAction(gui));
        menubar.add(mfile);
        this.f_main.setJMenuBar(menubar);
        this.f_main.validate();
    }

    public static HyperspaceSearchGUI getGUI() {
        return gui;
    }

}
