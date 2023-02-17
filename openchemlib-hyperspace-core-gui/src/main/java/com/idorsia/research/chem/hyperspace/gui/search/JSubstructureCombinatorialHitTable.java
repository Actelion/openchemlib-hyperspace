package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JProgressPanel;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.SubstanceFragmentRenderer;
import com.idorsia.research.chem.hyperspace.gui.process.HyperspaceProcessFromProgressController;
import org.apache.commons.lang3.tuple.Pair;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class JSubstructureCombinatorialHitTable extends JPanel {

    public static class JSubstructureCombinatorialHitTableModel extends AbstractTableModel
    {
        private int max_synthons = 3;

        private Map<String,String> fragid_to_bbid = new HashMap<>();

        public JSubstructureCombinatorialHitTableModel(Map<String,String> fragid_to_bbid, int max_synthons) {
            this.fragid_to_bbid = fragid_to_bbid;
            this.max_synthons = max_synthons;
        }

        List<SynthonSpace.CombinatorialHit> combi_hits = new ArrayList<>();

        public Map<String,String> getFragId2BBId() {
            return this.fragid_to_bbid;
        }

        public SynthonSpace.CombinatorialHit getHit(int idx) {
            return this.combi_hits.get(idx);
        }

        public void setResults(List<SynthonSpace.CombinatorialHit> combi_hits) {
            this.combi_hits = new ArrayList<>(combi_hits);
            fireTableDataChanged();
        }

        private String[] alphabeticLetters = new String[]{"A","B","C","D","E","F","G"};



        public long getNumberOfHitsExpanded() {
            if(this.combi_hits==null) {
                return 0;
            }
            long hits_expanded = this.combi_hits.stream().mapToLong(
                    hi -> hi.hit_fragments.values().stream().mapToLong(vi -> vi.size() ).reduce((x,y) -> x*y ).getAsLong() ).sum();
            return hits_expanded;
        }

        public int getNumberOfCombinatorialHits() {
            if(this.combi_hits==null) {
                return 0;
            }
            return this.combi_hits.size();
        }

        @Override
        public String getColumnName(int column) {
            switch(column) {
                case 0: return "SynthonRxn";
                case 1: return "Hits";
            }
            if(column>=2 && column<2+max_synthons){
                return ""+alphabeticLetters[column-2];
            }
            return alphabeticLetters[column-2-max_synthons];
        }

        @Override
        public int getRowCount() {
            return combi_hits.size();
        }

        @Override
        public int getColumnCount() {
            return 2+2*max_synthons;
        }

        /**
         * @return
         */
        public static List<SynthonSpace.FragType> getSortedFragTypesForHit(SynthonSpace.CombinatorialHit hi) {
            return hi.hit_fragments.keySet().stream().sorted().collect(Collectors.toList());
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SynthonSpace.CombinatorialHit hit = combi_hits.get(rowIndex);
            if(columnIndex==0) {
                return hit.rxn;
            }
            if(columnIndex==1){
                int results = hit.hit_fragments.values().stream().mapToInt( fi -> fi.size() ).reduce( (x,y)->x*y ).getAsInt();
                return results;
            }
            if(columnIndex>=2 && columnIndex<2+max_synthons){
                int frag_idx = columnIndex-2;
                List<SynthonSpace.FragType> sortedFragTypes = getSortedFragTypesForHit(hit);
                if(frag_idx<sortedFragTypes.size()) {
                    SynthonSpace.FragType fti = sortedFragTypes.get(frag_idx);
                    // check if we have this matched:
                    List<Map.Entry<Integer, Pair<SynthonSpace.FragType, BitSet>>> matched_i = hit.mapping.entrySet().stream().filter(xi -> xi.getValue().getLeft().equals(fti)).collect(Collectors.toList());
                    //if( hit.sri.fragments.length > (columnIndex-2) ) {
                    //    return hit.sri.fragments[columnIndex-2].getIDCode();
                    //}
                    if (matched_i.size() > 0) {
                        return hit.sri.fragments[matched_i.get(0).getKey()].getIDCode();
                    }
                    else {
                        return "";
                    }
                }
                return "";
            }
            // TODO: this is actually not correct. we somehow have to add information to the combinatorial hit, which
            //       split result fragment belongs to which frag type.
            if(columnIndex>=2+max_synthons && columnIndex<2+2*max_synthons){
                int frag_idx = (columnIndex-2)-max_synthons;
                List<SynthonSpace.FragType> sortedFragTypes = getSortedFragTypesForHit(hit);
                if(frag_idx<sortedFragTypes.size()) {
                    if (hit.sri.fragments.length > frag_idx) {
                        SynthonSpace.FragType fti = sortedFragTypes.get(frag_idx);
                        // check if we have this matched:
                        return ""+hit.hit_fragments.get(fti).size();
                    }
                }
                return "";
            }
//            if(columnIndex > 2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    List<SynthonSpace.FragType> fti = new ArrayList<>(hit.hit_fragments.keySet());
//                    fti.sort((x,y)->x.);
//                }
//            }

            return null;
        }


        public void exportExpandedToClipboardAsTSV(ProgressController pc) throws IOException {
            if(this.combi_hits==null) {
                System.out.println("[WARN] combi_hits is null, just returning");
                return;
            }

            Runnable ri_export = new Runnable() {
                @Override
                public void run() {
                    StringWriter str_out = new StringWriter();
                    BufferedWriter out = new BufferedWriter(str_out);
                    try {
                        SubstructureSearchHelper.exportHitsToTSV_enumerated(combi_hits,out,new StereoMolecule(),pc);
                        out.flush();
                        str_out.flush();
                        String myString = str_out.getBuffer().toString();
                        StringSelection stringSelection = new StringSelection(myString);
                        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                        clipboard.setContents(stringSelection, null);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            };
            Thread ti = new Thread(ri_export);
            ti.start();
        }

        public void exportExpandedToFileAsCSV(File fileToSave, HyperspaceProcessFromProgressController progress) throws IOException {
            FileWriter file_out = new FileWriter(fileToSave);
            BufferedWriter out = new BufferedWriter(file_out);

            if(this.combi_hits==null) {
                System.out.println("[WARN] combi_hits is null, just returning");
                out.close();
                file_out.close();
                return;
            }
            Runnable ri_export = new Runnable() {
                @Override
                public void run() {
                    try {
                        SubstructureSearchHelper.exportHitsToCSV_enumerated(combi_hits,out,new StereoMolecule(),progress,",");
                        out.flush();
                        file_out.flush();
                        file_out.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            };
            Thread ti = new Thread(ri_export);
            ti.start();
        }
    }


    /**
     * Model
     */
    JSubstructureCombinatorialHitTableModel tablemodel = null;


    /**
     * GUI
     *
     */
    JPanel jp_top;
    JPanel jp_top_left;
    JPanel jp_top_right;

    JLabel jl_hits;

    JSlider js_size;
    JMenuBar jmb_a;
    JMenu    jm_all;

    JScrollPane jsp_table;
    JTable jt_table;


    public JSubstructureCombinatorialHitTable(JSubstructureCombinatorialHitTableModel tablemodel) {
        this.tablemodel = tablemodel;

        reinitGUI();
        initPopupMenu();

        // default settings:
        this.js_size.setValue(64);
    }

    private void reinitGUI() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.jp_top = new JPanel();
        this.jp_top.setLayout(new BorderLayout());

        this.jp_top_left  = new JPanel();
        this.jp_top_right = new JPanel();
        this.jp_top_right.setLayout(new FlowLayout(FlowLayout.RIGHT));
        this.jp_top.add(this.jp_top_left,BorderLayout.WEST);
        this.jp_top.add(this.jp_top_right,BorderLayout.EAST);

        this.js_size = new JSlider(8,160);
        this.js_size.setPreferredSize(new Dimension(HiDPIHelper.scale(100),(int) this.js_size.getPreferredSize().getHeight()));
        this.jp_top_right.add(new JLabel("Size "));
        this.jp_top_right.add(this.js_size);

        this.jmb_a  = new JMenuBar();
        this.jp_top_right.add(this.jmb_a);
        initToolsMenu();

        this.jl_hits = new JLabel("Hits: NA Expanded: NA");
        this.jp_top_left.setLayout(new FlowLayout(FlowLayout.LEFT));
        this.jp_top_left.add(this.jl_hits);
        updateHitsPanel();

        this.jt_table = new JTable(this.tablemodel);

        // configure table correctly:
        initTableGraphics();

        this.jsp_table = new JScrollPane(this.jt_table);
        this.add(jsp_table,BorderLayout.CENTER);

        this.add(jp_top,BorderLayout.NORTH);

        // init listeners:
        this.js_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_table.setRowHeight(js_size.getValue());
            }
        });

        this.tablemodel.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent tableModelEvent) {
                updateHitsPanel();
            }
        });

    }

    private void initToolsMenu() {
        this.jm_all = new JMenu("Tools");
        JMenu export_expanded = new JMenu("Export expanded..");
        JMenuItem export_expanded_clipb_csv  = new JMenuItem("Copy TSV");
        JMenuItem export_expanded_save_dwar = new JMenuItem("Save CSV..");
        export_expanded.add(export_expanded_clipb_csv);
        export_expanded.add(export_expanded_save_dwar);
        this.jm_all.add(export_expanded);
        this.jmb_a.add(this.jm_all);

        export_expanded_clipb_csv.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                HyperspaceProcessFromProgressController progress = new HyperspaceProcessFromProgressController();
                try {
                    tablemodel.exportExpandedToClipboardAsTSV(progress);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                HyperspaceSearchGUI.getGUI().getProcessListModel().addProcess(progress);
            }
        });

        export_expanded_save_dwar.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // parent component of the dialog
                //JFrame parentFrame = (JFrame) SwingUtilities.getRoot(getThis());
                JFrame parentFrame = (JFrame) SwingUtilities.getRoot(HyperspaceSearchGUI.getGUI().getHyperspaceMainPanel());

                JFileChooser fileChooser = new JFileChooser();
                fileChooser.setDialogTitle("Specify a file to save");
                FileNameExtensionFilter csvFilter = new FileNameExtensionFilter("csv files (*.csv)", "csv");
                // add filters
                fileChooser.addChoosableFileFilter(csvFilter);
                fileChooser.setFileFilter(csvFilter);


                int userSelection = fileChooser.showSaveDialog(parentFrame);

                File fileToSave = null;
                if (userSelection == JFileChooser.APPROVE_OPTION) {
                    fileToSave = fileChooser.getSelectedFile();
                    System.out.println("Save as file: " + fileToSave.getAbsolutePath());
                }
                HyperspaceProcessFromProgressController progress = new HyperspaceProcessFromProgressController();
                try {
                    tablemodel.exportExpandedToFileAsCSV(fileToSave,progress);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
                HyperspaceSearchGUI.getGUI().getProcessListModel().addProcess(progress);
            }
        });
    }

    private JComponent getThis() {return this;}

    private void updateHitsPanel() {
        int  hits_combi = this.tablemodel.getNumberOfCombinatorialHits();
        long hits_exp   = this.tablemodel.getNumberOfHitsExpanded();
        this.jl_hits.setText("Hits: "+ ((hits_combi<=0)?"NA":hits_combi) +
                " / "+ ((hits_exp<=0)?"NA":hits_exp));
        this.revalidate();
    }

    private void initTableGraphics() {

        // init renderers and preferred width:
        int preferred_width_rxn_name = (new JLabel("long_reaction_name_a")).getPreferredSize().width;

        int preferred_width_long_number     = (new JLabel("8888888")).getPreferredSize().width;
        int preferred_width_normal_number   = (new JLabel("8888")).getPreferredSize().width;

        for(int zi=0;zi<this.tablemodel.getColumnCount();zi++) {
            if(zi==0) { // rxn name..
                this.jt_table.getColumnModel().getColumn(zi).setPreferredWidth( preferred_width_rxn_name );
                this.jt_table.getColumnModel().getColumn(zi).setMinWidth( preferred_width_rxn_name );
            }
            else {
                if(zi==1) { // total number
                    this.jt_table.getColumnModel().getColumn(zi).setPreferredWidth( preferred_width_long_number );
                    this.jt_table.getColumnModel().getColumn(zi).setMinWidth( preferred_width_long_number );
                }
                else {
                    if(zi-2 < this.tablemodel.max_synthons  ) { // bb structures
                        this.jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                        this.jt_table.getColumnModel().getColumn(zi).setPreferredWidth( 240 );
                    }
                    else { // bb count
                        this.jt_table.getColumnModel().getColumn(zi).setPreferredWidth( preferred_width_normal_number );
                        this.jt_table.getColumnModel().getColumn(zi).setMinWidth( preferred_width_normal_number );
                    }
                }
            }


        }

//        this.jt_table.getColumnModel().getColumn(2).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(3).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(4).setCellRenderer(new SubstanceFragmentRenderer());

        this.jt_table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }

    JPopupMenu jpm_popup;
    int tableMouseOverIdx = -1;

    private void initPopupMenu() {
        // to figure out the actual mouse position we have to add a mouse motion listener as well

        jt_table.addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                tableMouseOverIdx = jt_table.rowAtPoint(e.getPoint());
            }
        });

        jpm_popup = new JPopupMenu();
        jpm_popup.addPopupMenuListener(new PopupMenuListener() {

            @Override
            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                //int rowAtPoint = jt_table.rowAtPoint(SwingUtilities.convertPoint(jpm_popup, new Point(0, 0), jt_table));
                //int rowAtPoint = jt_table.rowAtPoint(SwingUtilities.convertPoint(jpm_popup, new Point(0, 0), jt_table));
                int rowAtPoint = tableMouseOverIdx;
                if (rowAtPoint > -1) {
                    jt_table.setRowSelectionInterval(rowAtPoint, rowAtPoint);
                }

//                SwingUtilities.invokeLater(new Runnable() {
//                    @Override
//                    public void run() {
//                        int rowAtPoint = jt_table.rowAtPoint(SwingUtilities.convertPoint(jpm_popup, new Point(0, 0), jt_table));
//                        if (rowAtPoint > -1) {
//                            jt_table.setRowSelectionInterval(rowAtPoint, rowAtPoint);
//                        }
//                    }
//                });
            }

            @Override
            public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {
                // TODO Auto-generated method stub
            }

            @Override
            public void popupMenuCanceled(PopupMenuEvent e) {
                // TODO Auto-generated method stub
            }
        });

        JMenuItem jmi_expand = new JMenuItem("Expand");
        jmi_expand.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent e) {
                System.out.println("expand row: ");
                fireExpandHitEvent();
            }
        });
        jpm_popup.add(jmi_expand);
        jt_table.setComponentPopupMenu(jpm_popup);

    }

    public static interface CombinatorialHitTableListener {
        public void expandHit(SynthonSpace.CombinatorialHit hit);
    }

    private List<CombinatorialHitTableListener> listeners = new ArrayList<>();

    public void addCombinatorialHitTableListener(CombinatorialHitTableListener li) {
        this.listeners.add(li);
    }

    public void fireExpandHitEvent() {
        int hit_idx = this.jt_table.getSelectedRow();
        for(CombinatorialHitTableListener li : listeners) {
            li.expandHit( this.tablemodel.getHit( hit_idx ) );
        }
    }



    private void initListeners() {

        this.js_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_table.setRowHeight(js_size.getValue());
            }
        });

    }

    public void setHits(List<SynthonSpace.CombinatorialHit> hits) {
        this.tablemodel.setResults(hits);
    }






}
