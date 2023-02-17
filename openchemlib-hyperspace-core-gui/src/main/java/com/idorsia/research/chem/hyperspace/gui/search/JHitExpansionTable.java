package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.SubstanceFragmentRenderer;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static com.idorsia.research.chem.hyperspace.SynthonAssembler.expandCombinatorialHit;

public class JHitExpansionTable extends JPanel {

    public static class HitExpansionTableModel  extends AbstractTableModel {

        private int max_synthons = 3;
        private boolean show_bb_ids = false;

        private Map<String,String> fragid_to_bbid = new HashMap<>();


        public HitExpansionTableModel(Map<String,String> fragid_to_bbid, int max_synthons, boolean show_bb_ids) {
            this.fragid_to_bbid = fragid_to_bbid;
            this.max_synthons = max_synthons;
            this.show_bb_ids = show_bb_ids;
        }

        public void setMaxSynthons() {
            this.max_synthons = max_synthons;
        }

        public int getMaxSynthons() {
            return this.max_synthons;
        }

        public void setShowBuildingBlockIDs(boolean show_bb_ids) {
            this.show_bb_ids = show_bb_ids;
        }

        SynthonSpace.CombinatorialHit hit = null;
        List<SynthonAssembler.ExpandedCombinatorialHit> expanded_hit = new ArrayList<>();

        public void setExpandedHit(List<SynthonAssembler.ExpandedCombinatorialHit> expanded_hit) {
            this.expanded_hit = expanded_hit;
            this.fireTableDataChanged();
        }

        public void setHit(SynthonSpace.CombinatorialHit combi_hit) {

            //List<Pair<StereoMolecule, StereoMolecule[]>> expanded = new ArrayList<>();
            Thread compute_expansion = new Thread() {
                @Override
                public void run() {
                    List<SynthonAssembler.ExpandedCombinatorialHit> expanded = expandCombinatorialHit(combi_hit);
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            setExpandedHit(expanded);
                        }
                    });
                }
            };
            compute_expansion.start();
        }

        private String[] alphabeticLetters = new String[]{"A","B","C","D","E","F","G"};

        @Override
        public String getColumnName(int column) {
            if(column==0) {return "SynthonRxn";}
            if(column<max_synthons+1) {
                return alphabeticLetters[column-1];
            }
            return alphabeticLetters[column-1-max_synthons];
        }

        @Override
        public int getRowCount() {
            return expanded_hit.size();
        }

        @Override
        public int getColumnCount() {
            if(show_bb_ids) {
                return 1 + 2 * max_synthons;
            }
            else {
                return 1 + max_synthons;
            }
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SynthonAssembler.ExpandedCombinatorialHit hit = expanded_hit.get(rowIndex);
            if (columnIndex == 0) {
                return hit.assembled_idcode;
            }
            if (columnIndex >= 1 && columnIndex < 1 + max_synthons) {
                int frag_idx = columnIndex - 1;
                if (hit.fragments.size() > frag_idx) {
                    return hit.fragments.get(frag_idx).idcode;
                } else {
                    return "";
                }
            }
            if (columnIndex >= 1 + max_synthons && columnIndex < 1 + 2 * max_synthons) {
                int frag_idx = (columnIndex - 1)-max_synthons;
                if (hit.fragments.size() > frag_idx) {
                    //String bb_id = fragid_to_bbid.get( hit.fragments.get(frag_idx).fragment_id );
                    String bb_id = hit.fragments.get(frag_idx).fragment_id;
                    if(bb_id != null) { return bb_id; }
                    return "<ERROR>";
                } else {
                    return "N/A";
                }
            }
//            if(columnIndex==1){
//                int results = hit.hit_fragments.values().stream().mapToInt( fi -> fi.size() ).reduce( (x,y)->x*y ).getAsInt();
//                return results;
//            }
//            if(columnIndex>=2 && columnIndex<=2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    hit.sri.fragments[columnIndex-2].getIDCode();
//                }
//            }
//            if(columnIndex > 2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    List<SynthonSpace.FragType> fti = new ArrayList<>(hit.hit_fragments.keySet());
//                    fti.sort((x,y)->x.);
//                }
//            }
            return null;
        }

    }


    /**
     * Model
     */
    HitExpansionTableModel tablemodel;


    /**
     * GUI
     *
     */
    JPanel jp_top;
    JPanel jp_top_left;
    JPanel jp_top_right;


    JSlider js_size;

    JScrollPane jsp_table;
    JTable jt_table;


    public JHitExpansionTable(HitExpansionTableModel tablemodel) {
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
        this.jp_top_right.setLayout(new FlowLayout(FlowLayout.LEFT));
        this.jp_top.add(this.jp_top_left,BorderLayout.WEST);
        this.jp_top.add(this.jp_top_right,BorderLayout.EAST);

        this.js_size = new JSlider(8,160);
        this.js_size.setPreferredSize(new Dimension(HiDPIHelper.scale(100),(int) this.js_size.getPreferredSize().getHeight()));

        this.jp_top_right.add(new JLabel("Size "));
        this.jp_top_right.add(this.js_size);


        this.jt_table = new JTable(this.tablemodel);

        // configure table correctly:
        initTableGraphics();

        this.jsp_table = new JScrollPane(this.jt_table);
        this.add(jsp_table,BorderLayout.CENTER);

        this.add(this.jp_top,BorderLayout.NORTH);


        // init listeners:
        this.js_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_table.setRowHeight(js_size.getValue());
            }
        });

    }

    private void initTableGraphics() {



        for(int zi=0;zi<this.tablemodel.getColumnCount();zi++) {
            if(zi==0) { // assembled structure
                jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                jt_table.getColumnModel().getColumn(zi).setPreferredWidth(HiDPIHelper.scale(400));
            }
            else if(zi-1<this.tablemodel.getMaxSynthons()) { // bbs
                jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                jt_table.getColumnModel().getColumn(zi).setPreferredWidth(HiDPIHelper.scale(300));
            }
            else  { // bb ids
                jt_table.getColumnModel().getColumn(zi).setPreferredWidth(100);
            }
        }


//        this.jt_table.getColumnModel().getColumn(2).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(3).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(4).setCellRenderer(new SubstanceFragmentRenderer());

//        this.jt_table.getColumnModel().getColumn(0).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(1).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_table.getColumnModel().getColumn(2).setCellRenderer(new SubstanceFragmentRenderer());
//
//        this.jt_table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
//
//        this.jt_table.setRowHeight(HiDPIHelper.scale(64));

    }

    JPopupMenu jpm_popup;

    private void initPopupMenu() {
        jpm_popup = new JPopupMenu();
        jpm_popup.addPopupMenuListener(new PopupMenuListener() {

            @Override
            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                int rowAtPoint = jt_table.rowAtPoint(SwingUtilities.convertPoint(jpm_popup, new Point(0, 0), jt_table));
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

        JMenuItem jmi_expand = new JMenuItem("Test");
        jmi_expand.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent e) {

                System.out.println("test clicked");
            }
        });
        jpm_popup.add(jmi_expand);
        jt_table.setComponentPopupMenu(jpm_popup);


    }


    private void initListeners() {

        this.js_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_table.setRowHeight(js_size.getValue());
            }
        });

    }



}
