package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpaceExplorer3;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.SubstanceFragmentRenderer;
import com.idorsia.research.chem.hyperspace.gui.process.HyperspaceProcessFromProgressController;
import org.apache.commons.lang3.tuple.Pair;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class JSimilarityHitTable extends JPanel {

    enum SimilarityHitTableColType {MOLECULE,STRING,DOUBLE};


    public static class JSimilarityHitTableTableModel extends AbstractTableModel {


        List<String> descriptorNames = new ArrayList<>();
        List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>> hits = new ArrayList<>();

        public JSimilarityHitTableTableModel(List<String> descriptorNames,List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>> hits_i) {
            this.descriptorNames = descriptorNames;
            this.hits = hits_i;
        }


        public void setHits(List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>> hits) {
            this.hits = hits;
        }


        @Override
        public int getRowCount() {
            return this.hits.size();
        }

        @Override
        public int getColumnCount() {
            return 2 + descriptorNames.size() + 6;
        }

        public SimilarityHitTableColType getColumnType(int columnIndex) {
            if(columnIndex==0) {return SimilarityHitTableColType.MOLECULE;}
            if(columnIndex==1) {return SimilarityHitTableColType.STRING;}
            int num_descriptors = this.descriptorNames.size();

            if(columnIndex>1 && (columnIndex-2) < this.descriptorNames.size() ) {
                return SimilarityHitTableColType.DOUBLE;
            }
            else {
                int index = columnIndex - (num_descriptors+2);

                if(index<3) {
                    return SimilarityHitTableColType.MOLECULE;
                }
                else if(index < 6) {
                    return SimilarityHitTableColType.STRING;
                }

            }
            return null;

        }


        /**
         * Column ordering:
         *
         * 0: idcode
         * 1: synthon rxn
         * 2-x : similarities, i.e. all provided double values
         * next three  : frag_idcodes
         * next three  : frag ids
         *
         * @param rowIndex
         * @param columnIndex
         * @return
         */
        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {

            Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]> hit = hits.get(rowIndex);

            if(columnIndex==0) {
                return hit.getLeft().mol.getIDCode();
            }
            if(columnIndex==1) {
                return hit.getLeft().assembly.getSimilarityHit3().getInitialHit().frag.rxn_id;
            }

            int num_descriptorvalues = this.descriptorNames.size();

            if(columnIndex>1 && (columnIndex-2) < this.descriptorNames.size() ) {
                return hit.getRight()[ columnIndex-2 ];
            }
            else {
                int index = columnIndex - (num_descriptorvalues+2);

                if(index<3) {
                    if( hit.getLeft().assembly.getFragmentHits().size()>index ) {
                        return hit.getLeft().assembly.getFragmentHits().get(index).frag.idcode;
                    }
                    else return "";
                }
                else if(index < 6) {
                    if( hit.getLeft().assembly.getFragmentHits().size()>index ) {
                        return hit.getLeft().assembly.getFragmentHits().get(index-3).frag.fragment_id;
                    }
                    else return "NaN";
                }

            }
            return null;
        }

        @Override
        public String getColumnName(int column) {
            if(column==0) {return "Structure";}
            if(column==1) {return "SynthonRxn";}
            if(column>1 && column-2 < this.descriptorNames.size() ) {
                return this.descriptorNames.get(column-2);
            }
            int idx_a = column - (2+descriptorNames.size());
            if(idx_a==0) { return "A"; }
            if(idx_a==1) { return "B"; }
            if(idx_a==2) { return "C"; }
            if(idx_a==3) { return "ID A"; }
            if(idx_a==4) { return "ID B"; }
            if(idx_a==5) { return "ID C"; }
            return "N/A";
        }

        public void exportExpandedToClipboardAsTSV(HyperspaceProcessFromProgressController progress) throws IOException {

            StringWriter str_out = new StringWriter();

            List<String> header_parts = new ArrayList<>();
            header_parts.add("Structure[idcode]");
            header_parts.add("SynthonRxn");
            for(String descr_i : this.descriptorNames) {
                header_parts.add("Sim_"+descr_i);
            }
            header_parts.add("FragA[idcode]");
            header_parts.add("FragB[idcode]");
            header_parts.add("FragC[idcode]");
            header_parts.add("ID_FragA");
            header_parts.add("ID_FragB");
            header_parts.add("ID_FragC");

            String header = String.join("\t",header_parts);

            str_out.write(header+"\n");

            int cnt_lines = 0;

            for(Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>  hit_i : this.hits) {
                List<String> line_parts = new ArrayList<>();
                line_parts.add(hit_i.getLeft().mol.getIDCode());
                line_parts.add(hit_i.getLeft().assembly.getSimilarityHit3().getInitialHit().frag.rxn_id);
                for(int zi=0;zi<this.descriptorNames.size();zi++) {
                    line_parts.add(""+hit_i.getRight()[zi]);
                }
                int num_frags = hit_i.getLeft().assembly.getFragmentHits().size();
                if(num_frags==2) {
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.idcode);
                    line_parts.add("");
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.fragment_id);
                    line_parts.add("NaN");
                }
                else if(num_frags==3) {
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(2).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(2).frag.fragment_id);
                }
                str_out.write(String.join("\t",line_parts)+"\n");
                cnt_lines++;
                if(cnt_lines%100==0) {
                    str_out.flush();
                    progress.startProgress("Copying", cnt_lines, this.hits.size());
                }
            }

            String myString = str_out.getBuffer().toString();
            StringSelection stringSelection = new StringSelection(myString);
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            clipboard.setContents(stringSelection, null);
            progress.setProgress(this.hits.size());
        }

        public void exportExpandedToFileAsCSV(File fileToSave, ProgressController progress) throws IOException {

            BufferedWriter out_file = new BufferedWriter(new FileWriter(fileToSave));
            progress.startProgress("Saving",0,this.hits.size());

            List<String> header_parts = new ArrayList<>();
            header_parts.add("Structure[idcode]");
            header_parts.add("SynthonRxn");
            for(String descr_i : this.descriptorNames) {
                header_parts.add("Sim_"+descr_i);
            }
            header_parts.add("FragA[idcode]");
            header_parts.add("FragB[idcode]");
            header_parts.add("FragC[idcode]");
            header_parts.add("ID_FragA");
            header_parts.add("ID_FragB");
            header_parts.add("ID_FragC");

            String header = String.join(",",header_parts);

            out_file.write(header+"\n");

            int cnt_lines = 0;

            for(Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>  hit_i : this.hits) {
                List<String> line_parts = new ArrayList<>();
                line_parts.add(hit_i.getLeft().mol.getIDCode());
                line_parts.add(hit_i.getLeft().assembly.getSimilarityHit3().getInitialHit().frag.rxn_id);
                for(int zi=0;zi<this.descriptorNames.size();zi++) {
                    line_parts.add(""+hit_i.getRight()[zi]);
                }
                int num_frags = hit_i.getLeft().assembly.getFragmentHits().size();
                if(num_frags==2) {
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.idcode);
                    line_parts.add("");
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.fragment_id);
                    line_parts.add("NaN");
                }
                else if(num_frags==3) {
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(2).frag.idcode);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(0).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(1).frag.fragment_id);
                    line_parts.add(hit_i.getLeft().assembly.getFragmentHits().get(2).frag.fragment_id);
                }
                out_file.write(String.join(",",line_parts)+"\n");
                cnt_lines++;
                if(cnt_lines%100==0) {
                    out_file.flush();
                    progress.startProgress("Saving", cnt_lines, this.hits.size());
                }
            }
            out_file.flush();
            out_file.close();
        }
    }

    JSimilarityHitTableTableModel tablemodel;


    /**
     * GUI
     *
     */
    JPanel jp_top;
    JPanel jp_top_left;
    JPanel jp_top_right;


    JSlider js_size;
    JMenuBar jmb_a;
    JMenu jm_all;

    JScrollPane jsp_table;
    JTable jt_table;









    public JSimilarityHitTable(List<String> descriptorNames) {
        this.tablemodel = new JSimilarityHitTableTableModel(descriptorNames,new ArrayList<>());

        reinitGUI();
        //initPopupMenu();

        // default settings:
        this.js_size.setValue(64);
    }

    public void setSimilarityHits(List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,double[]>> hits) {
        this.tablemodel.setHits(hits);
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

        this.jmb_a  = new JMenuBar();
        this.jp_top_right.add(this.jmb_a);



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




        initToolsMenu();

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


    private void initTableGraphics() {

        for(int zi=0;zi<this.tablemodel.getColumnCount();zi++) {
            if(zi==0) { // assembled structure
                jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                jt_table.getColumnModel().getColumn(zi).setPreferredWidth(HiDPIHelper.scale(400));
            }
            else if( zi>1 && (zi-2)<this.tablemodel.descriptorNames.size()) { // bbs
                //jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                jt_table.getColumnModel().getColumn(zi).setPreferredWidth(HiDPIHelper.scale(80));
            }
            else  {
                if(this.tablemodel.getColumnType(zi)==SimilarityHitTableColType.MOLECULE) {
                    jt_table.getColumnModel().getColumn(zi).setCellRenderer(new SubstanceFragmentRenderer());
                    jt_table.getColumnModel().getColumn(zi).setPreferredWidth(260);
                }
                else {
                    jt_table.getColumnModel().getColumn(zi).setPreferredWidth(100);
                }

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





}
