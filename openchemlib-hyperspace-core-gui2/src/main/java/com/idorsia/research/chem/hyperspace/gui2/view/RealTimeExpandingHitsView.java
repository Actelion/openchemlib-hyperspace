package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.gui2.model.RealTimeExpandingSearchResultModel;
import com.idorsia.research.chem.hyperspace.gui2.task.SubstructureSearchTask;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.HashSet;
import java.util.Set;

public class RealTimeExpandingHitsView extends JPanel {

    private RealTimeExpandingSearchResultModel resultModel;


    private Set<String> moleculesInTable = new HashSet<>();
    private JScrollPane scrollPane;
    private JTable      table;

    private JPanel      panelTop;
    private JPanel      panelTopLeft;
    private JPanel      panelTopRight;
    private JLabel      labelResultsInfo;

    private JCheckBox   checkBoxHighlight;
    private JCheckBox   checkBoxAlign;

    public RealTimeExpandingHitsView(RealTimeExpandingSearchResultModel resultModel) {
        this.resultModel = resultModel;
        this.reinit();
    }

    public void setModel(RealTimeExpandingSearchResultModel resultModel) {
        if(this.resultModel != resultModel) {
            this.resultModel.shutdownThreadpool();
            this.resultModel = resultModel;
            this.resultModel.restartThreadpool();
            this.reinit();
        }
    }

    public RealTimeExpandingSearchResultModel getModel() {
        return this.resultModel;
    }

    private FragmentHighlightChemistryCellRenderer ccr = new FragmentHighlightChemistryCellRenderer(new Dimension(120,120),null);

    protected void reinit() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.scrollPane = new JScrollPane();
        this.table      = new JTable(resultModel.getTableModel());
        this.table.setTableHeader(null);

        this.scrollPane.setViewportView(this.table);
        this.add(this.scrollPane,BorderLayout.CENTER);

        this.panelTop = new JPanel(); this.panelTop.setLayout(new BorderLayout());
        this.add(panelTop,BorderLayout.NORTH);
        this.panelTopLeft  = new JPanel(); this.panelTopLeft.setLayout(new BorderLayout());
        this.panelTopRight = new JPanel(); this.panelTopRight.setLayout(new BorderLayout());
        this.panelTop.add(panelTopLeft,BorderLayout.WEST);
        this.panelTop.add(panelTopRight,BorderLayout.CENTER);

        this.panelTopLeft.setLayout(new FlowLayout(FlowLayout.LEFT));
        this.labelResultsInfo = new JLabel(resultModel.getResultsInfoString());
        this.panelTopLeft.add(this.labelResultsInfo);


        this.checkBoxHighlight = new JCheckBox("Highlight");
        this.checkBoxAlign     = new JCheckBox("Align");

        this.checkBoxHighlight.setSelected(resultModel.isHighlightSubstructure());
        this.checkBoxAlign.setSelected(resultModel.isAlignSubstructure());

        this.panelTopLeft.add(new JLabel("  "));
        this.panelTopLeft.add(this.checkBoxHighlight);
        this.panelTopLeft.add(new JLabel(" "));
        this.panelTopLeft.add(this.checkBoxAlign);

        this.checkBoxHighlight.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                resultModel.setStructurePostprocessOptions(checkBoxHighlight.isSelected(),checkBoxAlign.isSelected());
                if(checkBoxHighlight.isSelected()) {
                    ccr.setHighlightedFragment(resultModel.getCombinatorialSearchResultModel().getQuery());
                }
                else {
                    ccr.setHighlightedFragment(new StereoMolecule());
                }
                SwingUtilities.updateComponentTreeUI(table);
            }
        });
        this.checkBoxAlign.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                resultModel.setStructurePostprocessOptions(checkBoxHighlight.isSelected(),checkBoxAlign.isSelected());
                SwingUtilities.updateComponentTreeUI(table);
            }
        });


        // configure scroll pane:
        //ChemistryCellRenderer ccr = new ChemistryCellRenderer(new Dimension(120,120));

        this.table.setRowHeight(120);
        if(this.resultModel!=null) {
            ccr.setHighlightedFragment(resultModel.getCombinatorialSearchResultModel().getQuery());
            this.table.getColumnModel().getColumn(0).setCellRenderer(ccr);
            this.resultModel.addListener(new RealTimeExpandingSearchResultModel.RealTimeExpandingSearchResultModelListener() {
                @Override
                public void resultsChanged() {
                    // we might arrive here from model-spawned thread, so to make sure..
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            //addAdditionalResults();
                            labelResultsInfo.setText(resultModel.getResultsInfoString());
                            repaint();
                        }
                    });
                }
            });
        }
        this.initMouseContextMenu();
    }

    public JPanel getPanelTopRight() {
        return this.panelTopRight;
    }

    private void initMouseContextMenu() {
        JPopupMenu popupMenu = new JPopupMenu();

        JMenuItem copyStructureSmiles = new JMenuItem("Copy structure as Smiles");
        copyStructureSmiles.addActionListener( (e) -> {
            int selectedRow = table.getSelectedRow();
            if (selectedRow != -1) {
                String ssd = resultModel.getTableModel().getStructureData(selectedRow);
                StereoMolecule si = new StereoMolecule();
                IDCodeParser icp = new IDCodeParser();
                icp.parse(si,ssd);
                Clipboard c = Toolkit.getDefaultToolkit().getSystemClipboard();
                IsomericSmilesCreator isc = new IsomericSmilesCreator(si);
                String smiles_i = isc.getSmiles();
                StringSelection mss = new StringSelection(smiles_i);
                c.setContents(mss,mss);
            }
        } );
        JMenuItem copyStructureIDCode = new JMenuItem("Copy structure as IDCode");
        copyStructureIDCode.addActionListener( (e) -> {
            int selectedRow = table.getSelectedRow();
            if (selectedRow != -1) {
                String ssd = resultModel.getTableModel().getStructureData(selectedRow);
                StringSelection mss = new StringSelection(ssd);
                Clipboard c = Toolkit.getDefaultToolkit().getSystemClipboard();
                c.setContents(mss,mss);
            }
        } );

        popupMenu.add(copyStructureSmiles);
        popupMenu.add(copyStructureIDCode);
        table.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseReleased(MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e) && table.rowAtPoint(e.getPoint()) != -1) {
                    table.setRowSelectionInterval(table.rowAtPoint(e.getPoint()), table.rowAtPoint(e.getPoint()));
                    popupMenu.show(e.getComponent(), e.getX(), e.getY());
                }
            }
        });
    }


}
