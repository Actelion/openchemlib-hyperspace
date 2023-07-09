package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.gui2.model.RealTimeExpandingSearchResultModel;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.util.HashSet;
import java.util.Set;

public class RealTimeExpandingHitsView extends JPanel {

    private RealTimeExpandingSearchResultModel resultModel;


    private Set<String> moleculesInTable = new HashSet<>();
    private JScrollPane scrollPane;
    private JTable      table;

    private JPanel      panelTop;
    private JLabel      labelResultsInfo;

    private JCheckBox   checkBoxHighlight;
    private JCheckBox   checkBoxAlign;

    public RealTimeExpandingHitsView(RealTimeExpandingSearchResultModel resultModel) {
        this.resultModel = resultModel;
        this.reinit();
    }

    public void setModel(RealTimeExpandingSearchResultModel resultModel) {
        this.resultModel = resultModel;
        this.reinit();
    }

    public RealTimeExpandingSearchResultModel getModel() {
        return this.resultModel;
    }

    private FragmentHighlightChemistryCellRenderer ccr = new FragmentHighlightChemistryCellRenderer(new Dimension(120,120),null);

    private void reinit() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.scrollPane = new JScrollPane();
        this.table      = new JTable(resultModel.getTableModel());

        this.scrollPane.setViewportView(this.table);
        this.add(this.scrollPane,BorderLayout.CENTER);

        this.panelTop = new JPanel();
        this.panelTop.setLayout(new FlowLayout(FlowLayout.LEFT));
        this.labelResultsInfo = new JLabel(resultModel.getResultsInfoString());
        this.panelTop.add(this.labelResultsInfo);
        this.add(panelTop,BorderLayout.NORTH);

        this.checkBoxHighlight = new JCheckBox("Highlight");
        this.checkBoxAlign     = new JCheckBox("Align");

        this.checkBoxHighlight.setSelected(resultModel.isHighlightSubstructure());
        this.checkBoxAlign.setSelected(resultModel.isAlignSubstructure());

        this.panelTop.add(new JLabel("  "));
        this.panelTop.add(this.checkBoxHighlight);
        this.panelTop.add(new JLabel(" "));
        this.panelTop.add(this.checkBoxAlign);

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
    }


}
