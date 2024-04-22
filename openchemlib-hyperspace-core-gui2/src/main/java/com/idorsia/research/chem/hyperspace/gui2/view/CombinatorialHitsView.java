package com.idorsia.research.chem.hyperspace.gui2.view;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.CombiHitTableModel;
import com.idorsia.research.chem.hyperspace.gui2.model.CombinatorialSearchResultModel;
import com.idorsia.research.chem.hyperspace.gui2.model.RealTimeExpandingSearchResultModel;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CombinatorialHitsView extends JPanel {

    private CombinatorialSearchResultModel resultModel;
    private RealTimeExpandingSearchResultModel expandingModel;



    private Set<String> moleculesInTable = new HashSet<>();

    private JSplitPane  splitPane;
    private JScrollPane scrollPane;
    private JTable      table;
    private RealTimeExpandingHitsView hitsView;

    public CombinatorialHitsView(CombinatorialSearchResultModel resultModel) {
        this.resultModel = resultModel;
        this.reinit();
        this.setExpandedHits(new ArrayList<>());
    }

    public void setExpandedHits(List<SynthonSpace.CombinatorialHit> hits) {
        CombinatorialSearchResultModel combiModel = new CombinatorialSearchResultModel(this.resultModel.getQuery());
        combiModel.addResults(hits);
        this.expandingModel = new RealTimeExpandingSearchResultModel(combiModel,4000);
        this.hitsView = new RealTimeExpandingHitsView(this.expandingModel);
        // Store the current divider location
        int dividerLocation = splitPane.getDividerLocation();
        this.splitPane.setBottomComponent(this.hitsView);
        this.splitPane.setDividerLocation(dividerLocation);
    }

    private void reinit() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        this.scrollPane = new JScrollPane();
        this.table      = new JTable(resultModel.getTableModel());
        this.table.getColumnModel().getColumn(1).setCellRenderer(new FragmentRenderer());
        this.table.getColumnModel().getColumn(2).setCellRenderer(new FragmentRenderer());
        this.table.getColumnModel().getColumn(3).setCellRenderer(new FragmentRenderer());
        this.table.setRowHeight(80);

        this.scrollPane.setViewportView(this.table);
        this.splitPane.setTopComponent(this.scrollPane);
        //this.add(this.scrollPane,BorderLayout.CENTER);
        this.add(this.splitPane,BorderLayout.CENTER);

        this.resultModel.addListener(new CombinatorialSearchResultModel.CombinatorialSearchResultModelListener() {
            @Override
            public void resultsChanged() {

            }
        });
        table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                List<SynthonSpace.CombinatorialHit> selectedHits = new ArrayList<>();
                int[] selectedRows = table.getSelectedRows();
                for( int zi=0;zi<selectedRows.length;zi++) {
                    selectedHits.add( ((CombiHitTableModel) table.getModel()).getHit( selectedRows[zi] ) );
                }
                setExpandedHits(selectedHits);
            }
        });
    }

    private void addAdditionalResults() {

    }

}
