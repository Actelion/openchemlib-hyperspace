package com.idorsia.research.chem.hyperspace.gui2.view;

import com.idorsia.research.chem.hyperspace.gui2.model.CombinatorialSearchResultModel;

import javax.swing.*;
import java.awt.*;
import java.util.HashSet;
import java.util.Set;

public class CombinatorialHitsView extends JPanel {

    private CombinatorialSearchResultModel resultModel;


    private Set<String> moleculesInTable = new HashSet<>();
    private JScrollPane scrollPane;
    private JTable      table;

    public CombinatorialHitsView(CombinatorialSearchResultModel resultModel) {
        this.resultModel = resultModel;
        this.reinit();
    }

    private void reinit() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.scrollPane = new JScrollPane();
        this.table      = new JTable(resultModel.getTableModel());

        this.scrollPane.setViewportView(this.table);
        this.add(this.scrollPane,BorderLayout.CENTER);

        this.resultModel.addListener(new CombinatorialSearchResultModel.CombinatorialSearchResultModelListener() {
            @Override
            public void resultsChanged() {

            }
        });
    }

    private void addAdditionalResults() {

    }

}
