package com.idorsia.research.chem.hyperspace.gui.search;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.UtilsGUI;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

public class JSubstructureSearchCombinatorialResultsPanel extends JPanel {

    public static class SubstructureSearchCombinatorialResultsPanelModel {

        private JSubstructureCombinatorialHitTable.JSubstructureCombinatorialHitTableModel combi_hits_table_model;
        private JHitExpansionTable.HitExpansionTableModel                                  expansion_table_model;

        public SubstructureSearchCombinatorialResultsPanelModel(
                JSubstructureCombinatorialHitTable.JSubstructureCombinatorialHitTableModel combi_hits_table_model ,
                JHitExpansionTable.HitExpansionTableModel                                  expansion_table_model) {
            this.combi_hits_table_model = combi_hits_table_model;
            this.expansion_table_model  = expansion_table_model;
        }

        public JSubstructureCombinatorialHitTable.JSubstructureCombinatorialHitTableModel getCombiHitsTableModel() {
            return this.combi_hits_table_model;
        }

        public JHitExpansionTable.HitExpansionTableModel getExpansionTableModel() {
            return this.expansion_table_model;
        }


    }


    JSubstructureSearchCombinatorialResultsPanel.SubstructureSearchCombinatorialResultsPanelModel model;

    public JSubstructureSearchCombinatorialResultsPanel(JSubstructureSearchCombinatorialResultsPanel.SubstructureSearchCombinatorialResultsPanelModel model) {
        this.model = model;
        reinitGUI();
    }

    JSplitPane  jsp_main = null;
    JPanel p_top     = null;
    JPanel p_bottom  = null;

    JSubstructureCombinatorialHitTable jt_combi_hits    = null;
    JHitExpansionTable                 jt_hit_expansion = null;

    public void reinitGUI() {
        this.removeAll();

        this.setLayout(new BorderLayout());

        this.jsp_main = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        this.add(this.jsp_main,BorderLayout.CENTER);

        this.p_top = new JPanel();
        this.p_bottom = new JPanel();
        //this.jsp_main.setLeftComponent(p_left);
        //this.jsp_main.setRightComponent(p_right);
        this.jsp_main.setTopComponent(p_top);
        this.jsp_main.setBottomComponent(p_bottom);

        this.jt_combi_hits = new JSubstructureCombinatorialHitTable(this.model.getCombiHitsTableModel());
        this.p_top.setLayout(new BorderLayout());
        this.p_top.add(this.jt_combi_hits,BorderLayout.CENTER);
        this.p_top.setBorder(new TitledBorder("Combinatorial Hits"));

        this.jt_hit_expansion = new JHitExpansionTable(this.model.getExpansionTableModel());
        this.p_bottom.setLayout(new BorderLayout());
        this.p_bottom.add(this.jt_hit_expansion,BorderLayout.CENTER);
        this.p_bottom.setBorder(new TitledBorder("Enumerated Hits"));

        Runnable setDefaultDivider = new Runnable(){
            @Override
            public void run() {
                jsp_main.setDividerLocation(0.68);
            }
        };
        SwingUtilities.invokeLater(setDefaultDivider);
        // init listeners..

        this.jt_combi_hits.addCombinatorialHitTableListener(new JSubstructureCombinatorialHitTable.CombinatorialHitTableListener() {
            @Override
            public void expandHit(SynthonSpace.CombinatorialHit hit) {
                model.getExpansionTableModel().setHit(hit);
            }
        });

    }

    public void setDefaultDivider() {
        //this.jsp_main.setDividerLocation(0.68);
        UtilsGUI.setDividerLocation(this.jsp_main,0.68);
    }

}
