package com.idorsia.research.chem.hyperspace.gui2.view;

import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;

import javax.swing.*;
import java.awt.*;

public class SynthonSpaceView extends JPanel {

    private LeetHyperspaceModel model;

    private JScrollPane scrollPane;
    private JTree       tree;

    public SynthonSpaceView(LeetHyperspaceModel model) {
        this.model = model;
        this.reinit();
    }

    private void reinit() {
        this.removeAll();
        this.setLayout(new BorderLayout());
        this.tree = new JTree(model.getSynthonSpacesTreeModel().getTreeModel());
        this.scrollPane = new JScrollPane(this.tree);
        this.add(this.scrollPane,BorderLayout.CENTER);
    }


}
