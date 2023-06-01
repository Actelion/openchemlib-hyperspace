package com.idorsia.research.chem.hyperspace.gui.actions;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceInit;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;

public class AddSubstructureSearchProviderAction extends AbstractAction {

    private HyperspaceSearchGUI gui;

    public AddSubstructureSearchProviderAction(HyperspaceSearchGUI gui) {
        super("Add Substructure Search Space..");
        this.gui = gui;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select a File");

        int selection = fileChooser.showOpenDialog(null);

        if (selection == JFileChooser.APPROVE_OPTION) {
            File selectedFile = fileChooser.getSelectedFile();
            System.out.println("Selected File: " + selectedFile.getAbsolutePath());
            //HyperspaceSubstructureSearch hss = new HyperspaceSubstructureSearch();
            HyperspaceSubstructureSearch.InitializationConfig conf = new HyperspaceSubstructureSearch.InitializationConfig(selectedFile.getName(),selectedFile.getName(),-1);
            conf.setServiceName("HSSS_"+selectedFile.getName());
            HyperspaceSubstructureSearch hss = new HyperspaceSubstructureSearch();
            hss.setConfigurationAndGUI(conf,gui);
            gui.addSearchProvider(hss);
        } else {
            System.out.println("No file selected.");
        }
    }
}
