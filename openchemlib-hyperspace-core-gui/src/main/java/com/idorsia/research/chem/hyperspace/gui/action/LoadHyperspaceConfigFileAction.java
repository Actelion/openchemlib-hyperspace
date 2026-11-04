package com.idorsia.research.chem.hyperspace.gui.action;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceInit;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.Map;
import java.util.stream.Collectors;

public class LoadHyperspaceConfigFileAction extends AbstractAction {

    private HyperspaceSearchGUI gui;

    public LoadHyperspaceConfigFileAction(HyperspaceSearchGUI gui) {
        super("Load config..");
        this.gui = gui;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select a File");

        int selection = fileChooser.showOpenDialog(null);

        if (selection == JFileChooser.APPROVE_OPTION) {
            gui.loadConfiguration(fileChooser.getSelectedFile());

        } else {
            System.out.println("No file selected.");
        }
    }
}
