package com.idorsia.research.chem.hyperspace.gui.action;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceInit;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;

public class SaveHyperspaceConfigFileAction extends AbstractAction {

    private HyperspaceSearchGUI gui;

    public SaveHyperspaceConfigFileAction(HyperspaceSearchGUI gui) {
        super("Save config..");
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
            //String filepath = selectedFile.getAbsolutePath();
            HyperspaceInit.saveSearchProviderInitFile(gui,selectedFile);
        } else {
            System.out.println("No file selected.");
        }
    }
}
