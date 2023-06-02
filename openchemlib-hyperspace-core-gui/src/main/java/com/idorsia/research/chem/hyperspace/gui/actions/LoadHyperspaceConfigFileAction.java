package com.idorsia.research.chem.hyperspace.gui.actions;

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
            File selectedFile = fileChooser.getSelectedFile();
            System.out.println("Selected File: " + selectedFile.getAbsolutePath());
            //String filepath = selectedFile.getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            String line = null;
            try {
                BufferedReader in = new BufferedReader(new FileReader(selectedFile));
                while ((line = in.readLine()) != null) {
                    sb.append(line);
                    sb.append("\n");
                }
            } catch (FileNotFoundException ex) {
                ex.printStackTrace();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            try {
                Map<String, AbstractSearchProvider> providers = HyperspaceInit.loadSearchProviderInitFile(gui, sb.toString());
                for(String service_name : providers.keySet().stream().sorted().collect(Collectors.toList()) ) {
                    AbstractSearchProvider spi = providers.get(service_name);
                    gui.addSearchProvider(spi);
                }
            } catch (Exception ex) {
                ex.printStackTrace();
            }

        } else {
            System.out.println("No file selected.");
        }
    }
}
