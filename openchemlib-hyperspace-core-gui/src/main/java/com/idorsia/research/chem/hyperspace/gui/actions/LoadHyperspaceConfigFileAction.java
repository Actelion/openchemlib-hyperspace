package com.idorsia.research.chem.hyperspace.gui.actions;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;

import javax.swing.*;
import java.awt.event.ActionEvent;

public class LoadHyperspaceConfigFileAction extends AbstractAction {

    private HyperspaceSearchGUI gui;

    public LoadHyperspaceConfigFileAction(HyperspaceSearchGUI gui) {
        super("Load Hyperspace Config File.. [disabled]");
        this.gui = gui;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

    }
}
