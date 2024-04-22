package com.idorsia.research.chem.hyperspace.gui2.action;

import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.task.LoadSynthonSpaceTask;
import com.idorsia.research.chem.hyperspace.gui2.view.AbstractLeetHyperspaceView;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.io.File;

public class LoadSynthonSpaceAction extends AbstractLeetHyperspaceAction {

    public LoadSynthonSpaceAction(LeetHyperspaceModel model, AbstractLeetHyperspaceView view) {
        super("Load Synthon Space..", null, model, view);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select a File");
        String userDirectory = System.getProperty("user.dir");
        if(userDirectory!=null) {
            fileChooser.setCurrentDirectory(new File(userDirectory));
        }
        int selection = fileChooser.showOpenDialog(null);

        if (selection == JFileChooser.APPROVE_OPTION) {
            File selectedFile = fileChooser.getSelectedFile();
            System.out.println("Selected File: " + selectedFile.getAbsolutePath());
            //HyperspaceSubstructureSearch hss = new HyperspaceSubstructureSearch();

            LoadSynthonSpaceTask task = new LoadSynthonSpaceTask(getModel(),selectedFile.getAbsolutePath());
            getModel().addTask(task);
            task.execute();
        } else {
            System.out.println("No file selected.");
        }

    }

}
