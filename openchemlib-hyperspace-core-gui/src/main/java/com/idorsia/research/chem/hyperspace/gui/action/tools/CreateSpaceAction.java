package com.idorsia.research.chem.hyperspace.gui.action.tools;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class CreateSpaceAction extends AbstractAction {

    private HyperspaceSearchGUI gui;

    public CreateSpaceAction(HyperspaceSearchGUI gui) {
        super("Import space..");
        this.gui = gui;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        File file_Input   = null;
        String spaceName  = "spaceA";

        FileSelectionDialog fds = new FileSelectionDialog();
        fds.setVisible(true);

        String spacename = fds.getSpaceName();
        String file_in   = fds.getSelectedFile();


    }


    private class FileSelectionDialog extends JDialog {
        private JTextField textField;
        private JButton selectButton;

        private String selectedFile = null;

        public FileSelectionDialog() {
            super(gui.getFrame(),true);
            setTitle("File Selection Dialog");
            setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            setLayout(new FlowLayout());

            // Create and configure the file chooser
            JFileChooser fileChooser = new JFileChooser();

            // Create the text field and label
            JLabel label = new JLabel("Space Name:");
            textField = new JTextField(20);

            // Create the select button
            selectButton = new JButton("Select File");

            // Add action listener to the select button
            selectButton.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    // Open the file chooser dialog
                    int returnVal = fileChooser.showOpenDialog(FileSelectionDialog.this);
                    if (returnVal == JFileChooser.APPROVE_OPTION) {
                        // Get the selected file
                        String selectedFile_a = fileChooser.getSelectedFile().getAbsolutePath();
                        // Do something with the selected file
                        System.out.println("Selected file: " + selectedFile_a);
                        selectedFile = selectedFile_a;
                    }
                }
            });

            // Add components to the frame
            add(label);
            add(textField);
            add(selectButton);

            pack();
            setLocationRelativeTo(gui.getHyperspaceMainPanel()); // Center the dialog on the screen
        }

        public String getSelectedFile() {return this.selectedFile;}
        public String getSpaceName() {return this.textField.getText();}

    }


}
