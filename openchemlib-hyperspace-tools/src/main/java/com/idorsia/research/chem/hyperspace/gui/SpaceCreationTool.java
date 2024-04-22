package com.idorsia.research.chem.hyperspace.gui;

import com.actelion.research.gui.VerticalFlowLayout;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class SpaceCreationTool {

    public static void main(String args[]) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
             public void run() {
                 FileSelectionDialog dialog = new FileSelectionDialog();
                 dialog.setVisible(true);
             }
        });
    }
    private static class FileSelectionDialog extends JFrame {
        private JTextField textField;
        private JButton selectButton;

        private String selectedFile = null;

        private JTextField textFieldThreads;

        public FileSelectionDialog() {
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
            this.setLayout(new VerticalFlowLayout());
            JPanel pa = new JPanel();
            pa.setLayout(new FlowLayout());
            pa.add(label);
            pa.add(textField);
            pa.add(selectButton);

            JPanel pa2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
            textFieldThreads = new JTextField(""+-1,4);
            JLabel labelThreads = new JLabel("Threads (-1 for #cores) ");
            pa2.add(labelThreads);
            pa2.add(textFieldThreads);


            JPanel pb = new JPanel();
            pb.setLayout(new FlowLayout(FlowLayout.LEFT));
            JButton jb_start = new JButton("Create space");
            pb.add(jb_start);
            this.add(pa);
            this.add(pa2);
            this.add(pb);
            pack();
            setLocationRelativeTo(null); // Center the dialog on the screen


            jb_start.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    int numThreads = -1;
                    try{
                        numThreads = Integer.parseInt(textFieldThreads.getText());
                    }
                    catch(Exception ex) {ex.printStackTrace();}
                    if(numThreads<=0) {
                        numThreads = Runtime.getRuntime().availableProcessors();
                    }
                    System.out.println("Start space creation, Threads = "+numThreads);
                    SynthonSpaceParser2.initREALSpace(selectedFile,getSpaceName(),"FragFp", numThreads);
                }
            });
        }

        public String getSelectedFile() {return this.selectedFile;}
        public String getSpaceName() {return this.textField.getText();}

//        public static void main(String[] args) {
//            SwingUtilities.invokeLater(new Runnable() {
//                @Override
//                public void run() {
//                    FileSelectionDialog dialog = new FileSelectionDialog();
//                    dialog.setVisible(true);
//                }
//            });
//        }
    }

}
