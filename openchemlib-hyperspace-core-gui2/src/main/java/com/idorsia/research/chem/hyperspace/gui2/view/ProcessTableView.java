package com.idorsia.research.chem.hyperspace.gui2.view;

import com.idorsia.research.chem.hyperspace.gui2.model.ProcessTableModel;
import com.idorsia.research.chem.hyperspace.gui2.task.LongRunningTask;

import javax.swing.*;
import java.awt.*;

public class ProcessTableView extends JPanel {

    private ProcessTableModel model;

    private JTable table;

    public ProcessTableView(ProcessTableModel model) {
        this.model = model;
        reinit();
    }

    private void reinit() {
        this.removeAll();
        setLayout(new BorderLayout());

        // Create the table with three columns: Name, Status, Progress
        String[] columnNames = {"Name", "Status", "Progress"};
        table = new JTable(model.getTableModel());
        JScrollPane scrollPane = new JScrollPane(table);
        add(scrollPane, BorderLayout.CENTER);

        table.getColumnModel().getColumn(1).setPreferredWidth(60);
        table.getColumnModel().getColumn(1).setMaxWidth(60);
        table.getColumnModel().getColumn(2).setPreferredWidth(60);
        table.getColumnModel().getColumn(2).setMaxWidth(60);
    }


    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Process Table Example");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setSize(400, 300);

            ProcessTableModel pmodel = new ProcessTableModel();
            ProcessTableView panel = new ProcessTableView(pmodel);
            frame.add(panel);

            // Create and add some tasks to the panel
            LongRunningTask task1 = new LongRunningTask("Task 1");
            LongRunningTask task2 = new LongRunningTask("Task 2");
            pmodel.addTask(task1);
            pmodel.addTask(task2);

            task1.execute();

            frame.setVisible(true);
        });
    }


}
