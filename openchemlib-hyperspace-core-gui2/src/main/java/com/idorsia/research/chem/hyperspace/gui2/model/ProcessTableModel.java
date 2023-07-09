package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.gui2.task.SubstructureSearchTask;

import javax.swing.*;
import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class ProcessTableModel {

    private TableModel tableModel = new TableModel();
    private TableModelOnlySearches tableModelSearches = new TableModelOnlySearches();
    private List<SwingWorker> taskList = new ArrayList<>();

    public void addTask(SwingWorker task) {
        this.taskList.add(task);
        task.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                tableModel.fireTableChanged(new TableModelEvent(tableModel));
                tableModelSearches.fireTableChanged(new TableModelEvent(tableModel));
            }
        });
        tableModel.fireTableChanged(new TableModelEvent(this.tableModel));
        tableModelSearches.fireTableChanged(new TableModelEvent(tableModel));
    }

    public TableModel getTableModel() {
        return this.tableModel;
    }

    public TableModelOnlySearches getPastSearchesTableModel() {return this.tableModelSearches;}

    public class TableModel extends AbstractTableModel {
        @Override
        public int getRowCount() {
            return taskList.size();
        }

        @Override
        public int getColumnCount() {
            return 3;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SwingWorker ti = taskList.get(rowIndex);
            switch(columnIndex) {
                case 0: return "Task";
                case 1: return ti.getState();
                case 2: return ti.getProgress();
            }
            return null;
        }
    }

    public List<SubstructureSearchTask> getPastSearches() {
        return new ArrayList<>( this.taskList.stream().filter( xi -> xi instanceof SubstructureSearchTask ).map( xi -> (SubstructureSearchTask) xi).collect(Collectors.toList()) );
    }

    public class TableModelOnlySearches extends AbstractTableModel {
        @Override
        public int getRowCount() {
            return getPastSearches().size();
        }

        @Override
        public int getColumnCount() {
            return 3;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SubstructureSearchTask ti = getPastSearches().get(rowIndex);
            switch(columnIndex) {
                case 0: return ti.getQuery();
                case 1: return ti.getState();
                case 2: return ti.getResultsModel().getNumberOfResultsString();
            }
            return null;
        }
    }

}
