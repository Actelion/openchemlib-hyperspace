package com.idorsia.research.chem.hyperspace.gui2.model;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class LeetHyperspaceModel {

    private List<SynthonSpace> spaces = new ArrayList<>();

    private ProcessTableModel  processTableModel;

    private SynthonSpacesTreeModel synthonSpacesTreeModel;

    private List<LeetHyperspaceModelListener> listeners;

    private StereoMolecule query;

    public LeetHyperspaceModel() {
        processTableModel = new ProcessTableModel();
        listeners = new ArrayList<>();
        synthonSpacesTreeModel = new SynthonSpacesTreeModel(this);
    }

    public void addSynthonSpace(SynthonSpace si) {
        this.spaces.add(si);
        fireSynthonSpacesChanged();
    }


    public List<SynthonSpace> getSynthonSpaces() {
        return this.spaces;
    }

    public SynthonSpacesTreeModel getSynthonSpacesTreeModel() {
        return this.synthonSpacesTreeModel;
    }

    public void addTask(SwingWorker task) {
        this.processTableModel.addTask(task);
    }

    public void setQuery(StereoMolecule qi) {
        this.query = qi;
        fireQueryMoleculeChanged();
    }

    public StereoMolecule getQuery() {
        return this.query;
    }

    public ProcessTableModel getProcessTableModel() {
        return this.processTableModel;
    }

    public void addListener(LeetHyperspaceModelListener li) {
        this.listeners.add(li);
    }

    public boolean removeListener(LeetHyperspaceModelListener li) {
        return this.listeners.remove(li);
    }

    private void fireSynthonSpacesChanged() {
        for(LeetHyperspaceModelListener li : listeners) {
            li.synthonSpacesChanged();
        }
    }
    private void fireQueryMoleculeChanged() {
        for(LeetHyperspaceModelListener li : listeners) {
            li.queryMoleculeChanged();
        }
    }

    public static interface LeetHyperspaceModelListener {
        default public void synthonSpacesChanged() {
        }

        default public void queryMoleculeChanged() {
        }
    }
}
