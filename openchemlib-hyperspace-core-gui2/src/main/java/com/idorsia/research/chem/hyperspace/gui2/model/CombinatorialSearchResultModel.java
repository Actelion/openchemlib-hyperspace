package com.idorsia.research.chem.hyperspace.gui2.model;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.List;

public class CombinatorialSearchResultModel {

    private StereoMolecule query;
    private List<SynthonSpace.CombinatorialHit> hits = new ArrayList<>();
    private List<CombinatorialSearchResultModelListener> listeners = new ArrayList<>();

    public CombinatorialSearchResultModel(StereoMolecule query) {
        this.query = query;
    }

    public StereoMolecule getQuery() {
        return this.query;
    }

    public synchronized void addResults(List<SynthonSpace.CombinatorialHit> hits1) {
        this.hits.addAll(hits1);
        System.out.println("[CombinatorialSearchResultModel] hits: "+hits1.size());
        this.tableModel.addHits(hits1);
        fireResultsChanged();
    }

    public void addListener(CombinatorialSearchResultModelListener li) {
        this.listeners.add(li);
    }

    public boolean removeListener(CombinatorialSearchResultModelListener li) {
        return this.listeners.remove(li);
    }

    private void fireResultsChanged() {
        for(CombinatorialSearchResultModelListener li : listeners) {
            li.resultsChanged();
        }
    }

    public static interface CombinatorialSearchResultModelListener {
        public void resultsChanged();
    }

    private CombiHitTableModel tableModel = new CombiHitTableModel(new ArrayList<>());

    public CombiHitTableModel getTableModel() {
        return tableModel;
    }

    public String getNumberOfResultsString() {
        List<SynthonSpace.CombinatorialHit> hits_i = new ArrayList<>(hits);
        long sum_i = hits_i.stream().mapToLong( xi -> xi.hit_fragments.values().stream().mapToLong( yi -> yi.size() ).reduce((x,y)->x*y).getAsLong()).sum();
        //if(sum_i < 4096) {return ""+sum_i;}
        if(sum_i < 20000) {return ""+sum_i;}
        if(sum_i < 100000) {return " > 10000";}
        if(sum_i < 1000000) {return " > 100000";}
        if(sum_i < 4000000) {return " > 1000000";}
        return " >> 1000000";
    }

    public synchronized List<SynthonSpace.CombinatorialHit> getHits() {
        return new ArrayList<>(this.hits);
    }

}
