package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.table.AbstractTableModel;
import java.util.List;

public class CombiHitTableModel extends AbstractTableModel {

    private List<SynthonSpace.CombinatorialHit> hits;

    public CombiHitTableModel(List<SynthonSpace.CombinatorialHit> hits) {
        this.hits = hits;
    }

    public void setHits() {
        this.hits = hits;
        fireTableDataChanged();
    }

    public void addHits(List<SynthonSpace.CombinatorialHit> newHits) {
        this.hits.addAll(newHits);
        fireTableRowsInserted(0,this.hits.size());
        fireTableDataChanged();
    }

    @Override
    public int getRowCount() {
        return hits.size();
    }

    @Override
    public int getColumnCount() {
        return 3;
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        return this.hits.get(rowIndex).rxn;
    }
}
