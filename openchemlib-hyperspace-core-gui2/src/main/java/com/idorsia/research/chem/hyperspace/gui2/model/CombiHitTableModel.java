package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import org.apache.commons.lang3.tuple.Pair;

import javax.swing.table.AbstractTableModel;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
        return 4;
    }

    public SynthonSpace.CombinatorialHit getHit(int idx) {
        return this.hits.get(idx);
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        SynthonSpace.CombinatorialHit hit = this.hits.get(rowIndex);
        List<SynthonSpace.FragType> sortedFragTypes = hit.hit_fragments.keySet().stream().sorted().collect(Collectors.toList());
        switch(columnIndex) {
            case 0:
                long[] numSynthons = hit.hit_fragments.values().stream().mapToLong( xi -> xi.size()).toArray();
                String hitsStr = "["+ Arrays.stream(numSynthons).mapToObj(xi -> ""+xi).collect(Collectors.joining("x"))+"]";
                String hitsStr2 = ""+Arrays.stream(numSynthons).reduce((x,y)->x*y).getAsLong();
                return hit.rxn+" Hits: "+hitsStr2+ " "+hitsStr;
            default:
                int frag_idx = columnIndex-1;
                if(frag_idx < sortedFragTypes.size()) {
                    SynthonSpace.FragType fti = sortedFragTypes.get(frag_idx);
                    List<Map.Entry<Integer, Pair<SynthonSpace.FragType, BitSet>>> matched_i = hit.mapping.entrySet().stream().filter(xi -> xi.getValue().getLeft().equals(fti)).collect(Collectors.toList());
                    if (matched_i.size() > 0) {
                        return hit.sri.fragments[matched_i.get(0).getKey()].getIDCode();
                    } else {
                        return "";
                    }
                }
                else {
                    return "";
                }
        }
    }
}
