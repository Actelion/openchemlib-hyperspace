package com.idorsia.research.synthonbrowser;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;

import javax.swing.table.AbstractTableModel;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

final class SynthonTableModel extends AbstractTableModel {

    record Row(String reactionId, int setIndex, RawSynthon synthon) {}

    private static final String[] COLUMNS = {
            "Structure[idcode]",
            "Fragment ID",
            "Connectors",
            "IDCode",
            "Reaction",
            "Set"
    };

    private final List<Row> rows = new ArrayList<>();

    @Override
    public int getRowCount() {
        return rows.size();
    }

    @Override
    public int getColumnCount() {
        return COLUMNS.length;
    }

    @Override
    public String getColumnName(int column) {
        return COLUMNS[column];
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        if (columnIndex == 5) {
            return Integer.class;
        }
        return String.class;
    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return false;
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        return getValueAt(rows.get(rowIndex), columnIndex);
    }

    Object getValueAt(Row row, int columnIndex) {
        RawSynthon synthon = row.synthon();
        return switch (columnIndex) {
            case 0 -> synthon.getIdcode();
            case 1 -> synthon.getFragmentId();
            case 2 -> formatConnectors(synthon.getConnectors());
            case 3 -> synthon.getIdcode();
            case 4 -> row.reactionId();
            case 5 -> row.setIndex();
            default -> "";
        };
    }

    void setRows(List<Row> newRows) {
        rows.clear();
        rows.addAll(newRows);
        fireTableDataChanged();
    }

    Row getRow(int modelRow) {
        return rows.get(modelRow);
    }

    static String formatConnectors(BitSet connectors) {
        if (connectors == null || connectors.isEmpty()) {
            return "-";
        }
        List<Integer> values = new ArrayList<>();
        for (int bit = connectors.nextSetBit(0); bit >= 0; bit = connectors.nextSetBit(bit + 1)) {
            values.add(bit);
        }
        return values.stream().map(String::valueOf).collect(Collectors.joining(","));
    }
}
