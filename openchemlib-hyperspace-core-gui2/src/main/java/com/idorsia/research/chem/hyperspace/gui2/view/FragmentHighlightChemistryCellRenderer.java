package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.table.ChemistryCellRenderer;
import com.actelion.research.gui.table.ChemistryRenderPanel;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import javax.swing.*;
import java.awt.*;

public class FragmentHighlightChemistryCellRenderer extends ChemistryCellRenderer {

    public FragmentHighlightChemistryCellRenderer(Dimension preferredSize, StereoMolecule frag) {
        super(preferredSize);
        this.frag = frag;
    }

    private StereoMolecule frag = null;

    public void setHighlightedFragment(StereoMolecule frag) {
        this.frag = frag;
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {
        if(frag==null) {
            return super.getTableCellRendererComponent(table, value, false, hasFocus, row, col);
        }

        StereoMolecule mi = HyperspaceUtils.parseIDCode((String) value);
        HyperspaceUtils.setHighlightedSubstructure(mi,this.frag,false);

        return super.getTableCellRendererComponent(table, mi, false, hasFocus, row, col);
    }
}
