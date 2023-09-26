package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.*;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import java.awt.*;

public class HitDataRenderer extends DefaultTableCellRenderer {

    //private JComponent pi = null;
    private IDCodeParser icp = new IDCodeParser();

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        //pi = (JComponent) super.getTableCellRendererComponent(table, null, isSelected, hasFocus, row, column);
        JPanel pi = new JPanel();

        // Set background and foreground colors based on selection
        if (isSelected) {
            pi.setBackground(table.getSelectionBackground());
            pi.setForeground(table.getSelectionForeground());
        } else {
            pi.setBackground(table.getBackground());
            pi.setForeground(table.getForeground());
        }

        pi.removeAll();
        pi.setLayout(new BorderLayout());
        if(value instanceof SynthonAssembler.ExpandedCombinatorialHit) {
            SynthonAssembler.ExpandedCombinatorialHit hit = (SynthonAssembler.ExpandedCombinatorialHit) value;
            JLabel plabel = new JLabel(hit.fragments.get(0).rxn_id);
            pi.add(plabel, BorderLayout.NORTH);

            JPanel pfrags = new JPanel();
            pfrags.setLayout(new GridLayout(1,hit.fragments.size()));
            for(int zi=0;zi<hit.fragments.size();zi++) {
                JPanel pfrag = new JPanel(); pfrag.setLayout(new BorderLayout());
                JPanel pfragbottom = new JPanel(); pfragbottom.setLayout(new FlowLayout(FlowLayout.CENTER));
                JLabel fragIdLabel = new JLabel( hit.fragments.get(zi).fragment_id );
                //fragIdLabel.setFont();
                pfragbottom.add( fragIdLabel );
                JPanel pfragcenter = new JPanel(); pfragcenter.setLayout(new BorderLayout());
                pfrag.add(pfragbottom,BorderLayout.SOUTH);
                pfrag.add(pfragcenter,BorderLayout.CENTER);
                StereoMolecule mi = new StereoMolecule();
                icp.parse(mi,hit.fragments.get(zi).idcode);
                JStructureView jsv = new JStructureView(mi);
                pfragcenter.add(jsv, BorderLayout.CENTER);
                pfrags.add(pfrag);
            }
            pi.add(pfrags, BorderLayout.CENTER);
            return pi;
        }
        else {
            return pi;
        }
    }
}
