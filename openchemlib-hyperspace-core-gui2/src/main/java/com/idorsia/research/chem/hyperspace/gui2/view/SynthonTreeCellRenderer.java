package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.gui.table.ChemistryCellRenderer;
import com.actelion.research.gui.table.ChemistryRenderPanel;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.tree.SynthonTreeNode;

import javax.swing.*;
import javax.swing.tree.DefaultTreeCellRenderer;
import java.awt.*;

public class SynthonTreeCellRenderer extends DefaultTreeCellRenderer {
    @Override
    public Component getTreeCellRendererComponent(JTree tree, Object value,
                                                  boolean selected, boolean expanded,
                                                  boolean leaf, int row, boolean hasFocus) {

        if(value instanceof SynthonTreeNode) {
            SynthonSpace.FragId node = ((SynthonTreeNode) value).getSynthon();
            // Call super method to ensure default rendering behavior
            super.getTreeCellRendererComponent(tree, value, selected, expanded, leaf, row, hasFocus);
            // Customize rendering based on your requirements
            // Customize rendering based on the node's properties
            setText(node.fragment_id);

            // Custom painting for the icon
            ChemistryRenderPanel renderer = new ChemistryRenderPanel();
            renderer.setChemistry(HyperspaceUtils.parseIDCode(node.idcode));
            renderer.setSize(80, 80);
            Icon customIcon = new Icon() {
                @Override
                public void paintIcon(Component c, Graphics g, int x, int y) {
                    // Custom painting logic using Graphics object 'g'
                    // For example, draw a rectangle
                    g.setColor(Color.BLUE);
                    g.fillRect(x, y, getIconWidth(), getIconHeight());
                    renderer.paint(g);
                }

                @Override
                public int getIconWidth() {
                    return 20; // Specify the width of the custom icon
                }

                @Override
                public int getIconHeight() {
                    return 20; // Specify the height of the custom icon
                }
            };
            setIcon(customIcon);

            // Apply custom font, color, etc. if needed
            setFont(new Font("Arial", Font.BOLD, 12));
            setForeground(Color.RED);
            //}
        }
        else {
            super.getTreeCellRendererComponent(tree, value, selected, expanded, leaf, row, hasFocus);
        }

        return this;
    }
}