package com.idorsia.research.chem.hyperspace.gui;

import com.actelion.research.chem.*;
import com.actelion.research.gui.generic.GenericRectangle;
import com.idorsia.research.chem.hyperspace.fragment.FragmentGenerator;
import org.pushingpixels.substance.api.renderer.SubstanceDefaultTableCellRenderer;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;

public class SubstanceFragmentRenderer extends SubstanceDefaultTableCellRenderer {

    JPanel panel;
    StructureIcon icon;
    //JStructureView view;

    public SubstanceFragmentRenderer() {
        this.panel = new JPanel();
        this.panel.setLayout(new BorderLayout());
        this.icon = new StructureIcon(new StereoMolecule());
        this.setIcon(this.icon);
    }

    @Override
    public void setBounds(int x, int y, int width, int height) {
        super.setBounds(x, y, width, height);
        // ensure we update the size of the wrapped bar component
        icon.setSize(width, height);
    }

    @Override
    public void setBounds(Rectangle r) {
        super.setBounds(r);
        // ensure we update the size of the wrapped bar component
        icon.setSize((int) r.getWidth(), (int) r.getHeight());
    }

    IDCodeParser icp = new IDCodeParser();

    @Override
    protected void setValue(Object value) {
        try {
            if (value != null) {
                StereoMolecule mol = new StereoMolecule();
                if (value instanceof FragmentGenerator.Frag) {
                    // update the value of the view
                    //view.setIDCode(( (FragmentGenerator.Frag) value ).frag_idcode );
                    mol = new StereoMolecule(((FragmentGenerator.Frag) value).getFrag());
                    mol.ensureHelperArrays(Molecule.cHelperNeighbours);
                }
                if (value instanceof String) {
                    icp.parse(mol, (String) value);
                }
                if (value instanceof StereoMolecule) {
                    mol = (StereoMolecule) value;
                    mol.ensureHelperArrays(Molecule.cHelperNeighbours);
                }

                for (int zb = 0; zb < mol.getBonds(); zb++) {
                    mol.setBondBackgroundHiliting(zb, false);
                }
                for (int zb = 0; zb < mol.getBonds(); zb++) {
                    mol.setBondForegroundHiliting(zb, false);
                }
                icon.setMolecule(mol);

            }
        }
        catch(Exception ex) {
            System.out.println("exception in drawing code");
            ex.printStackTrace();
        }
    }


    /**
     * Simple custom component which renders a molecule.
     * Implements the Icon interface so it can be used where an Icon can.
     */
    public class StructureIcon implements Icon {

        private Depictor2D depictor;
        private String     text;
        private double value = 0.5d;

        public StructureIcon(StereoMolecule mi) {
            //this.depictor = new Depictor2D(mi);
            this.setMolecule(mi);
            this.text = "";
        }

        public void setText(String text) {
            this.text = text;
        }

        public void setMolecule(StereoMolecule mi) {
            this.depictor = new Depictor2D(mi);
            Color foreground = UIManager.getLookAndFeelDefaults().getColor( "Table.foreground" );
            if(foreground==null) {
                foreground = new Color(200,220,220);
            }
            this.depictor.setForegroundColor(foreground,new Color(0,0,0,0));
        }

        @Override
        public void paintIcon(Component c, Graphics g, int x, int y) {
//            Graphics2D g2 = (Graphics2D) g;
//            g2.setColor(Color.GREEN);
//            g2.fillRect(0, 0, (int) (getIconWidth() * value), getIconHeight());
            ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
            ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_RENDERING,RenderingHints.VALUE_RENDER_QUALITY);
            //depictor.validateView((Graphics2D)g,new Rectangle2D.Double(x,y,this.w,this.h), AbstractDepictor.cModeInflateToMaxAVBL);
            depictor.validateView((Graphics2D)g,new GenericRectangle(x,y,this.w,this.h), AbstractDepictor.cModeInflateToMaxAVBL);
            depictor.paint((Graphics2D)g);

            // draw text:
            if(this.text!=null) {
                if(this.text.length()>0) {
                    Font fa = g.getFont(); float fontsize_a = 10.0f;
                    for(int za=0;za<5;za++) {
                        Font fb = fa.deriveFont(fontsize_a);
                        FontMetrics fmi = g.getFontMetrics(fb);
                        int hgt = fmi.getHeight();
                        if( (1.0*hgt) / this.h <=0.3){break;}
                    }
                    g.drawString(this.text, (int)( x+4 ) , (int) (y + 0.65 * this.h) );
                }
            }

        }

        private int w;
        private int h;

        public void setSize(int w, int h) {
            this.w = w;
            this.h = h;
            //double ar = (1.0*w) / h;
            //double structure_ar = this.depictor.getBoundingRect().width / this.depictor.getBoundingRect().height;
        }

        @Override
        public int getIconWidth() {
            return this.w;
        }

        @Override
        public int getIconHeight() {
            return this.h;
        }

    }

}
