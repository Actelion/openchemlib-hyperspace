package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JStructureView;
import com.idorsia.research.chem.hyperspace.fragment.FragmentGenerator;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.event.CellEditorListener;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.util.EventObject;

public class FragmentRenderer implements TableCellRenderer, TableCellEditor, ListCellRenderer<String> {

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {

        // resolve object:
        String idcode = null;

        if(value instanceof String) {
            idcode = (String) value;
        }
        else {
            if (value instanceof FragmentGenerator.Frag) {
                FragmentGenerator.Frag fi = (FragmentGenerator.Frag) value;
                idcode = fi.frag_idcode;
            }
            if (value instanceof String) {
                idcode = (String) value;
            }
            if (value instanceof StereoMolecule) {
                idcode = ((StereoMolecule) value).getIDCode();
            }

            if (idcode != null) {
                return new JLabel("<<wrong object>>");
            }
        }


        //StereoMolecule mi = new StereoMolecule( fi.parent );
        JStructureView view = new JStructureView();
        view.setIDCode(idcode);

        if(hasFocus) {
            view.setBackground(Color.WHITE);
        }
        if(isSelected) {
            view.setBorder(new LineBorder(Color.red,1));
        }

        return view;
    }

    @Override
    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
        return this.getTableCellRendererComponent(table,value,isSelected,true,row,column);
    }

    @Override
    public Object getCellEditorValue() {
        return null;
    }

    @Override
    public boolean isCellEditable(EventObject anEvent) {
        return true;
    }

    @Override
    public boolean shouldSelectCell(EventObject anEvent) {
        return false;
    }

    @Override
    public boolean stopCellEditing() {
        return false;
    }

    @Override
    public void cancelCellEditing() {

    }

    @Override
    public void addCellEditorListener(CellEditorListener l) {

    }

    @Override
    public void removeCellEditorListener(CellEditorListener l) {

    }

    @Override
    public Component getListCellRendererComponent(JList<? extends String> list, String value, int index, boolean isSelected, boolean cellHasFocus) {
        StereoMolecule mi = new StereoMolecule();
        JStructureView view = new JStructureView();

        //view.setOpaqueBackground(false);
        view.setOpaque(false);

        JComponent ci = new JPanel(){
            @Override
            protected void paintComponent(Graphics g) {
                //super.paintComponent(g);
                GradientPaint gp = new GradientPaint(0,0, new Color(28,16,16) , this.getWidth(),this.getHeight() ,new Color(42,32,32) );
                Graphics2D g2 = (Graphics2D) g;
                g2.setPaint(gp);
                g2.fillRect(0,0,this.getWidth(),this.getHeight());
            }
        };
        ci.setOpaque(true);
        ci.setLayout(new BorderLayout());
        ci.add(view,BorderLayout.CENTER);

//        MoleculeDragAdapter drag_adapter = new MoleculeDragAdapter(view){
//            @Override
//            public Transferable getTransferable(Point point) {
//                return new MoleculeTransferable(mi);
//            }
//        };

        view.setIDCode(value);
        if(isSelected) {
            view.setForeground(Color.black.brighter());
            view.setBackground(new Color(220,240,255));
            view.setBorder(new LineBorder(Color.red.darker().darker(),1));
            //view.setBorder(new LineBorder(Color.cyan.brighter().brighter(),1));
        }
        else {
            view.setBorder( new LineBorder(new Color(200,76,10),1));
            //view.setForeground(Color.yellow.brighter());
            //view.setBackground(Color.black.brighter().brighter());
            //view.setBorder(new LineBorder(Color.cyan,1));
        }

        //return view;
        return ci;
    }



}
