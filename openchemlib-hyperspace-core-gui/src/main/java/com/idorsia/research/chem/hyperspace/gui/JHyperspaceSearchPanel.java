package com.idorsia.research.chem.hyperspace.gui;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.JDrawArea;
import com.actelion.research.gui.JDrawPanel;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.editor.GenericEditorArea;
import com.actelion.research.gui.editor.SwingEditorArea;
import com.actelion.research.gui.editor.SwingEditorPanel;
import com.idorsia.research.chem.hyperspace.gui.search.JSearchProviderListPanel;

import javax.swing.*;
import java.awt.*;

public class JHyperspaceSearchPanel extends JPanel {

    JSplitPane jsp_main;


    JPanel jp_draw;
    JExtendedDrawPanel jdp_draw;

    JPanel jp_left;
    JSearchProviderListPanel jp_searchproviders;

    public JHyperspaceSearchPanel() {
        initGUI();
    }

    private void initGUI() {
        this.removeAll();
        this.setLayout(new BorderLayout());
        jsp_main = new JSplitPane();
        this.add(jsp_main,BorderLayout.CENTER);

        this.jp_draw = new JPanel();
        this.jp_draw.setLayout(new BorderLayout());

        //this.jdp_draw = new JDrawPanel(new StereoMolecule());
        this.jdp_draw = new JExtendedDrawPanel();
        this.jp_draw.add(this.jdp_draw,BorderLayout.CENTER);
        this.jdp_draw.getDrawPanel().getDrawArea().setClipboardHandler(new ClipboardHandler());
        //this.jdp_draw.getDrawArea().setAllowQueryFeatures(true);
        //this.jdp_draw.getDrawArea().;

        // init paste idcode, paste smiles
        initDrawContextMenu(this.jdp_draw.getDrawPanel().getDrawArea());

        // init right panel and add search provider list
        this.jp_left = new JPanel();
        this.jp_left.setLayout(new BorderLayout());

        this.jp_searchproviders = new JSearchProviderListPanel();
        this.jp_left.add(this.jp_searchproviders,BorderLayout.CENTER);

        this.jsp_main.setRightComponent(this.jp_draw);
        this.jsp_main.setLeftComponent(this.jp_left);
        //this.jsp_main.setRightComponent(this.jp_searchproviders);

        //this.setBorder(BorderFactory.createTitledBorder("Editor"));
        this.jp_draw.setBorder(BorderFactory.createTitledBorder("Editor"));
    }

    public JSearchProviderListPanel getSearchProviderListPanel() {
        return this.jp_searchproviders;
    }

    public SwingEditorPanel getDrawPanel() {
        return this.jdp_draw.getDrawPanel();
    }


    public static void initDrawContextMenu(GenericEditorArea da) {
        //da.setComponentPopupMenu(new );
        //JPopupMenu pm_a = new JPopupMenu();
        //pm_a
    }


    public void restoreLayout() {

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                //jsp_main.setDividerLocation(0.2);
                UtilsGUI.setDividerLocation(jsp_main,0.4);
            }
        });
    }
}
