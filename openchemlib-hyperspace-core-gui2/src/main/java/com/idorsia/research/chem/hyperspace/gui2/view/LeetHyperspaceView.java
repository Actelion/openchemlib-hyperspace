package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.gui.editor.EditorEvent;
import com.actelion.research.gui.generic.GenericEventListener;
import com.idorsia.research.chem.hyperspace.gui2.view.JExtendedDrawPanel;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;

import javax.swing.*;
import java.awt.*;

public class LeetHyperspaceView extends AbstractLeetHyperspaceView {


    private LeetHyperspaceModel model;

    private JPanel pLeft;
    private JPanel pRight;
    private JTabbedPane tpLeftTop;
    private JPanel pLeftBottom;

    private JPanel pSynthonSpaceView;
    private JPanel pProcssTableView;
    private JPanel pEditorView;
    private JTabbedPane tpCombiResultView;
    private JPanel pCombiResultViewA;
    private JPanel pCombiResultViewB;


    public LeetHyperspaceView(LeetHyperspaceModel model) {
        this.model = model;

        reinit();
    }

    private void reinit() {
        this.removeAll();
        this.setLayout(new GridLayout(1,2));
        //this.setLayout(new BorderLayout());


        pLeft = new JPanel(); pLeft.setLayout(new BorderLayout());
        pRight = new JPanel(); pRight.setLayout(new BorderLayout());
        tpLeftTop = new JTabbedPane();
        pLeftBottom = new JPanel(); pLeftBottom.setLayout(new BorderLayout());

        this.add(pLeft);
        this.add(pRight);

        pLeft.add(tpLeftTop,BorderLayout.NORTH);
        pLeft.add(pLeftBottom,BorderLayout.CENTER);

        pSynthonSpaceView = new JPanel(); pSynthonSpaceView.setLayout(new BorderLayout());
        pProcssTableView  = new JPanel(); pProcssTableView.setLayout(new BorderLayout());
        pEditorView       = new JPanel(); pEditorView.setLayout(new BorderLayout());
        tpCombiResultView = new JTabbedPane();
        pCombiResultViewA  = new JPanel(); pCombiResultViewA.setLayout(new BorderLayout());
        pCombiResultViewB  = new JPanel(); pCombiResultViewB.setLayout(new BorderLayout());
        tpCombiResultView.addTab("Combinatorial",pCombiResultViewA);
        tpCombiResultView.addTab("Enumerated",pCombiResultViewB);

        tpLeftTop.addTab("Editor",pEditorView);
        tpLeftTop.addTab("Synthon Space",pSynthonSpaceView);
        pLeftBottom.add(pProcssTableView,BorderLayout.CENTER);
        //pRight.add(pCombiResultView,BorderLayout.CENTER);
        pRight.add(tpCombiResultView,BorderLayout.CENTER);

        //this.add(pSynthonSpaceView);
        //this.add(pProcssTableView);
        //this.add(pEditorView);
        //this.add(pCombiResultView);

        ProcessTableView processTableView = new ProcessTableView(model.getProcessTableModel());
        this.pProcssTableView.add(processTableView,BorderLayout.CENTER);

        SynthonSpaceView synthonSpaceView = new SynthonSpaceView(model);
        this.pSynthonSpaceView.add(synthonSpaceView,BorderLayout.CENTER);

        JExtendedDrawPanel pDrawPanel = new JExtendedDrawPanel();
        this.pEditorView.add(pDrawPanel,BorderLayout.CENTER);
        pDrawPanel.getDrawPanel().getDrawArea().addDrawAreaListener(new GenericEventListener<EditorEvent>() {
            @Override
            public void eventHappened(EditorEvent editorEvent) {
                model.setQuery(pDrawPanel.getDrawPanel().getDrawArea().getMolecule());
            }
        });

        SwingUtilities.updateComponentTreeUI(this);
    }

    public void addResultsView(CombinatorialHitsView view) {
        this.pCombiResultViewA.removeAll();
        this.pCombiResultViewA.add(view, BorderLayout.CENTER);
        SwingUtilities.updateComponentTreeUI(this);
    }

    public void addResultsView(RealTimeExpandingHitsView view) {
        this.pCombiResultViewB.removeAll();
        this.pCombiResultViewB.add(view, BorderLayout.CENTER);
        SwingUtilities.updateComponentTreeUI(this);
    }


}
