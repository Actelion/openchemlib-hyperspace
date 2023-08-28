package com.idorsia.research.chem.hyperspace.gui2.view;

import com.actelion.research.gui.editor.EditorEvent;
import com.actelion.research.gui.generic.GenericEventListener;
import com.idorsia.research.chem.hyperspace.gui2.view.JExtendedDrawPanel;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;

import javax.swing.*;
import java.awt.*;

public class LeetHyperspaceView extends AbstractLeetHyperspaceView {


    private LeetHyperspaceModel model;

    private JPanel pSynthonSpaceView;
    private JPanel pProcssTableView;
    private JPanel pEditorView;
    private JPanel pCombiResultView;


    public LeetHyperspaceView(LeetHyperspaceModel model) {
        this.model = model;

        reinit();
    }

    private void reinit() {
        this.removeAll();
        this.setLayout(new GridLayout(2,2));

        pSynthonSpaceView = new JPanel(); pSynthonSpaceView.setLayout(new BorderLayout());
        pProcssTableView  = new JPanel(); pProcssTableView.setLayout(new BorderLayout());
        pEditorView       = new JPanel(); pEditorView.setLayout(new BorderLayout());
        pCombiResultView  = new JPanel(); pCombiResultView.setLayout(new BorderLayout());

        this.add(pSynthonSpaceView);
        this.add(pProcssTableView);
        this.add(pEditorView);
        this.add(pCombiResultView);

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
        this.pCombiResultView.removeAll();
        this.pCombiResultView.add(view, BorderLayout.CENTER);
        SwingUtilities.updateComponentTreeUI(this);
    }

    public void addResultsView(RealTimeExpandingHitsView view) {
        this.pCombiResultView.removeAll();
        this.pCombiResultView.add(view, BorderLayout.CENTER);
        SwingUtilities.updateComponentTreeUI(this);
    }


}
