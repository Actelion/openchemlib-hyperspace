package com.idorsia.research.chem.hyperspace.gui2.action;

import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.view.AbstractLeetHyperspaceView;

import javax.swing.*;

public abstract class AbstractLeetHyperspaceAction extends AbstractAction  {

    private LeetHyperspaceModel model;
    private AbstractLeetHyperspaceView  view;

    public AbstractLeetHyperspaceAction(String name, Icon icon, LeetHyperspaceModel model, AbstractLeetHyperspaceView view) {
        super(name, icon);
        this.model = model;
        this.view = view;
    }

    public LeetHyperspaceModel getModel() {
        return model;
    }

    public AbstractLeetHyperspaceView getView() {
        return view;
    }
}
