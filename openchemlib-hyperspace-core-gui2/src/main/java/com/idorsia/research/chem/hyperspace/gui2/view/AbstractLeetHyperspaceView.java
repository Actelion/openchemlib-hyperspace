package com.idorsia.research.chem.hyperspace.gui2.view;

import javax.swing.*;

public abstract class AbstractLeetHyperspaceView extends JPanel {
    public abstract void addResultsView(CombinatorialHitsView view);
    public abstract void addResultsView(RealTimeExpandingHitsView view);
}
