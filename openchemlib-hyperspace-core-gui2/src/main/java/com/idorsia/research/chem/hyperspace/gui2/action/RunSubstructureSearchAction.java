package com.idorsia.research.chem.hyperspace.gui2.action;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.CombinatorialSearchResultModel;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.model.RealTimeExpandingSearchResultModel;
import com.idorsia.research.chem.hyperspace.gui2.task.SubstructureSearchTask;
import com.idorsia.research.chem.hyperspace.gui2.view.AbstractLeetHyperspaceView;
import com.idorsia.research.chem.hyperspace.gui2.view.CombinatorialHitsView;
import com.idorsia.research.chem.hyperspace.gui2.view.RealTimeExpandingHitsView;

import java.awt.event.ActionEvent;
import java.util.function.Supplier;

public class RunSubstructureSearchAction extends AbstractLeetHyperspaceAction {

    private SynthonSpace space;
    private Supplier<StereoMolecule> query;

    protected CombinatorialSearchResultModel resultModel;

    public RunSubstructureSearchAction(LeetHyperspaceModel model, AbstractLeetHyperspaceView view, SynthonSpace space, String spaceName, Supplier<StereoMolecule> query) {
        super("Substructure ["+spaceName+"]", null, model, view);
        this.space = space;
        this.query = query;
        this.resultModel = new CombinatorialSearchResultModel(new StereoMolecule());
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        StereoMolecule mq = new StereoMolecule(this.query.get());
        mq.ensureHelperArrays(Molecule.cHelperCIP);
        this.resultModel = new CombinatorialSearchResultModel(mq);
        RealTimeExpandingSearchResultModel expandingResultModel = new RealTimeExpandingSearchResultModel(resultModel,32000);
        CombinatorialHitsView resultView = new CombinatorialHitsView(resultModel);
        RealTimeExpandingHitsView resultView2 = new RealTimeExpandingHitsView(expandingResultModel);
        getView().addResultsView(resultView);
        getView().addResultsView(resultView2);
        SubstructureSearchTask ti = new SubstructureSearchTask(space,mq,resultModel);
        getModel().addTask(ti);
        ti.execute();
    }
}
