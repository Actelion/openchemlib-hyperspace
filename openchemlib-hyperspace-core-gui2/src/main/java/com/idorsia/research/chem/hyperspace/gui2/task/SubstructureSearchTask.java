package com.idorsia.research.chem.hyperspace.gui2.task;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.CombinatorialSearchResultModel;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class SubstructureSearchTask extends SwingWorker<List<SynthonSpace.CombinatorialHit>,List<SynthonSpace.CombinatorialHit>> {


    private SynthonSpace space;
    private StereoMolecule query;
    private CombinatorialSearchResultModel resultsModel;

    public SubstructureSearchTask(SynthonSpace space, StereoMolecule query, CombinatorialSearchResultModel resultsModel) {
        this.space = space;
        this.query = query;
        this.resultsModel = resultsModel;
    }

    @Override
    protected List<SynthonSpace.CombinatorialHit> doInBackground() throws Exception {
        SynthonSpace.CombinatorialHitReceiver receiver = new SynthonSpace.CombinatorialHitReceiver() {
            @Override
            public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hi, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                publish(hi);
            }
        };

        CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
        SubstructureSearchHelper.run_substructure_search_streaming_01_withBridgedBondsExpansion(this.space,cdp,query,Runtime.getRuntime().availableProcessors(),true,false,receiver);

        return new ArrayList<>();
    }

    private int nextChunk = 0;
    @Override
    protected void process(List<List<SynthonSpace.CombinatorialHit>> chunks) {
        for(int zi=nextChunk;zi<chunks.size();zi++) {
            List<SynthonSpace.CombinatorialHit> c_last = chunks.get(zi);
            if (!c_last.isEmpty()) {
                resultsModel.addResults(c_last);
            }
        }
        nextChunk = chunks.size();
    }

    public StereoMolecule getQuery() {
        return this.query;
    }

    public CombinatorialSearchResultModel getResultsModel() {
        return this.resultsModel;
    }

}
