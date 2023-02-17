package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpaceExplorer3;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */


public class HyperspaceSimilaritySearchProcess extends AbstractHyperspaceSearchProcess<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> implements AbstractHyperspaceProcess.HasProgress {


    private SynthonSimilaritySpace3 space3;

    private HyperspaceSimilaritySearch.SimilaritySearchConfiguration config = null;
    //private StereoMolecule query;

    private List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> results;

    public HyperspaceSimilaritySearchProcess(HyperspaceSimilaritySearch search, SynthonSimilaritySpace3 space3) {
        super(search);
        this.space3 = space3;
    }


    @Override
    public double getProgress() {
        return this.progress;
    }



    @Override
    public void setSearchConfiguration(AbstractSearchProvider.SearchConfiguration config) {
        this.config = (HyperspaceSimilaritySearch.SimilaritySearchConfiguration) config;
        //this.query = config.getQueryMolecules().get(0);
    }

    @Override
    public void startSearchAsync() {
        super.startSearchAsync();
        setProcessStatus(ProcessStatus.COMPUTING);

        SimilaritySearchThread search_thread = new SimilaritySearchThread();
        search_thread.start();
    }

    @Override
    public AbstractSearchProvider.SearchConfiguration getSearchConfiguration() {
        return this.config;
    }

    @Override
    public List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> getSearchResults() {
        return this.results;
    }

    @Override
    public String getName() {
        return null;
    }


    public class SimilaritySearchThread extends Thread {
        @Override
        public void run() {

            int number_of_threads = getSearchProvider().getSearchConfiguration().getNumberOfThreads();
            // check if this gets overruled by the max number of threads that we can use
            if( getSearchProvider().getMaxNumberOfThreads() > 0 ) {
                number_of_threads = Math.min(number_of_threads, getSearchProvider().getMaxNumberOfThreads());
            }
            if(number_of_threads<=0) { number_of_threads = Runtime.getRuntime().availableProcessors(); }



            CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
            SynthonSimilaritySpaceExplorer3 explorer = new SynthonSimilaritySpaceExplorer3(space3,cdp);
            SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3 search_config = new SynthonSimilaritySpaceExplorer3.SimilaritySearchConfig3(number_of_threads,3 ,3 ,0.8 ,4000);

            Consumer<Double> f_progress = new Consumer<Double>() {
                @Override
                public void accept(Double aDouble) {
                    setProgressAndFireUpdate(0.95*aDouble);
                }
            };

            StereoMolecule query = getSearchConfiguration().getQueryMolecules().get(0);
            List<Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,Double>> next_results = explorer.runSearch(query,search_config,f_progress);

            List<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult> results_a = new ArrayList<>();
            for(Pair<SynthonSimilaritySpaceExplorer3.FinalSimilarityResult,Double> fri : next_results) {
                results_a.add(fri.getLeft());
            }

            results = results_a;

            setProcessStatusMessage("");
            setProcessStatus(ProcessStatus.DONE);
            setProgressAndFireUpdate(1.0);
        }
    }

    private double progress = 0.0;

    private void setProgressAndFireUpdate(double progress) {
        this.progress = progress;
        this.fireProcessStatusChanged();
    }


}
