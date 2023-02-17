package com.idorsia.research.chem.hyperspace.gui.process;

import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;

import java.util.List;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



public abstract class AbstractHyperspaceSearchProcess<T> extends AbstractHyperspaceProcess {

    private AbstractSearchProvider<T> mSearchProvider = null;

    public AbstractHyperspaceSearchProcess(AbstractSearchProvider<T> search_provider) {
        this.mSearchProvider = search_provider;
    }

    public AbstractSearchProvider<T> getSearchProvider() {return this.mSearchProvider;}

    public abstract void setSearchConfiguration(AbstractSearchProvider.SearchConfiguration config);

    public abstract AbstractSearchProvider.SearchConfiguration getSearchConfiguration();

    /**
     * This function should start the search in a separate thread,
     * Implementations should call the super function that adds a
     * process status listener that automatically calls the
     * presentSearchResults function in case that this search
     * provider has a GUI.
     */
    public void startSearchAsync() {
        this.addSearchProviderListener(new AbstractHyperspaceProcess.HyperspaceProcessListener() {
            @Override
            public void processStatusChanged() {
                if(getProcessStatus()== AbstractHyperspaceProcess.ProcessStatus.DONE) {
                    if( mSearchProvider.hasGUI() ) {
                        mSearchProvider.presentSearchResults(getSearchResults());
                    }
                    else {
                        System.out.println("Completed search request");
                        System.out.println("Request: \n"+ getSearchConfiguration().serializeToJSON());
                        System.out.println("Results: "+getSearchResults().size());
                    }
                }
            }
        });
    }

    public abstract List<T> getSearchResults();

}
