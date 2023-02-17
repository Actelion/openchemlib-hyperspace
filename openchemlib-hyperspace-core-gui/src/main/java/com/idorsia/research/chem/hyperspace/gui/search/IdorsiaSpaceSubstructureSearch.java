package com.idorsia.research.chem.hyperspace.gui.search;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;

import javax.swing.*;
import java.io.IOException;
import java.util.List;

public class IdorsiaSpaceSubstructureSearch extends AbstractSearchProvider<SynthonSpace.CombinatorialHit> {


    HyperspaceSubstructureSearch.InitializationConfig config_initialization;
    HyperspaceSynthonSpaceInitializationProcess process_initialization;



    HyperspaceSearchGUI gui;

    public IdorsiaSpaceSubstructureSearch(HyperspaceSearchGUI gui, HyperspaceSubstructureSearch.InitializationConfig config) {
        this.config_initialization = config;
        this.gui = gui;
    }

    @Override
    public String getSearchName() {
        return "Hyperspace SS";
    }

    @Override
    public String getSpaceName() {
        return "Idorsia VCS";
    }

    @Override
    public JPanel getSearchConfigurationView() {
        return null;
    }

    @Override
    public boolean hasGUI() {
        return true;
    }

    @Override
    public SearchConfiguration getSearchConfiguration() {
        return null;
    }

    @Override
    public int getMaxNumberOfThreads() {
        return -1;
    }

    @Override
    public Class getSearchConfigurationClass() {
        return null;
    }

    @Override
    public void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui) {

    }

    @Override
    public AbstractHyperspaceProcess startInitialization() {
        process_initialization = new HyperspaceSynthonSpaceInitializationProcess(this,config_initialization.getFile()); //new HyperspaceSubstructureSearch.InitializeHyperspaceSynthonSpaceProcess(config_initialization.getFile());
        process_initialization.startInitializationAsync();
        return process_initialization;
    }

    @Override
    public AbstractHyperspaceSearchProcess<SynthonSpace.CombinatorialHit> runSearch(SearchConfiguration configuration) throws IncompatibleSearchConfigurationException {
        return null;
    }

    @Override
    public void presentSearchResults(List<SynthonSpace.CombinatorialHit> results) {

    }

    @Override
    public String encodeResultsToString(List<SynthonSpace.CombinatorialHit> results) throws IOException {
        return null;
    }

    @Override
    public RemoteSearchProviderAdapter<SynthonSpace.CombinatorialHit> getRemoteSearchProvider() {
        return null;
    }

    @Override
    public void setRemoteSearchProvider(RemoteSearchProviderAdapter<SynthonSpace.CombinatorialHit> remote_adapter) {

    }

    @Override
    public List<SynthonSpace.CombinatorialHit> decodeResultsFromString(String results_string) throws IOException {
        return null;
    }

    @Override
    public SearchProviderStatus getStatus() {
        return null;
    }

    @Override
    public double getStatusOfInitialization() {
        return 0;
    }
}
