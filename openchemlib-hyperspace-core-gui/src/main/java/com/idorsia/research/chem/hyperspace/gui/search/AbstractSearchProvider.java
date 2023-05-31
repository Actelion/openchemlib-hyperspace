package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;

import javax.swing.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



public abstract class AbstractSearchProvider<T> {

    public static interface SearchConfiguration {
        public List<StereoMolecule> getQueryMolecules();
        public String serializeToJSON();
        public void   deserializeFromJSON(String json) throws Exception;
        public int    getNumberOfThreads();
        public void   setNumberOfThreads(int threads);
    }

    public static interface SearchProviderConfiguration {
        public String getServiceName();
        public void setServiceName(String serviceName);
        public String serializeToJSON();
        public void deserializeFromJSON(String json) throws Exception;
    }

    public static class IncompatibleSearchConfigurationException extends Exception {
    }

    public abstract void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui);

    public enum SearchProviderStatus { REMOTE, NOT_INITIALIZED , INITIALIZING , READY , ERROR };


    public abstract String getSearchName();

    public abstract String getSpaceName();

    public abstract int getMaxNumberOfThreads();

    public abstract JPanel getSearchConfigurationView();

    public abstract SearchConfiguration getSearchConfiguration();

    public abstract Class getSearchConfigurationClass();

    public abstract AbstractHyperspaceProcess startInitialization();

    public abstract AbstractHyperspaceSearchProcess<T> runSearch(SearchConfiguration configuration ) throws IncompatibleSearchConfigurationException;

    public abstract void presentSearchResults(List<T> results);

    public abstract SearchProviderStatus getStatus();

    public abstract double getStatusOfInitialization();

    public abstract RemoteSearchProviderAdapter<T> getRemoteSearchProvider();

    public abstract void setRemoteSearchProvider(RemoteSearchProviderAdapter<T> remote_adapter);

    public abstract String encodeResultsToString(List<T> results) throws IOException;

    public abstract List<T> decodeResultsFromString(String results_string) throws IOException;

    public abstract boolean hasGUI();

    public static interface SearchProviderListener {
        public void statusChanged();
        public void requestViewUpdate();
    }

    List<SearchProviderListener> listeners = new ArrayList<>();

    protected void fireSearchStatusChanged() {
        for(SearchProviderListener li : listeners) {
            li.statusChanged();
        }
    }

    public abstract SearchProviderConfiguration getSearchProviderConfiguration();

    public void addSearchProviderListener(SearchProviderListener li) {
        this.listeners.add(li);
    }
    public void removeSearchProviderListener(SearchProviderListener li) {
        this.listeners.remove(li);
    }

}
