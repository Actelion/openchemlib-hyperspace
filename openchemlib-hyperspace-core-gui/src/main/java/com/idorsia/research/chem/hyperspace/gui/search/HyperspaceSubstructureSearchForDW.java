
package com.idorsia.research.chem.hyperspace.gui.search;

import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.actelion.research.chem.hyperspace.SimpleCombinatorialHit;
import com.actelion.research.chem.hyperspace.SimpleSynthon;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;


import javax.swing.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Base64;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Idorsia Pharmaceuticals Ltd. 2022
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */



/**
 * This class mimics the HyperspaceSubstructureSearch, but it uses a simplified
 * combinatorial hit object for easy transfer to DataWarrior.
 *
 * Some remarks on implementing DW services that encapsulate existing search providers:
 * 1. for the search process implemented via member variable, remember to listen and forward state changes!
 *
 */
public class HyperspaceSubstructureSearchForDW extends AbstractSearchProvider<SimpleCombinatorialHit> {


    public static List<com.actelion.research.chem.hyperspace.SimpleCombinatorialHit> decodeSimpleCombinatorialHits(String hits_string) throws IOException {
        byte[] data_serialization = Base64.getDecoder().decode(hits_string);
        ObjectInputStream in      = new ObjectInputStream( new ByteArrayInputStream(data_serialization));
        List<SimpleCombinatorialHit> results = null;
        try {
            results = (List<SimpleCombinatorialHit>)  in.readObject();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            System.out.println("[decodeResultsFromString] :: Class not found while deserializing..");
            return null;
        }
        return results;
    }
    public static String encodeSimpleCombinatorialHits(List<SimpleCombinatorialHit> hits) throws IOException {
        ByteArrayOutputStream out_bytes = new ByteArrayOutputStream();
        ObjectOutputStream out = new ObjectOutputStream( new BufferedOutputStream( out_bytes ) );
        out.writeObject(hits);
        out.flush();
        out_bytes.flush();
        byte[] data_serialization = out_bytes.toByteArray();
        String b64_data           = Base64.getEncoder().encodeToString(data_serialization);
        return b64_data;
    }


    // This is the object that we use to perform all computations
    private HyperspaceSubstructureSearch search = new HyperspaceSubstructureSearch();
    public HyperspaceSubstructureSearch getSubstructureSearch() {
        return this.search;
    }

    private RemoteSearchProviderAdapter<SimpleCombinatorialHit> remoteSearchProviderAdapter = null;

    public static class HyperspaceSubstructureSearchProcessForDW extends AbstractHyperspaceSearchProcess<SimpleCombinatorialHit> implements AbstractHyperspaceProcess.HasProgress {
        HyperspaceSubstructureSearchProcess process;

        public HyperspaceSubstructureSearchProcessForDW(HyperspaceSubstructureSearchForDW search, SynthonSpace space, CachedDescriptorProvider cdp) {
            super(search);
            process = new HyperspaceSubstructureSearchProcess(search.getSubstructureSearch(),space,cdp);
            process.addSearchProviderListener(new HyperspaceProcessListener() {
                @Override
                public void processStatusChanged() {
                    fireProcessStatusChanged();
                }
            });
        }

        @Override
        public String getName() {
            return process.getName()+"_forDW";
        }

        @Override
        public double getProgress() {
            return process.getProgress();
        }

        @Override
        public ProcessStatus getProcessStatus() {
            return process.getProcessStatus();
        }

        @Override
        public String getProcessStatusMessage() {
            return process.getProcessStatusMessage();
        }

        @Override
        public void setSearchConfiguration(SearchConfiguration config) {
            process.setSearchConfiguration(config);
        }

        @Override
        public SearchConfiguration getSearchConfiguration() {
            return process.getSearchConfiguration();
        }

        @Override
        public List<SimpleCombinatorialHit> getSearchResults() {
            // here we have to convert the results:
            List<SynthonSpace.CombinatorialHit> hits = process.getSearchResults();
            return new ArrayList<>( hits.stream().map( ci -> createSimpleCombinatorialHit(ci) ).collect(Collectors.toList()));
        }

        @Override
        public void startSearchAsync() {
            process.startSearchAsync();
        }
    }

    public static SimpleCombinatorialHit createSimpleCombinatorialHit(SynthonSpace.CombinatorialHit ci) {
        SimpleSynthon[][] simple_synthons = new SimpleSynthon[ci.hit_fragments.keySet().size()][];
        List<List<SynthonSpace.FragId>>  frags = new ArrayList<>(ci.hit_fragments.values());
        for(int zi=0;zi<frags.size();zi++) {
            simple_synthons[zi] = new SimpleSynthon[frags.get(zi).size()];
            for(int zj=0;zj<frags.get(zi).size();zj++) {
                SynthonSpace.FragId fid = frags.get(zi).get(zj);
                simple_synthons[zi][zj] = new SimpleSynthon( fid.idcode , fid.fragment_id, fid.rxn_id , ""+fid.frag );
            }
        }
        return new SimpleCombinatorialHit(ci.rxn,simple_synthons);
    }

    @Override
    public void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui) {
        search.setConfigurationAndGUI(config,gui);
    }

    @Override
    public String getSearchName() {
        return search.getSearchName()+"_ForDataWarrior";
    }

    @Override
    public String getSpaceName() {
        return search.getSpaceName();
    }

    @Override
    public int getMaxNumberOfThreads() {
        return search.getMaxNumberOfThreads();
    }

    @Override
    public JPanel getSearchConfigurationView() {
        return search.getSearchConfigurationView();
    }

    @Override
    public SearchConfiguration getSearchConfiguration() {
        return search.getSearchConfiguration();
    }

    @Override
    public SearchProviderConfiguration getSearchProviderConfiguration() {
        return this.remoteSearchProviderAdapter.getSearchProviderConfiguration();
    }

    @Override
    public Class getSearchConfigurationClass() {
        return search.getSearchConfigurationClass();
    }

    @Override
    public AbstractHyperspaceProcess startInitialization() {
        return search.startInitialization();
    }

    @Override
    public AbstractHyperspaceSearchProcess<SimpleCombinatorialHit> runSearch(SearchConfiguration configuration) throws IncompatibleSearchConfigurationException {
        CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
//        if(gui != null) {
//            p_results = new JPanel();
//            gui.getHyperspaceMainPanel().addAnalysisTab("Substructure Search", p_results);
//        }
        HyperspaceSubstructureSearchProcessForDW search_process = new HyperspaceSubstructureSearchProcessForDW(this, this.getSubstructureSearch().space, cdp);
        search_process.setSearchConfiguration(configuration);
        search_process.startSearchAsync();

        return search_process;
    }

    @Override
    public void presentSearchResults(List<SimpleCombinatorialHit> results) {
        // nothing to implement, this is just a server service provider
    }

    @Override
    public SearchProviderStatus getStatus() {
        return this.search.getStatus();
    }

    @Override
    public double getStatusOfInitialization() {
        return this.search.getStatusOfInitialization();
    }

    @Override
    public RemoteSearchProviderAdapter<SimpleCombinatorialHit> getRemoteSearchProvider() {
        return this.remoteSearchProviderAdapter;
    }

    @Override
    public void setRemoteSearchProvider(RemoteSearchProviderAdapter<SimpleCombinatorialHit> remote_adapter) {
        this.remoteSearchProviderAdapter = remote_adapter;
    }

    @Override
    public String encodeResultsToString(List<SimpleCombinatorialHit> results) throws IOException {
        return encodeSimpleCombinatorialHits(results);
    }

    @Override
    public List<SimpleCombinatorialHit> decodeResultsFromString(String results_string) throws IOException {
        return decodeSimpleCombinatorialHits(results_string);
    }

    @Override
    public boolean hasGUI() {
        return false;
    }
}
