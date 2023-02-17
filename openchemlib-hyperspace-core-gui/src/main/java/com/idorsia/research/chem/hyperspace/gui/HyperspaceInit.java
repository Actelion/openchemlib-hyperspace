package com.idorsia.research.chem.hyperspace.gui;

import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.search.*;
import org.apache.commons.lang3.tuple.Pair;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class HyperspaceInit {


    /**
     *
     *
     * @param gui can be null
     * @param name
     * @param config
     * @return
     */
    public static AbstractSearchProvider resolveSearchProvider(HyperspaceSearchGUI gui, String name, String config) throws Exception {


        if(name.equals( "HyperspaceSSS" ) ) {
            HyperspaceSubstructureSearch.InitializationConfig init_config = new HyperspaceSubstructureSearch.InitializationConfig();
            init_config.deserializeFromJSON(config);
            HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSSSRemote")) {
            RemoteSearchProviderAdapter.RemoteSearchProviderConfig remote_config = new RemoteSearchProviderAdapter.RemoteSearchProviderConfig();
            remote_config.deserializeFromJSON(config);
            //HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch(gui,null);
            HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch();
            hsss.setConfigurationAndGUI(null,gui);
            RemoteSearchProviderAdapter provider = new RemoteSearchProviderAdapter(hsss,remote_config);
            return provider;
        }

        if(name.equals("HyperspaceSSSForDW")) {
            HyperspaceSubstructureSearch.InitializationConfig init_config = new HyperspaceSubstructureSearch.InitializationConfig();
            init_config.deserializeFromJSON(config);
            HyperspaceSubstructureSearchForDW hsss = new HyperspaceSubstructureSearchForDW();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSimilarity")) {
            HyperspaceSubstructureSearch.InitializationConfig init_config = new HyperspaceSubstructureSearch.InitializationConfig();
            init_config.deserializeFromJSON(config);
            HyperspaceSimilaritySearch hsss = new HyperspaceSimilaritySearch();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSimilarityRemote")) {
            RemoteSearchProviderAdapter.RemoteSearchProviderConfig remote_config = new RemoteSearchProviderAdapter.RemoteSearchProviderConfig();
            remote_config.deserializeFromJSON(config);
            //HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch(gui,null);
            HyperspaceSimilaritySearch hsss = new HyperspaceSimilaritySearch();
            hsss.setConfigurationAndGUI(null,gui);
            RemoteSearchProviderAdapter provider = new RemoteSearchProviderAdapter(hsss,remote_config);
            return provider;
        }


        throw new Exception("Could not find service with name "+name);
    }


    /**
     *
     * @param gui can be null
     * @param json_init
     * @return
     * @throws Exception
     */
    public static Map<String,AbstractSearchProvider> loadSearchProviderInitFile(HyperspaceSearchGUI gui, String json_init) throws Exception {

        Map<String,AbstractSearchProvider> search_providers = new HashMap<>();

        JSONObject jo = new JSONObject(json_init);
        JSONArray j_search_providers = jo.getJSONArray("ServiceProviders");
        Iterator<Object> sp_it = j_search_providers.iterator();
        while(sp_it.hasNext()) {
            JSONObject ji = (JSONObject) sp_it.next();

            String search_provider_name = ji.getString("ServiceProvider");
            JSONObject config           = ji.getJSONObject("Config");
            String service_name = ji.getString("ServiceName");

            AbstractSearchProvider provider = resolveSearchProvider(gui,search_provider_name,config.toString());
            search_providers.put(service_name,provider);
        }

        return search_providers;
    }


    /**
     *
     *
     * @param args
     */
    public static void main(String args[]) {

    }

}
