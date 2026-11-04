package com.idorsia.research.chem.hyperspace.gui;

import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.search.*;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class HyperspaceInit {
    private static final Map<AbstractSearchProvider, String> providerTypes =
            java.util.Collections.synchronizedMap(new java.util.WeakHashMap<>());

    public static Map<String, AbstractSearchProvider> loadSearchProviderInitFile(HyperspaceSearchGUI gui, File file) throws Exception {
        if (!file.isFile()) throw new IOException("Configuration file not found: " + file.getAbsolutePath());
        JSONObject json = new JSONObject(java.nio.file.Files.readString(file.toPath()));
        JSONArray providers = json.getJSONArray("ServiceProviders");
        for (int i = 0; i < providers.length(); i++) {
            JSONObject provider = providers.getJSONObject(i);
            String type = provider.getString("ServiceProvider");
            if (type.equals("HyperspaceSSS") || type.equals("HyperspaceSSSForDW") || type.equals("HyperspaceSimilarity")) {
                JSONObject config = provider.getJSONObject("Config");
                java.nio.file.Path path = file.toPath().toAbsolutePath().getParent().resolve(config.getString("File")).normalize();
                if (!java.nio.file.Files.isRegularFile(path) || !java.nio.file.Files.isReadable(path))
                    throw new IOException("Space file is missing or unreadable: " + path);
                config.put("File", path.toString());
            }
        }
        return loadSearchProviderInitFile(gui, json.toString());
    }



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
            init_config.setServiceName(name);
            HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSSSRemote")) {
            RemoteSearchProviderAdapter.RemoteSearchProviderConfig remote_config = new RemoteSearchProviderAdapter.RemoteSearchProviderConfig();
            remote_config.deserializeFromJSON(config);
            remote_config.setServiceName(name);
            //HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch(gui,null);
            HyperspaceSubstructureSearch hsss = new HyperspaceSubstructureSearch();
            hsss.setConfigurationAndGUI(null,gui);
            RemoteSearchProviderAdapter provider = new RemoteSearchProviderAdapter(hsss,remote_config);
            return provider;
        }

        if(name.equals("HyperspaceSSSForDW")) {
            HyperspaceSubstructureSearch.InitializationConfig init_config = new HyperspaceSubstructureSearch.InitializationConfig();
            init_config.deserializeFromJSON(config);
            init_config.setServiceName(name);
            HyperspaceSubstructureSearchForDW hsss = new HyperspaceSubstructureSearchForDW();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSimilarity")) {
            HyperspaceSubstructureSearch.InitializationConfig init_config = new HyperspaceSubstructureSearch.InitializationConfig();
            init_config.deserializeFromJSON(config);
            init_config.setServiceName(name);
            HyperspaceSimilaritySearch hsss = new HyperspaceSimilaritySearch();
            hsss.setConfigurationAndGUI(init_config,gui);
            return hsss;
        }

        if(name.equals("HyperspaceSimilarityRemote")) {
            RemoteSearchProviderAdapter.RemoteSearchProviderConfig remote_config = new RemoteSearchProviderAdapter.RemoteSearchProviderConfig();
            remote_config.deserializeFromJSON(config);
            remote_config.setServiceName(name);
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
            provider.getSearchProviderConfiguration().setServiceName(service_name);
            providerTypes.put(provider, search_provider_name);
            search_providers.put(service_name,provider);
        }

        return search_providers;
    }


    public static void saveSearchProviderInitFile(HyperspaceSearchGUI gui, File output) {
        List<AbstractSearchProvider> sps = gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getSearchProviderListPanel().getSearchProviders();

        saveSearchProviders(sps, output);
    }

    public static void saveSearchProviders(List<AbstractSearchProvider> providers, File output) {
        JSONArray array = new JSONArray();
        for (AbstractSearchProvider provider : providers) {
            String type = providerTypes.get(provider);
            if (type == null) {
                if (provider instanceof HyperspaceSimilaritySearch) type = "HyperspaceSimilarity";
                else if (provider instanceof HyperspaceSubstructureSearchForDW) type = "HyperspaceSSSForDW";
                else if (provider instanceof HyperspaceSubstructureSearch) type = "HyperspaceSSS";
                else throw new IllegalArgumentException("Cannot serialize provider: " + provider.getClass().getName());
            }
            JSONObject config = new JSONObject(provider.getSearchProviderConfiguration().serializeToJSON());
            if (config.has("File")) config.put("File", new File(config.getString("File")).getAbsolutePath());
            array.put(new JSONObject().put("ServiceName", provider.getSearchProviderConfiguration().getServiceName())
                    .put("ServiceProvider", type).put("Config", config));
        }
        try {
            java.nio.file.Files.writeString(output.toPath(), new JSONObject().put("ServiceProviders", array).toString(2));
        } catch (IOException e) { throw new java.io.UncheckedIOException(e); }
    }

    /**
     *
     *
     * @param args
     */
    public static void main(String args[]) {

    }

}
