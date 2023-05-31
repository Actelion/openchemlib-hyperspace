package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceServer;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.util.EntityUtils;
import org.json.JSONObject;

import javax.swing.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Base64;
import java.util.List;

public class RemoteSearchProviderAdapter<T> extends AbstractSearchProvider<T> {

    private final AbstractSearchProvider<T> searchProvider;
    private RemoteSearchProviderConfig remoteConfig;

    public static class RemoteSearchProviderConfig implements SearchProviderConfiguration{
        public String server;
        public int    port;
        public  String searchServiceName;
        public String guiServiceName;

        public RemoteSearchProviderConfig(){}
        public RemoteSearchProviderConfig(String server, int port, String searchServiceName) {
            this.server = server;
            this.port = port;
            this.searchServiceName = searchServiceName;
        }
        public String serializeToJSON() {
            JSONObject jo = new JSONObject();
            jo.put("Server",this.server);
            jo.put("Port",this.port);
            jo.put("SearchServiceName",this.searchServiceName);
            jo.put("GUIServiceName",this.guiServiceName);
            return jo.toString();
        }
        public void deserializeFromJSON(String json) {
            JSONObject jo = new JSONObject(json);
            this.server = jo.getString("Server");
            this.port   = jo.getInt("Port");
            this.searchServiceName = jo.getString("SearchServiceName");
            this.guiServiceName = jo.getString("GUIServiceName");
        }

        private String serviceName;
        @Override
        public String getServiceName() {
            return this.serviceName;
        }

        @Override
        public void setServiceName(String serviceName) {
            this.serviceName = serviceName;
        }

    }

    private String guiServiceName = "NA";
    public void setGUIServiceName(String serviceName) {
        this.guiServiceName = serviceName;
    }

    public RemoteSearchProviderAdapter (AbstractSearchProvider<T> provider, RemoteSearchProviderConfig r_conf){
        this.searchProvider = provider;
        this.remoteConfig   = r_conf;

        this.setGUIServiceName(r_conf.guiServiceName);
        provider.setRemoteSearchProvider(this);
    }

    @Override
    public String getSearchName() {
        return this.guiServiceName;
        //return this.searchProvider.getSearchName();
    }

    @Override
    public String getSpaceName() {
        return this.guiServiceName;
        //return this.getSearchName();
    }

    @Override
    public int getMaxNumberOfThreads() {
        return this.searchProvider.getMaxNumberOfThreads();
    }

    @Override
    public void setConfigurationAndGUI(SearchProviderConfiguration config, HyperspaceSearchGUI gui) {

    }

    @Override
    public JPanel getSearchConfigurationView() {
        return this.searchProvider.getSearchConfigurationView();
    }

    @Override
    public SearchConfiguration getSearchConfiguration() {
        return this.searchProvider.getSearchConfiguration();
    }

    @Override
    public Class getSearchConfigurationClass() {
        return RemoteSearchProviderConfig.class;
    }

    @Override
    public SearchProviderConfiguration getSearchProviderConfiguration() {
        // todo implement..
        return this.remoteConfig;
    }

    @Override
    public AbstractHyperspaceProcess startInitialization() {
        return new AbstractHyperspaceProcess() {
            @Override
            public String getName() {
                return "Init Remote";
            }
            @Override
            public ProcessStatus getProcessStatus() {
                return ProcessStatus.DONE;
            }
        };
    }

    @Override
    public AbstractHyperspaceSearchProcess<T> runSearch(SearchConfiguration configuration) throws IncompatibleSearchConfigurationException {

        RemoteSearchProcess<T> search_process = new RemoteSearchProcess<>(this,"",this.searchProvider,configuration,remoteConfig);
        search_process.startSearchAsync();
        return search_process;
    }

    @Override
    public void presentSearchResults(List<T> results) {
        this.searchProvider.presentSearchResults(results);
    }

    @Override
    public SearchProviderStatus getStatus() {
        return SearchProviderStatus.REMOTE;
    }

    @Override
    public boolean hasGUI() {
        return true;
    }

    @Override
    public double getStatusOfInitialization() {
        return 1.0;
    }

    @Override
    public String encodeResultsToString(List<T> results) throws IOException {
        return searchProvider.encodeResultsToString(results);
    }

    @Override
    public List<T> decodeResultsFromString(String results_string) throws IOException {
        return searchProvider.decodeResultsFromString(results_string);
    }


    @Override
    public RemoteSearchProviderAdapter<T> getRemoteSearchProvider() {
        System.out.println("[WARN] :: ::getRemoteSearchProvider :: does not do anything..");
        return null;
    }

    @Override
    public void setRemoteSearchProvider(RemoteSearchProviderAdapter<T> remote_adapter) {
        System.out.println("[WARN] :: ::setRemoteSearchProvider :: does not do anything..");
    }

    public static class RemoteSearchProcess<T> extends AbstractHyperspaceSearchProcess<T> {
        String name;
        SearchConfiguration config;
        RemoteSearchProviderConfig  remoteConfig;
        private final AbstractSearchProvider<T> searchProvider;

        public RemoteSearchProcess(AbstractSearchProvider<T> remote_search_provider, String name,AbstractSearchProvider<T> searchProvider,  SearchConfiguration config, RemoteSearchProviderConfig remoteConfig) {
            super(remote_search_provider);
            this.name = name;
            this.searchProvider = searchProvider;
            this.config = config;
            this.remoteConfig = remoteConfig;

            this.setProcessStatus(ProcessStatus.WAITING);
        }

        @Override
        public void setSearchConfiguration(SearchConfiguration config) {
            this.config = config;
        }

        @Override
        public SearchConfiguration getSearchConfiguration() {
            return this.config;
        }

        @Override
        public void startSearchAsync() {
            super.startSearchAsync();
            this.setProcessStatus(ProcessStatus.COMPUTING);

            Runnable run_search = new Runnable() {
                @Override
                public void run() {
                    int timeout_sec = 3600;
                    RequestConfig config_request = RequestConfig.custom()
                            .setConnectTimeout(timeout_sec * 1000)
                            .setConnectionRequestTimeout(timeout_sec * 1000)
                            .setSocketTimeout(timeout_sec * 1000).build();

                    HttpClient httpClient = HttpClientBuilder.create().setDefaultRequestConfig(config_request).build();
                    try {

                        //StringEntity params = new StringEntity("details={\"name\":\"xyz\",\"age\":\"20\"} ");
                        // create full config..
                        JSONObject full_query = HyperspaceServer.createFullJSONQuery(remoteConfig.searchServiceName, new JSONObject(config.serializeToJSON()) );

                        String json_data = full_query.toString();//config.serializeToJSON();

                        //String base64_json_data = Base64.getEncoder().encodeToString(json_data.getBytes());
                        String base64_json_data = Base64.getUrlEncoder().encodeToString(json_data.getBytes());
                        //request.addHeader("content-type", "application/xml");
                        //request.set
                        //HttpGet request = new HttpGet("http://"+remoteConfig.server+":"+remoteConfig.port+"?data="+base64_json_data);
                        HttpGet request = new HttpGet("http://"+remoteConfig.server+":"+remoteConfig.port+"/search?data="+base64_json_data);

                        setProcessStatusMessage("Waiting for Result..");
                        HttpResponse response  = httpClient.execute(request);
                        HttpEntity resp_entity =  response.getEntity();
                        if (resp_entity != null) {
                            // return it as a String
                            String result = EntityUtils.toString(resp_entity,"UTF-8");
                            // for some reason this adds a newline to the end, we have to remove this.
                            //result = result.substring(0,result.length()-1);
                            //System.out.println(result);

                            //byte[] result_data = Base64.getDecoder().decode(result);
                            List<T> decoded_search_results = searchProvider.decodeResultsFromString(result);
                            search_results = decoded_search_results;

                            setProcessStatusMessage("Done");
                            setProcessStatus(ProcessStatus.DONE);
                        }
                        else {
                            System.out.println("[RemoteSearchProviderAdapter::run_search] :: null");
                            setProcessStatusMessage("HttpResponse is null");
                            setProcessStatus(ProcessStatus.FAILED);
                        }
                    } catch (Exception ex) {
                        ex.printStackTrace();
                        setProcessStatusMessage("Exception");
                        setProcessStatus(ProcessStatus.FAILED);
                    } finally {
                        // @Deprecated httpClient.getConnectionManager().shutdown();
                    }

                }
            };

            Thread ri = new Thread(run_search);
            ri.start();
        }

        List<T> search_results = null;

        //public static List<T>

        @Override
        public List<T> getSearchResults() {
            return search_results;
        }

        @Override
        public String getName() {
            return this.name;
        }
    }

}
