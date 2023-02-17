package com.idorsia.research.chem.hyperspace.gui;


import com.idorsia.research.chem.hyperspace.HyperspaceEngine;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceSearchProcess;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearchProcess;
import com.idorsia.research.chem.hyperspace.gui.search.RemoteSearchProviderAdapter;
import org.apache.commons.lang3.tuple.Pair;
import org.json.JSONObject;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.web.embedded.tomcat.TomcatServletWebServerFactory;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.RestController;
//import org.openmolecules.comm.ServerCommunicator;
//import org.openmolecules.comm.ServerTask;
//import org.openmolecules.comm.ServerTaskFactory;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;


/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-Server
 */



/**
 * NOTE: if the main function is called with a single parameter, then
 *       it assumes, that it is the path to a input file that contains
 *       the actual comand line parameters, i.e. the input file
 *       is parsed, and the main function is called with the
 *       parsed arguments.
 *
 *
 */
@SpringBootApplication
@RestController
public class HyperspaceServer {

    private static final int DEFAULT_PORT = 8090;
    private static final int DEFAULT_MAX_QUERIES_COUNT = 4;
    private static final int DEFAULT_THREADS_PER_QUERY = 4;

    private static int maxQueries      = 4;
    private static int port            = 8090;
    public  static int threadsPerQuery = 8;


    static Map<String,AbstractSearchProvider> searchProviders = new HashMap<>();

    //private static ArrayList<String> mSpacesSim          = new ArrayList<String>();
    //private static ArrayList<String> mSpacesSS           = new ArrayList<String>();

    //private static HyperspaceEngine mEngine              = new HyperspaceEngine();

    public static void main(String args[]) {
        System.out.println("Hyperspace Server 2.1.1 (Java HTTP version)");
        System.out.println("Time: "+getDateAndTime());

        HyperspaceEngine.setLogLevelToQuiet();

        if (!parseArguments(args)) {
            showUsage();
        }

        initServer();
        System.out.println("Server is up and running!");
    }

    private static void initServer() {
        System.out.println("Start initializing Analysis Providers..");
        // start init processes:
        for(String service_name : searchProviders.keySet()) {
            System.out.println("Init "+service_name+"\n");
            AbstractHyperspaceProcess proc_init = searchProviders.get(service_name).startInitialization();
            while(searchProviders.get(service_name).getStatus() != AbstractSearchProvider.SearchProviderStatus.READY &&
                    searchProviders.get(service_name).getStatus() != AbstractSearchProvider.SearchProviderStatus.ERROR) {
                double init_status = searchProviders.get(service_name).getStatusOfInitialization();
                System.out.println( String.format("initializing: %.2f",init_status));
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            System.out.println("Done.. Service is: "+searchProviders.get(service_name).getStatus());
        }

        System.out.println("Done initializing all providers!");

        // set port and max threads:
        SpringApplication app = new SpringApplication(HyperspaceServer.class);
        app.setDefaultProperties(Collections.singletonMap("server.port", port));
        String args2[] = new String[0];
        app.run(args2);

        // TODO: implement maxQueries for default tomcat server..

        //TomcatServletWebServerFactory tomcatFactory = (TomcatServletWebServerFactory) app.getWebServer().getFactory();
        //tomcatFactory.addConnectorCustomizers(connector -> {
        //    connector.setMaxThreads(200);
        //});

        // OLD CODE FROM SIMPLE SERVER:
//        ServerTaskFactory factory = new ServerTaskFactory() {
//            public ServerTask createServerTask() {
//                return new HyperspaceTask();
//            }
//        };
//
//        try {
//            ServerCommunicator.initialize(factory, maxQueries, port);
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
    }


    private static boolean parseArguments(String args[]) {
        HyperspaceServer.port = DEFAULT_PORT;
        int maxQueries = DEFAULT_MAX_QUERIES_COUNT;
        int TpQ = DEFAULT_THREADS_PER_QUERY;
        String settingsFile = null;
        //boolean processMMP = false;
        boolean verbose = false;

        for (int i=0; i<args.length; i++) {
            if (args[i].equals("-f") && args.length > i + 1) {
                settingsFile = args[++i];
                continue;
            }
            if (args[i].equals("-p") && args.length > i + 1) {
                try {
                    HyperspaceServer.port = Integer.parseInt(args[++i]);
                    continue;
                } catch (NumberFormatException e) {
                    return false;
                }
            }
            if (args[i].equals("-s") && args.length > i + 1) {
                try {
                    HyperspaceServer.maxQueries = Integer.parseInt(args[++i]);
                    continue;
                } catch (NumberFormatException e) {
                    return false;
                }
            }
            if (args[i].equals("-t") && args.length > i + 1) {
                try {
                    TpQ = Integer.parseInt(args[++i]);
                    HyperspaceServer.threadsPerQuery = TpQ;
                    continue;
                } catch (NumberFormatException e) {
                    return false;
                }
            }
            if (args[i].equals("--init") && args.length > i + 1 ) {
                String init_file = args[++i];
                String init_file_txt = HyperspaceUtils.readFileIntoString(init_file);
                try{
                    Map<String,AbstractSearchProvider> providers = HyperspaceInit.loadSearchProviderInitFile(null,init_file_txt);
                    searchProviders.putAll(providers);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

//            if (args[i].equals("--sss") && args.length > i + 1) {
//                try {
//                    String config_file_path = args[++i];
//                    mEngine.loadSpace_SS(config_file_path);
//                    continue;
//                } catch (NumberFormatException e) {
//                    return false;
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
//            if (args[i].equals("--sim") && args.length > i + 1) {
//                try {
//                    String config_file_path = args[++i];
//                    mEngine.loadSpace_Similarity(config_file_path);
//                    continue;
//                } catch (NumberFormatException e) {
//                    return false;
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
        }
        return true;
    }

    private static void showUsage() {
        System.out.println("Usage from console:");
        System.out.println("hyperspace_server [-p port] [-s maxQueries] [-t tpq] --init <init_file>");
        System.out.println("  -p  Default port is "+DEFAULT_PORT+". Use option -p to choose a different port.");
        System.out.println("  -s  Maximum number of simultaneously handled queries (default:"+DEFAULT_MAX_QUERIES_COUNT+").");
        System.out.println("  -t  Threads per handled queries (default:"+DEFAULT_THREADS_PER_QUERY+").");
    }

    private static String getDateAndTime() {
        return new SimpleDateFormat("d-MMM-yyyy HH:mm:ss").format(Calendar.getInstance().getTime());
    }


    public static String createFullJSONQuery(RemoteSearchProviderAdapter.RemoteSearchProviderConfig remote_config, AbstractSearchProvider.SearchConfiguration config) {
        JSONObject full_query = new JSONObject();
        full_query.put("service",remote_config.searchServiceName);
        full_query.put("query",config.serializeToJSON());
        return full_query.toString();
    }


    /**
     *
     * @return (service,query)
     */
    public static Pair<String,JSONObject> extractFullJSONQuery(String json_full_query) {
        JSONObject full_query = new JSONObject(json_full_query);
        String      service_str  = full_query.getString("service");
        JSONObject  query_str    = full_query.getJSONObject("query");
        return Pair.of(service_str,query_str);
    }

    public static JSONObject createFullJSONQuery(String service, JSONObject query) {
        JSONObject full_query = new JSONObject();
        full_query.put("service",service);
        full_query.put("query",query);
        return full_query;
    }

    @GetMapping("/search")
    public ResponseEntity<String> call(@RequestParam(value = "data") String val){
        System.out.println("Request-data= "+val);
        String what = new String( Base64.getUrlDecoder().decode(val) );


    //// TODO: different handling of different subtasks, e.g. AbstractSubstructuer Search, etc.. (i.e. introduce reasonable intermediate abstractions..)
    //public static class HyperspaceTask extends ServerTask {
    //    @Override
    //    public void performTask() {
    //        String what_encoded = getRequestText("data");
    //
    //        if (what_encoded == null) {
    //            createErrorResponse("Request missing");
    //            return;
    //        }
    //
    //        String what = new String( Base64.getUrlDecoder().decode(what_encoded) );
            //HyperspaceSubstructureSearch.SubstructureSearchConfiguration search_config = new HyperspaceSubstructureSearch.SubstructureSearchConfiguration();
            AbstractSearchProvider.SearchConfiguration search_config = null;

            Pair<String,JSONObject> full_query = extractFullJSONQuery(what);

            // get service and process..
            if(!searchProviders.containsKey(full_query.getLeft())) {
                //createErrorResponse("Service not found: "+full_query.getLeft());
                return new ResponseEntity<>("Service not found: "+full_query.getLeft(),HttpStatus.NOT_FOUND);
                //return;
            }

            try {
                // reflection
                search_config = (AbstractSearchProvider.SearchConfiguration) searchProviders.get(full_query.getLeft()).getSearchConfigurationClass().getConstructor().newInstance();
                search_config.deserializeFromJSON(full_query.getRight().toString());
                //search_config.deserializeFromJSON(full_query.getRight().toString());
            } catch (Exception e) {
                e.printStackTrace();
                //createErrorResponse("Problem decoding request..");
                return new ResponseEntity<>("Problem decoding request..",HttpStatus.BAD_REQUEST);
                //return;
            }


            AbstractHyperspaceProcess proc_search = null;
            try {
                search_config.setNumberOfThreads( threadsPerQuery );
                proc_search = searchProviders.get(full_query.getLeft()).runSearch(search_config);
            } catch (AbstractSearchProvider.IncompatibleSearchConfigurationException e) {
                e.printStackTrace();
                return new ResponseEntity<>("Incompatible search configuration..",HttpStatus.BAD_REQUEST);
                //createErrorResponse("Incompatible search configuration..");
                //return;
            }

            boolean done_or_failed = false;
            try {
                done_or_failed = proc_search.waitUntilDoneOrFailed(3600000);
            } catch (InterruptedException e) {
                e.printStackTrace();
                return new ResponseEntity<>("Interrupted exception during computation..",HttpStatus.INTERNAL_SERVER_ERROR);
                //createErrorResponse("Interrupted exception during computation..");
                //return;
            }

            if(!done_or_failed) {
                System.out.println("[HyperspaceTask::performTask] :: timeout..");
                return new ResponseEntity<>("Interrupted exception during computation..",HttpStatus.INTERNAL_SERVER_ERROR);
                //createErrorResponse("Reached timeout during computation..");
                //return;
            }

            if(proc_search.getProcessStatus() == AbstractHyperspaceProcess.ProcessStatus.DONE) {
                // all good, return results
                //List<SynthonSpace.CombinatorialHit> search_results = ((HyperspaceSubstructureSearchProcess)proc_search).getSearchResults();
                List search_results = ((AbstractHyperspaceSearchProcess)proc_search).getSearchResults();

                try {
                    //createPlainResponse( searchProviders.get(full_query.getLeft()).encodeResultsToString(search_results) );
                    String result = searchProviders.get(full_query.getLeft()).encodeResultsToString(search_results);
                    return new ResponseEntity<>( result , HttpStatus.OK);
                } catch (IOException ex) {
                    System.out.println("[HyperspaceTask::performTask] :: IOException while encoding results");
                    return new ResponseEntity<>("IOException while encoding results..",HttpStatus.INTERNAL_SERVER_ERROR);
                    //createErrorResponse("IOException while encoding results..");
                    //ex.printStackTrace();
                    //return;
                }
            }

            // We should not end here
            return null;
        }
}
