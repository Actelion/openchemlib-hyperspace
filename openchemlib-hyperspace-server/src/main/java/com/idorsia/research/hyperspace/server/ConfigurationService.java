package com.idorsia.research.hyperspace.server;


import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import org.apache.commons.io.input.CountingInputStream;
import org.springframework.core.io.Resource;
import org.springframework.core.io.ResourceLoader;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.GZIPInputStream;

@Service
public class ConfigurationService {

    private final ResourceLoader resourceLoader;
    private List<ServiceConfig.ServiceProvider> serviceProviders = Collections.emptyList();

    private Map<ServiceConfig.ServiceProvider, SynthonSpace> loadedSpaces = new ConcurrentHashMap<>();

    public ConfigurationService(ResourceLoader resourceLoader) {
        this.resourceLoader = resourceLoader;
    }

    @PostConstruct
    public void loadConfiguration() throws IOException {
        Resource resource = resourceLoader.getResource("classpath:config.json");
        String json = new String(Files.readAllBytes(Paths.get(resource.getURI())));

        ObjectMapper objectMapper = new ObjectMapper();
        ServiceConfig serviceConfig = objectMapper.readValue(json, ServiceConfig.class);
        this.serviceProviders = new ArrayList<>();
        // TODO: initialize service providers..
        for(ServiceConfig.ServiceProvider sp : serviceConfig.getServiceProviders()) {
            loadServiceProvider(sp);
        }
    }

    public void loadServiceProvider(ServiceConfig.ServiceProvider providerConf) {
        Runnable ri = new Runnable() {
            @Override
            public void run() {
                try {
                    File fi_file = new File(providerConf.getConfig().getFile());
                    long file_size_in_bytes = fi_file.length();

                    FileInputStream file_in = new FileInputStream(fi_file);
                    CountingInputStream counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
                    ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(counting_in_stream));

                    System.out.println("Start loading space from file: "+fi_file);
                    SynthonSpace space_a = null;
                    space_a = (SynthonSpace) in.readObject();
                    space_a.initAfterJavaDeserialization();
                    //space = space_a;

                    System.out.println("Loaded space: " + space_a.getSpaceInfoString());

                    loadedSpaces.put(providerConf,space_a);
                } catch (FileNotFoundException e) {
                    throw new RuntimeException(e);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                } catch (ClassNotFoundException e) {
                    throw new RuntimeException(e);
                }
            }
        };
        Thread ti = new Thread(ri);
        ti.start();
    }

    public List<ServiceConfig.ServiceProvider> getServiceProviders() {
        return serviceProviders;
    }

    public Map<ServiceConfig.ServiceProvider, SynthonSpace> getLoadedSpaces() {
        return this.loadedSpaces;
    }

}
