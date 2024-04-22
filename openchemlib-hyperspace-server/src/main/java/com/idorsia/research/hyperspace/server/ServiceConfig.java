package com.idorsia.research.hyperspace.server;


import com.fasterxml.jackson.annotation.JsonProperty;

import java.util.List;

public class ServiceConfig {

    @JsonProperty("ServiceProviders")
    private List<ServiceProvider> serviceProviders;

    public List<ServiceProvider> getServiceProviders() {
        return serviceProviders;
    }

    public void setServiceProviders(List<ServiceProvider> serviceProviders) {
        this.serviceProviders = serviceProviders;
    }

    public static class ServiceProvider {
        @JsonProperty("ServiceName")
        private String serviceName;

        @JsonProperty("Config")
        private ServiceDetails config;

        public String getServiceName() {
            return serviceName;
        }

        public void setServiceName(String serviceName) {
            this.serviceName = serviceName;
        }

        public ServiceDetails getConfig() {
            return config;
        }

        public void setConfig(ServiceDetails config) {
            this.config = config;
        }
    }

    public static class ServiceDetails {
        @JsonProperty("SpaceName")
        private String spaceName;

        @JsonProperty("File")
        private String file;

        public String getSpaceName() {
            return spaceName;
        }

        public void setSpaceName(String spaceName) {
            this.spaceName = spaceName;
        }

        public String getFile() {
            return file;
        }

        public void setFile(String file) {
            this.file = file;
        }
    }
}
