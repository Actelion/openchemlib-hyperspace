package com.idorsia.research.hyperspace.server;


import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RestController;

import java.util.List;
import java.util.stream.Collectors;

@RestController
public class ServiceProviderController {

    private final ConfigurationService configurationService;

    public ServiceProviderController(ConfigurationService configurationService) {
        this.configurationService = configurationService;
    }

    @GetMapping("/serviceProviders")
    public List<String> getAvailableServiceProviders() {
        return configurationService.getServiceProviders().stream()
                .map(ServiceConfig.ServiceProvider::getServiceName)
                .collect(Collectors.toList());
    }

}
