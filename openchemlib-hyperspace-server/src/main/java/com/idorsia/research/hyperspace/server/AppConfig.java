package com.idorsia.research.hyperspace.server;


import org.apache.tomcat.util.threads.ThreadPoolExecutor;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.scheduling.concurrent.ThreadPoolTaskExecutor;

@Configuration
public class AppConfig {

    @Bean(name = "hyperspaceServer")
    public ThreadPoolTaskExecutor computationExecutor() {
        ThreadPoolTaskExecutor executor = new ThreadPoolTaskExecutor();
        executor.setCorePoolSize(8); // Number of threads to keep in the pool
        executor.setMaxPoolSize(8); // Maximum allowed number of threads
        executor.setQueueCapacity(100); // Size of the queue to hold tasks before they are executed
        executor.setThreadNamePrefix("Computation-");
        executor.initialize();
        //executor.setRejectedExecutionHandler(new ThreadPoolExecutor.AbortPolicy());
        return executor;
    }
}
