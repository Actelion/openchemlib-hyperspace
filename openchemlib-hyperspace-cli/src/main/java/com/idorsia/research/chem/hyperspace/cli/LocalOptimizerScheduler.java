package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.localopt.LocalBeamOptimizer;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResult;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResultWriter;
import com.idorsia.research.chem.hyperspace.localopt.SeedAssembly;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

class LocalOptimizerScheduler {

    private final int threads;
    private final LocalBeamOptimizer optimizer;
    private final LocalOptimizationRequest request;
    private final Path outputPath;

    LocalOptimizerScheduler(int threads,
                            LocalBeamOptimizer optimizer,
                            LocalOptimizationRequest request,
                            Path outputPath) {
        this.threads = Math.max(1, threads);
        this.optimizer = optimizer;
        this.request = request;
        this.outputPath = outputPath;
    }

    void run(List<SeedAssembly> seeds) throws IOException, InterruptedException {
        if (threads == 1) {
            try (LocalOptimizationResultWriter writer = new LocalOptimizationResultWriter(outputPath)) {
                for (SeedAssembly seed : seeds) {
                    LocalOptimizationResult result = optimizer.optimize(seed, request);
                    writer.write(result);
                }
            }
            return;
        }

        ExecutorService executor = Executors.newFixedThreadPool(threads);
        CompletionService<LocalOptimizationResult> completionService = new ExecutorCompletionService<>(executor);
        List<Future<LocalOptimizationResult>> futures = new ArrayList<>();
        for (SeedAssembly seed : seeds) {
            futures.add(completionService.submit(new SeedTask(seed)));
        }

        executor.shutdown();

        try (LocalOptimizationResultWriter writer = new LocalOptimizationResultWriter(outputPath)) {
            for (int i = 0; i < futures.size(); i++) {
                Future<LocalOptimizationResult> future = completionService.take();
                try {
                    writer.write(future.get());
                } catch (Exception e) {
                    executor.shutdownNow();
                    throw new IOException("Seed optimization failed", e);
                }
            }
        } finally {
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        }
    }

    private class SeedTask implements Callable<LocalOptimizationResult> {
        private final SeedAssembly seed;

        private SeedTask(SeedAssembly seed) {
            this.seed = seed;
        }

        @Override
        public LocalOptimizationResult call() {
            return optimizer.optimize(seed, request);
        }
    }
}
