package com.idorsia.research.chem.hyperspace.screening;

import com.idorsia.research.chem.hyperspace.localopt.LocalBeamOptimizer;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResult;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResultWriter;
import com.idorsia.research.chem.hyperspace.localopt.PheSAAssemblyScorer;
import com.idorsia.research.chem.hyperspace.localopt.SeedAssembly;
import com.idorsia.research.chem.hyperspace.localopt.SynthonSetAccessor;
import com.actelion.research.chem.phesa.PheSAMolecule;

import java.io.IOException;
import java.nio.file.Path;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Handles asynchronous execution of the full local optimization stage.
 */
public final class FullOptimizerManager implements AutoCloseable {

    private final ExecutorService executor;
    private final CompletionService<OptimizationOutcome> completionService;
    private final AtomicBoolean accepting = new AtomicBoolean(true);
    private final AtomicInteger inflight = new AtomicInteger();
    private final LocalBeamOptimizer optimizer;
    private final LocalOptimizationRequest request;
    private final LocalOptimizationResultWriter writer;
    private final ScreeningMetrics metrics;
    private final Thread drainerThread;

    public FullOptimizerManager(SynthonSetAccessor accessor,
                                PheSAMolecule queryDescriptor,
                                LocalOptimizationRequest request,
                                int threads,
                                Path outputPath,
                                ScreeningMetrics metrics) throws IOException {
        this.optimizer = new LocalBeamOptimizer(accessor, new PheSAAssemblyScorer(queryDescriptor, request.getMinScoreThreshold()));
        this.request = request;
        this.metrics = metrics;
        this.writer = new LocalOptimizationResultWriter(outputPath);
        this.executor = Executors.newFixedThreadPool(Math.max(1, threads));
        this.completionService = new ExecutorCompletionService<>(executor);
        this.drainerThread = new Thread(this::drainLoop, "FullOptimizerDrainer");
        this.drainerThread.start();
    }

    public void submitCandidate(ScreeningCandidate candidate) {
        if (!accepting.get()) {
            return;
        }
        metrics.incrementSubmitted();
        inflight.incrementAndGet();
        completionService.submit(new OptimizationJob(candidate));
    }

    private void drainLoop() {
        try {
            while (accepting.get() || inflight.get() > 0) {
                OptimizationOutcome outcome;
                try {
                    outcome = completionService.take().get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    break;
                } catch (ExecutionException e) {
                    inflight.decrementAndGet();
                    e.getCause().printStackTrace();
                    continue;
                }
                inflight.decrementAndGet();
                if (outcome == null) {
                    continue;
                }
                try {
                    writer.write(outcome.result());
                } catch (IOException e) {
                    throw new RuntimeException("Unable to write optimization result", e);
                }
                if (!outcome.result().getBeamEntries().isEmpty()) {
                    metrics.recordHit(outcome.result().getReactionId());
                }
            }
        } finally {
            try {
                writer.close();
            } catch (IOException ignored) {
            }
        }
    }

    @Override
    public void close() {
        accepting.set(false);
        executor.shutdownNow();
        drainerThread.interrupt();
        try {
            drainerThread.join();
        } catch (InterruptedException ignored) {
        }
    }

    private final class OptimizationJob implements Callable<OptimizationOutcome> {
        private final ScreeningCandidate candidate;

        private OptimizationJob(ScreeningCandidate candidate) {
            this.candidate = candidate;
        }

        @Override
        public OptimizationOutcome call() {
            SeedAssembly seed = new SeedAssembly(candidate.getReactionId(),
                    candidate.getFragmentIds(),
                    candidate.getSimilarity());
            LocalOptimizationResult result = optimizer.optimize(seed, request);
            return new OptimizationOutcome(candidate, result);
        }
    }

    private record OptimizationOutcome(ScreeningCandidate candidate,
                                       LocalOptimizationResult result) {
    }
}
