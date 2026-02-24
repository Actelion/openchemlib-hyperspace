package com.idorsia.research.chem.hyperspace.screening;

import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.localopt.LocalBeamOptimizer;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResult;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResultWriter;
import com.idorsia.research.chem.hyperspace.localopt.PheSAAssemblyScorer;
import com.idorsia.research.chem.hyperspace.localopt.SeedAssembly;
import com.idorsia.research.chem.hyperspace.localopt.SynthonSetAccessor;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Chains sampling, optional micro optimization, and full local optimization.
 */
public final class ContinuousScreeningOrchestrator implements AutoCloseable {

    private final CandidateSampler sampler;
    private final ReactionScheduler scheduler;
    private final DuplicateFilter duplicateFilter;
    private final ScreeningMetrics metrics = new ScreeningMetrics();
    private final LocalBeamOptimizer fullOptimizer;
    private final LocalOptimizationRequest fullRequest;
    private final ExecutorService executor;
    private final LocalOptimizationResultWriter writer;
    private final ScheduledExecutorService progressExecutor;
    private final int progressIntervalSeconds;
    private final AtomicLong jobCounter = new AtomicLong();
    private final long baseSeed;
    private final AtomicBoolean writerClosed = new AtomicBoolean(false);
    private final AtomicBoolean progressStarted = new AtomicBoolean(false);
    private volatile long targetIterations = -1;
    private volatile long startTimeNanos = 0L;
    private final boolean microEnabled;
    private final LocalBeamOptimizer microOptimizer;
    private final LocalOptimizationRequest microRequest;

    public ContinuousScreeningOrchestrator(Config config) throws IOException {
        DownsampledSynthonSpace downsampledView = config.downsampledRaw.toDownsampledSpace();
        this.sampler = new CandidateSampler(downsampledView,
                config.queryDescriptor,
                config.samplerConfig,
                metrics.candidateComparisonCounter(),
                new CandidateSampler.ReactionSamplingTracker() {
                    @Override
                    public void recordAttempt(String reactionId) {
                        metrics.recordSamplingAttempt(reactionId);
                    }

                    @Override
                    public void recordSuccess(String reactionId) {
                        metrics.recordSamplingSuccess(reactionId);
                    }
                });
        this.scheduler = new ReactionScheduler(computeReactionSizes(downsampledView),
                config.reactionMinWeight,
                config.reactionWeightExponent);
        this.duplicateFilter = new DuplicateFilter(config.duplicateCacheSize);
        SynthonSetAccessor fullAccessor = new SynthonSetAccessor(config.fullRaw);
        this.fullOptimizer = new LocalBeamOptimizer(fullAccessor,
                new PheSAAssemblyScorer(config.queryDescriptor,
                        config.fullOptimizationRequest.getMinScoreThreshold(),
                        metrics.fullOptComparisonCounter()));
        this.fullRequest = config.fullOptimizationRequest;
        this.writer = new LocalOptimizationResultWriter(config.hitOutput);
        this.progressIntervalSeconds = config.progressIntervalSeconds;
        int threads = Math.max(1, config.workerThreads);
        int queueCapacity = Math.max(1, config.queueCapacity);
        this.executor = new ThreadPoolExecutor(
                threads,
                threads,
                0L,
                TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<>(queueCapacity),
                new ThreadPoolExecutor.CallerRunsPolicy());
        this.progressExecutor = java.util.concurrent.Executors.newSingleThreadScheduledExecutor(r -> {
            Thread t = new Thread(r, "ScreeningProgress");
            t.setDaemon(true);
            return t;
        });
        this.baseSeed = config.randomSeed;
        this.microEnabled = config.microEnabled;
        if (microEnabled) {
            SynthonSetAccessor downsampledAccessor = new SynthonSetAccessor(downsampledView);
            this.microOptimizer = new LocalBeamOptimizer(downsampledAccessor,
                    new PheSAAssemblyScorer(config.queryDescriptor,
                            config.microOptimizationRequest.getMinScoreThreshold(),
                            metrics.microOptComparisonCounter()));
            this.microRequest = config.microOptimizationRequest;
        } else {
            this.microOptimizer = null;
            this.microRequest = null;
        }
    }

    private Map<String, List<Integer>> computeReactionSizes(DownsampledSynthonSpace view) {
        Map<String, List<Integer>> sizes = new HashMap<>();
        view.getDownsampledSets().forEach((rxnId, fragMap) -> {
            List<Integer> counts = new ArrayList<>();
            fragMap.values().forEach(list -> counts.add(list.size()));
            sizes.put(rxnId, counts);
        });
        return sizes;
    }

    public void run(long iterations) {
        long remaining = iterations < 0 ? Long.MAX_VALUE : iterations;
        try {
            startProgressReporter(iterations);
            while (remaining > 0) {
                remaining--;
                long jobId = jobCounter.getAndIncrement();
                submitJob(jobId);
            }
            shutdownAndAwait();
        } finally {
            stopProgressReporter();
            closeWriter();
        }
    }

    private void submitJob(long jobId) {
        try {
            executor.execute(new SampleAndOptimizeJob(jobId));
        } catch (RejectedExecutionException ignored) {
            // executor is shutting down; ignore new submissions
        }
    }

    private void shutdownAndAwait() {
        executor.shutdown();
        try {
            while (!executor.awaitTermination(1, TimeUnit.MINUTES)) {
                // keep waiting
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            executor.shutdownNow();
        }
    }

    private ScreeningCandidate maybeMicroOptimize(ScreeningCandidate candidate) {
        if (!microEnabled) {
            return candidate;
        }
        SeedAssembly seed = new SeedAssembly(candidate.getReactionId(), candidate.getFragmentIds(), candidate.getSimilarity());
        LocalOptimizationResult result = microOptimizer.optimize(seed, microRequest);
        if (result.getBeamEntries().isEmpty()) {
            return null;
        }
        LocalOptimizationResult.BeamEntry best = result.getBeamEntries().get(0);
        List<String> fragmentIds = new ArrayList<>();
        best.getFragments().forEach(frag -> fragmentIds.add(frag.fragment_id));
        metrics.incrementMicroOptimized();
        return new ScreeningCandidate(candidate.getReactionId(),
                fragmentIds,
                best.getIdcode(),
                best.getScore(),
                best.getAtomCount(),
                best.getRotatableBonds());
    }

    @Override
    public void close() {
        executor.shutdownNow();
        stopProgressReporter();
        closeWriter();
    }

    public ScreeningMetrics getMetrics() {
        return metrics;
    }

    private void closeWriter() {
        if (!writerClosed.compareAndSet(false, true)) {
            return;
        }
        synchronized (writer) {
            try {
                writer.close();
            } catch (IOException ignored) {
            }
        }
    }

    private void startProgressReporter(long iterations) {
        if (!progressStarted.compareAndSet(false, true)) {
            return;
        }
        if (progressIntervalSeconds <= 0) {
            return;
        }
        this.targetIterations = iterations;
        this.startTimeNanos = System.nanoTime();
        progressExecutor.scheduleAtFixedRate(this::reportProgress,
                progressIntervalSeconds,
                progressIntervalSeconds,
                TimeUnit.SECONDS);
    }

    private void stopProgressReporter() {
        progressExecutor.shutdownNow();
    }

    private void reportProgress() {
        try {
            long elapsedSeconds = Math.max(1L, (System.nanoTime() - startTimeNanos) / 1_000_000_000L);
            long sampled = metrics.getSampled();
            long submitted = metrics.getSubmitted();
            long hits = metrics.getHits();
            long duplicateSeeds = metrics.getDuplicateSeeds();
            long candidateComparisons = metrics.getCandidateComparisons();
            long microComparisons = metrics.getMicroOptComparisons();
            long fullComparisons = metrics.getFullOptComparisons();

            StringBuilder sb = new StringBuilder();
            sb.append("[Screening] elapsed=").append(elapsedSeconds).append("s ");
            sb.append("sampled=").append(sampled).append(" submitted=").append(submitted)
                    .append(" hits=").append(hits).append(" dupSeeds=").append(duplicateSeeds);
            if (targetIterations > 0) {
                int pct = (int) Math.min(100L, Math.round(sampled * 100.0 / targetIterations));
                sb.append(" progress=").append(pct).append("%");
            }
            sb.append(" comps=").append(candidateComparisons + microComparisons + fullComparisons)
                    .append(" (cand=").append(candidateComparisons)
                    .append(", micro=").append(microComparisons)
                    .append(", full=").append(fullComparisons).append(")");
            sb.append(" rate=").append(String.format("%.1f", sampled / (double) elapsedSeconds)).append("/s");
            System.out.println(sb);
        } catch (Exception ignored) {
        }
    }

    private final class SampleAndOptimizeJob implements Runnable {
        private final long jobId;

        private SampleAndOptimizeJob(long jobId) {
            this.jobId = jobId;
        }

        @Override
        public void run() {
            Random rng = new Random(baseSeed ^ jobId);
            metrics.incrementSampled();
            String reactionId = scheduler.pick(rng);
            ScreeningCandidate candidate = sampler.sample(reactionId, rng);
            if (candidate == null) {
                return;
            }
            if (duplicateFilter.markIfDuplicate(candidate.getAssembledIdcode())) {
                metrics.incrementDuplicateSeeds();
                return;
            }
            candidate = maybeMicroOptimize(candidate);
            if (candidate == null) {
                return;
            }
            metrics.incrementSubmitted();
            SeedAssembly seed = new SeedAssembly(candidate.getReactionId(),
                    candidate.getFragmentIds(),
                    candidate.getSimilarity());
            LocalOptimizationResult result = fullOptimizer.optimize(seed, fullRequest);
            if (!result.getBeamEntries().isEmpty()) {
                metrics.recordHit(result.getReactionId());
            }
            synchronized (writer) {
                try {
                    writer.write(result);
                } catch (IOException e) {
                    throw new RuntimeException("Unable to write optimization result", e);
                }
            }
        }
    }

    public static final class Config {
        private RawSynthonSpace fullRaw;
        private RawSynthonSpace downsampledRaw;
        private PheSAMolecule queryDescriptor;
        private CandidateSampler.Config samplerConfig;
        private boolean microEnabled;
        private LocalOptimizationRequest microOptimizationRequest;
        private LocalOptimizationRequest fullOptimizationRequest;
        private int workerThreads = 4;
        private int queueCapacity = 1000;
        private int progressIntervalSeconds = 60;
        private double reactionWeightExponent = 1.0;
        private double reactionMinWeight = 0.01;
        private Path hitOutput;
        private int duplicateCacheSize = 200_000;
        private long randomSeed = 13L;

        public Config withFullRaw(RawSynthonSpace fullRaw) {
            this.fullRaw = fullRaw;
            return this;
        }

        public Config withDownsampledRaw(RawSynthonSpace downsampledRaw) {
            this.downsampledRaw = downsampledRaw;
            return this;
        }

        public Config withQueryDescriptor(PheSAMolecule descriptor) {
            this.queryDescriptor = descriptor;
            return this;
        }

        public Config withSamplerConfig(CandidateSampler.Config config) {
            this.samplerConfig = config;
            return this;
        }

        public Config withMicroEnabled(boolean enabled) {
            this.microEnabled = enabled;
            return this;
        }

        public Config withMicroOptimizationRequest(LocalOptimizationRequest request) {
            this.microOptimizationRequest = request;
            return this;
        }

        public Config withFullOptimizationRequest(LocalOptimizationRequest request) {
            this.fullOptimizationRequest = request;
            return this;
        }

        public Config withFullOptimizerThreads(int threads) {
            this.workerThreads = threads;
            return this;
        }

        public Config withQueueCapacity(int capacity) {
            this.queueCapacity = capacity;
            return this;
        }

        public Config withProgressIntervalSeconds(int seconds) {
            this.progressIntervalSeconds = seconds;
            return this;
        }

        public Config withReactionWeightExponent(double exponent) {
            this.reactionWeightExponent = exponent;
            return this;
        }

        public Config withReactionMinWeight(double minWeight) {
            this.reactionMinWeight = minWeight;
            return this;
        }

        public Config withHitOutput(Path output) {
            this.hitOutput = output;
            return this;
        }

        public Config withDuplicateCacheSize(int size) {
            this.duplicateCacheSize = size;
            return this;
        }

        public Config withRandomSeed(long seed) {
            this.randomSeed = seed;
            return this;
        }

        public void validate() {
            if (fullRaw == null || downsampledRaw == null) {
                throw new IllegalArgumentException("Both full and downsampled raw spaces are required");
            }
            if (queryDescriptor == null) {
                throw new IllegalArgumentException("Query descriptor missing");
            }
            if (samplerConfig == null) {
                throw new IllegalArgumentException("Sampler config missing");
            }
            if (fullOptimizationRequest == null) {
                throw new IllegalArgumentException("Full optimization request missing");
            }
            if (hitOutput == null) {
                throw new IllegalArgumentException("Hit output path missing");
            }
            if (microEnabled && microOptimizationRequest == null) {
                throw new IllegalArgumentException("Micro optimization request missing");
            }
            if (workerThreads <= 0) {
                throw new IllegalArgumentException("Worker threads must be positive");
            }
            if (queueCapacity <= 0) {
                throw new IllegalArgumentException("Queue capacity must be positive");
            }
            if (progressIntervalSeconds < 0) {
                throw new IllegalArgumentException("Progress interval must be >= 0");
            }
            if (reactionWeightExponent <= 0) {
                throw new IllegalArgumentException("Reaction weight exponent must be positive");
            }
            if (reactionMinWeight < 0) {
                throw new IllegalArgumentException("Reaction min weight must be >= 0");
            }
        }
    }
}
