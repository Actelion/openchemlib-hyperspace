package com.idorsia.research.chem.hyperspace.screening;

import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.localopt.LocalBeamOptimizer;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationResult;
import com.idorsia.research.chem.hyperspace.localopt.PheSAAssemblyScorer;
import com.idorsia.research.chem.hyperspace.localopt.SeedAssembly;
import com.idorsia.research.chem.hyperspace.localopt.SynthonSetAccessor;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Chains sampling, optional micro optimization, and full local optimization.
 */
public final class ContinuousScreeningOrchestrator implements AutoCloseable {

    private final CandidateSampler sampler;
    private final ReactionScheduler scheduler;
    private final DuplicateFilter duplicateFilter;
    private final ScreeningMetrics metrics = new ScreeningMetrics();
    private final FullOptimizerManager fullManager;
    private final Random random;
    private final boolean microEnabled;
    private final LocalBeamOptimizer microOptimizer;
    private final LocalOptimizationRequest microRequest;

    public ContinuousScreeningOrchestrator(Config config) throws IOException {
        DownsampledSynthonSpace downsampledView = config.downsampledRaw.toDownsampledSpace();
        this.sampler = new CandidateSampler(downsampledView, config.queryDescriptor, config.samplerConfig);
        this.scheduler = new ReactionScheduler(computeReactionSizes(downsampledView));
        this.duplicateFilter = new DuplicateFilter(config.duplicateCacheSize);
        SynthonSetAccessor fullAccessor = new SynthonSetAccessor(config.fullRaw);
        this.fullManager = new FullOptimizerManager(fullAccessor,
                config.queryDescriptor,
                config.fullOptimizationRequest,
                config.fullOptimizerThreads,
                config.hitOutput,
                metrics);
        this.random = new Random(config.randomSeed);
        this.microEnabled = config.microEnabled;
        if (microEnabled) {
            SynthonSetAccessor downsampledAccessor = new SynthonSetAccessor(downsampledView);
            this.microOptimizer = new LocalBeamOptimizer(downsampledAccessor,
                    new PheSAAssemblyScorer(config.queryDescriptor, config.microOptimizationRequest.getMinScoreThreshold()));
            this.microRequest = config.microOptimizationRequest;
        } else {
            this.microOptimizer = null;
            this.microRequest = null;
        }
    }

    private Map<String, Integer> computeReactionSizes(DownsampledSynthonSpace view) {
        Map<String, Integer> sizes = new HashMap<>();
        view.getDownsampledSets().forEach((rxnId, fragMap) -> {
            int count = fragMap.values().stream().mapToInt(List::size).sum();
            sizes.put(rxnId, count);
        });
        return sizes;
    }

    public void run(long iterations) {
        long remaining = iterations < 0 ? Long.MAX_VALUE : iterations;
        while (remaining > 0) {
            remaining--;
            metrics.incrementSampled();
            String reactionId = scheduler.pick(random);
            ScreeningCandidate candidate = sampler.sample(reactionId, random);
            if (candidate == null) {
                continue;
            }
            if (duplicateFilter.markIfDuplicate(candidate.getAssembledIdcode())) {
                metrics.incrementDuplicates();
                continue;
            }
            candidate = maybeMicroOptimize(candidate);
            if (candidate == null) {
                continue;
            }
            if (duplicateFilter.markIfDuplicate(candidate.getAssembledIdcode())) {
                metrics.incrementDuplicates();
                continue;
            }
            fullManager.submitCandidate(candidate);
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
        List<String> fragmentIds = best.getFragments().stream()
                .map(frag -> frag.fragment_id)
                .toList();
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
        fullManager.close();
    }

    public ScreeningMetrics getMetrics() {
        return metrics;
    }

    public static final class Config {
        private RawSynthonSpace fullRaw;
        private RawSynthonSpace downsampledRaw;
        private PheSAMolecule queryDescriptor;
        private CandidateSampler.Config samplerConfig;
        private boolean microEnabled;
        private LocalOptimizationRequest microOptimizationRequest;
        private LocalOptimizationRequest fullOptimizationRequest;
        private int fullOptimizerThreads = 4;
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
            this.fullOptimizerThreads = threads;
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
        }
    }
}
