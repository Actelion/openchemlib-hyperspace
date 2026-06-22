package com.idorsia.research.chem.hyperspace.cli;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationLogLevel;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.screening.CandidateSampler;
import com.idorsia.research.chem.hyperspace.screening.ReactionScheduler;

import java.io.IOException;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

public final class ContinuousScreeningConfig {

    private InputPaths inputs = new InputPaths();
    private QueryInput query = new QueryInput();
    private String queryFile;
    private CandidateSampling sampling = new CandidateSampling();
    private FullOptimization fullOptimization = new FullOptimization();
    private MicroOptimization microOptimization = new MicroOptimization();
    private Orchestration orchestration = new Orchestration();
    private RunControl run = new RunControl();
    private OutputPaths output = new OutputPaths();

    public static ContinuousScreeningConfig fromJson(Path path) throws IOException {
        ObjectMapper mapper = new ObjectMapper()
                .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, true);
        ContinuousScreeningConfig config = mapper.readValue(path.toFile(), ContinuousScreeningConfig.class);
        config.applyDefaults();
        config.validate();
        return config;
    }

    private void applyDefaults() {
        if (inputs == null) {
            inputs = new InputPaths();
        }
        if (query == null) {
            query = new QueryInput();
        }
        if (sampling == null) {
            sampling = new CandidateSampling();
        }
        if (fullOptimization == null) {
            fullOptimization = new FullOptimization();
        }
        if (fullOptimization.request == null) {
            fullOptimization.request = OptimizationRequestSettings.fullDefaults();
        }
        if (microOptimization == null) {
            microOptimization = new MicroOptimization();
        }
        if (microOptimization.request == null) {
            microOptimization.request = OptimizationRequestSettings.microDefaults();
        }
        if (orchestration == null) {
            orchestration = new Orchestration();
        }
        if (run == null) {
            run = new RunControl();
        }
        if (output == null) {
            output = new OutputPaths();
        }
    }

    private void validate() {
        requireNonBlank(inputs.rawFull, "inputs.rawFull");
        requireNonBlank(inputs.rawDownsampled, "inputs.rawDownsampled");
        requireNonBlank(output.hitsTsv, "output.hitsTsv");
        requireFinite(output.minReportedSimilarity, "output.minReportedSimilarity");

        boolean hasInlineSmiles = hasText(query.smiles);
        boolean hasInlineIdcode = hasText(query.idcode);
        boolean hasInlineSdfFile = hasText(query.sdfFile);
        boolean hasInlinePhesaFile = hasText(query.phesaFile);
        int inlineSourceCount = (hasInlineSmiles ? 1 : 0)
                + (hasInlineIdcode ? 1 : 0)
                + (hasInlineSdfFile ? 1 : 0)
                + (hasInlinePhesaFile ? 1 : 0);
        boolean hasInlineQuery = inlineSourceCount > 0;
        boolean hasExternalQuery = hasText(queryFile);
        if (hasInlineQuery == hasExternalQuery) {
            throw new IllegalArgumentException("Specify exactly one query source: inline query or queryFile");
        }
        if (hasInlineQuery && inlineSourceCount != 1) {
            throw new IllegalArgumentException("Specify exactly one of query.smiles, query.idcode, query.sdfFile, or query.phesaFile");
        }
        if (hasInlineSdfFile && query.sdfRecordIndex < 0) {
            throw new IllegalArgumentException("query.sdfRecordIndex must be >= 0");
        }

        if (run.iterations != null && run.iterations < -1) {
            throw new IllegalArgumentException("run.iterations must be >= -1");
        }
        if (hasText(run.maxRuntime) && run.maxRuntimeSeconds != null) {
            throw new IllegalArgumentException("Specify only one of run.maxRuntime or run.maxRuntimeSeconds");
        }
        run.resolveMaxRuntime();

        if (sampling.attemptsPerReaction <= 0) {
            throw new IllegalArgumentException("sampling.attemptsPerReaction must be positive");
        }
        if (sampling.minAtoms <= 0) {
            throw new IllegalArgumentException("sampling.minAtoms must be positive");
        }
        if (sampling.maxAtoms < 0) {
            throw new IllegalArgumentException("sampling.maxAtoms must be >= 0 (0 means unlimited)");
        }
        if (sampling.maxRotatableBonds <= 0) {
            throw new IllegalArgumentException("sampling.maxRotatableBonds must be positive");
        }
        requireFinite(sampling.minSimilarity, "sampling.minSimilarity");

        fullOptimization.request.validate("fullOptimization.request");

        if (microOptimization.enabled) {
            microOptimization.request.validate("microOptimization.request");
        }

        if (orchestration.queueCapacity <= 0) {
            throw new IllegalArgumentException("orchestration.queueCapacity must be positive");
        }
        if (orchestration.workerThreads <= 0) {
            throw new IllegalArgumentException("orchestration.workerThreads must be positive");
        }
        if (orchestration.progressIntervalSeconds < 0) {
            throw new IllegalArgumentException("orchestration.progressIntervalSeconds must be >= 0");
        }
        if (orchestration.duplicateCacheSize <= 0) {
            throw new IllegalArgumentException("orchestration.duplicateCacheSize must be positive");
        }
        if (orchestration.reactionWeightExponent <= 0.0) {
            throw new IllegalArgumentException("orchestration.reactionWeightExponent must be positive");
        }
        if (orchestration.reactionMinWeight < 0.0) {
            throw new IllegalArgumentException("orchestration.reactionMinWeight must be >= 0");
        }
        if (orchestration.reactionWeighting != null) {
            orchestration.reactionWeighting.validate("orchestration.reactionWeighting");
        }
    }

    public Path resolveRawFull(Path configPath) {
        return resolveRelativeToConfig(inputs.rawFull, configPath);
    }

    public Path resolveRawDownsampled(Path configPath) {
        return resolveRelativeToConfig(inputs.rawDownsampled, configPath);
    }

    public Path resolveOutputHits(Path configPath) {
        return resolveRelativeToConfig(output.hitsTsv, configPath);
    }

    public double getMinReportedSimilarity() {
        return output.minReportedSimilarity;
    }

    public QueryInput resolveQuery(Path configPath) throws IOException {
        if (hasText(queryFile)) {
            Path resolvedQueryPath = resolveRelativeToConfig(queryFile, configPath);
            ObjectMapper mapper = new ObjectMapper()
                    .configure(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, true);
            QueryInput loadedQuery = mapper.readValue(resolvedQueryPath.toFile(), QueryInput.class);
            QueryInput resolvedQuery = resolveQueryPaths(loadedQuery, resolvedQueryPath);
            validateResolvedQuery(resolvedQuery, "queryFile");
            return resolvedQuery;
        }
        QueryInput resolvedQuery = resolveQueryPaths(query, configPath);
        validateResolvedQuery(resolvedQuery, "query");
        return resolvedQuery;
    }

    public CandidateSampler.Config toSamplerConfig() {
        return new CandidateSampler.Config(
                sampling.attemptsPerReaction,
                sampling.minAtoms,
                sampling.maxAtoms,
                sampling.maxRotatableBonds,
                sampling.minSimilarity
        );
    }

    public LocalOptimizationRequest toFullOptimizationRequest() {
        return fullOptimization.request.toLocalOptimizationRequest(run.getEffectiveRandomSeed());
    }

    public LocalOptimizationRequest toMicroOptimizationRequest() {
        return microOptimization.request.toLocalOptimizationRequest(run.getEffectiveRandomSeed());
    }

    public ReactionScheduler.Weighting toReactionWeighting() {
        if (orchestration.reactionWeighting == null) {
            return null;
        }
        return orchestration.reactionWeighting.toSchedulerWeighting();
    }

    private static Path resolveRelativeToConfig(String value, Path configPath) {
        Path p = Path.of(value);
        if (p.isAbsolute()) {
            return p.normalize();
        }
        Path parent = configPath.toAbsolutePath().getParent();
        if (parent == null) {
            return p.toAbsolutePath().normalize();
        }
        return parent.resolve(p).normalize();
    }

    private static void requireFinite(double value, String name) {
        if (!Double.isFinite(value)) {
            throw new IllegalArgumentException(name + " must be finite");
        }
    }

    private static void requireNonBlank(String value, String name) {
        if (value == null || value.isBlank()) {
            throw new IllegalArgumentException(name + " is required");
        }
    }

    private static boolean hasText(String value) {
        return value != null && !value.isBlank();
    }

    private static void validateResolvedQuery(QueryInput resolvedQuery, String sourceName) {
        if (resolvedQuery == null) {
            throw new IllegalArgumentException(sourceName + " did not contain a query object");
        }
        boolean hasSmiles = hasText(resolvedQuery.smiles);
        boolean hasIdcode = hasText(resolvedQuery.idcode);
        boolean hasSdfFile = hasText(resolvedQuery.sdfFile);
        boolean hasPhesaFile = hasText(resolvedQuery.phesaFile);
        int sourceCount = (hasSmiles ? 1 : 0)
                + (hasIdcode ? 1 : 0)
                + (hasSdfFile ? 1 : 0)
                + (hasPhesaFile ? 1 : 0);
        if (sourceCount != 1) {
            throw new IllegalArgumentException(sourceName + " must specify exactly one of smiles, idcode, sdfFile, or phesaFile");
        }
        if (hasSdfFile && resolvedQuery.sdfRecordIndex < 0) {
            throw new IllegalArgumentException(sourceName + ".sdfRecordIndex must be >= 0");
        }
    }

    private static QueryInput resolveQueryPaths(QueryInput sourceQuery, Path sourcePath) {
        if (sourceQuery == null) {
            return null;
        }
        QueryInput resolved = new QueryInput();
        resolved.setSmiles(sourceQuery.smiles);
        resolved.setIdcode(sourceQuery.idcode);
        resolved.setSdfRecordIndex(sourceQuery.sdfRecordIndex);
        if (hasText(sourceQuery.sdfFile)) {
            resolved.setSdfFile(resolveRelativeToConfig(sourceQuery.sdfFile, sourcePath).toString());
        }
        if (hasText(sourceQuery.phesaFile)) {
            resolved.setPhesaFile(resolveRelativeToConfig(sourceQuery.phesaFile, sourcePath).toString());
        }
        return resolved;
    }

    public static final class InputPaths {
        private String rawFull;
        private String rawDownsampled;

        public String getRawFull() {
            return rawFull;
        }

        public void setRawFull(String rawFull) {
            this.rawFull = rawFull;
        }

        public String getRawDownsampled() {
            return rawDownsampled;
        }

        public void setRawDownsampled(String rawDownsampled) {
            this.rawDownsampled = rawDownsampled;
        }
    }

    public static final class QueryInput {
        private String smiles;
        private String idcode;
        private String sdfFile;
        private String phesaFile;
        private int sdfRecordIndex = 0;

        public String getSmiles() {
            return smiles;
        }

        public void setSmiles(String smiles) {
            this.smiles = smiles;
        }

        public String getIdcode() {
            return idcode;
        }

        public void setIdcode(String idcode) {
            this.idcode = idcode;
        }

        public String getSdfFile() {
            return sdfFile;
        }

        public void setSdfFile(String sdfFile) {
            this.sdfFile = sdfFile;
        }

        public String getPhesaFile() {
            return phesaFile;
        }

        public void setPhesaFile(String phesaFile) {
            this.phesaFile = phesaFile;
        }

        public int getSdfRecordIndex() {
            return sdfRecordIndex;
        }

        public void setSdfRecordIndex(int sdfRecordIndex) {
            this.sdfRecordIndex = sdfRecordIndex;
        }
    }

    public static final class CandidateSampling {
        private long attemptsPerReaction = 20;
        private int minAtoms = 1;
        private int maxAtoms = 0;
        private int maxRotatableBonds = 15;
        private double minSimilarity = 0.5;

        public long getAttemptsPerReaction() {
            return attemptsPerReaction;
        }

        public void setAttemptsPerReaction(long attemptsPerReaction) {
            this.attemptsPerReaction = attemptsPerReaction;
        }

        public int getMinAtoms() {
            return minAtoms;
        }

        public void setMinAtoms(int minAtoms) {
            this.minAtoms = minAtoms;
        }

        public int getMaxAtoms() {
            return maxAtoms;
        }

        public void setMaxAtoms(int maxAtoms) {
            this.maxAtoms = maxAtoms;
        }

        public int getMaxRotatableBonds() {
            return maxRotatableBonds;
        }

        public void setMaxRotatableBonds(int maxRotatableBonds) {
            this.maxRotatableBonds = maxRotatableBonds;
        }

        public double getMinSimilarity() {
            return minSimilarity;
        }

        public void setMinSimilarity(double minSimilarity) {
            this.minSimilarity = minSimilarity;
        }
    }

    public static final class FullOptimization {
        private OptimizationRequestSettings request = OptimizationRequestSettings.fullDefaults();

        public OptimizationRequestSettings getRequest() {
            return request;
        }

        public void setRequest(OptimizationRequestSettings request) {
            this.request = request;
        }
    }

    public static final class MicroOptimization {
        private boolean enabled = false;
        private OptimizationRequestSettings request = OptimizationRequestSettings.microDefaults();

        public boolean isEnabled() {
            return enabled;
        }

        public void setEnabled(boolean enabled) {
            this.enabled = enabled;
        }

        public OptimizationRequestSettings getRequest() {
            return request;
        }

        public void setRequest(OptimizationRequestSettings request) {
            this.request = request;
        }
    }

    public static final class Orchestration {
        private int workerThreads = 4;
        private int queueCapacity = 1000;
        private int progressIntervalSeconds = 20;
        private double reactionWeightExponent = 1.0;
        private double reactionMinWeight = 0.01;
        private ReactionWeightingSettings reactionWeighting;
        private int duplicateCacheSize = 200_000;

        public int getWorkerThreads() {
            return workerThreads;
        }

        public void setWorkerThreads(int workerThreads) {
            this.workerThreads = workerThreads;
        }

        public int getQueueCapacity() {
            return queueCapacity;
        }

        public void setQueueCapacity(int queueCapacity) {
            this.queueCapacity = queueCapacity;
        }

        public int getProgressIntervalSeconds() {
            return progressIntervalSeconds;
        }

        public void setProgressIntervalSeconds(int progressIntervalSeconds) {
            this.progressIntervalSeconds = progressIntervalSeconds;
        }

        public double getReactionWeightExponent() {
            return reactionWeightExponent;
        }

        public void setReactionWeightExponent(double reactionWeightExponent) {
            this.reactionWeightExponent = reactionWeightExponent;
        }

        public double getReactionMinWeight() {
            return reactionMinWeight;
        }

        public void setReactionMinWeight(double reactionMinWeight) {
            this.reactionMinWeight = reactionMinWeight;
        }

        public ReactionWeightingSettings getReactionWeighting() {
            return reactionWeighting;
        }

        public void setReactionWeighting(ReactionWeightingSettings reactionWeighting) {
            this.reactionWeighting = reactionWeighting;
        }

        public int getDuplicateCacheSize() {
            return duplicateCacheSize;
        }

        public void setDuplicateCacheSize(int duplicateCacheSize) {
            this.duplicateCacheSize = duplicateCacheSize;
        }
    }

    public static final class ReactionWeightingSettings {
        private String mode = ReactionScheduler.Mode.EXPONENT.name();
        private List<ReactionWeightBucketSettings> buckets = new ArrayList<>();

        private void validate(String prefix) {
            ReactionScheduler.Mode parsedMode = parseMode(prefix + ".mode");
            if (parsedMode == ReactionScheduler.Mode.EXPONENT) {
                return;
            }
            toSchedulerWeighting();
        }

        private ReactionScheduler.Weighting toSchedulerWeighting() {
            ReactionScheduler.Mode parsedMode = parseMode("reactionWeighting.mode");
            if (parsedMode == ReactionScheduler.Mode.EXPONENT) {
                return null;
            }
            List<ReactionScheduler.Bucket> schedulerBuckets = new ArrayList<>();
            if (buckets != null) {
                for (ReactionWeightBucketSettings bucket : buckets) {
                    if (bucket == null) {
                        schedulerBuckets.add(null);
                    } else {
                        schedulerBuckets.add(new ReactionScheduler.Bucket(bucket.maxProductExclusive, bucket.weight));
                    }
                }
            }
            return ReactionScheduler.Weighting.bucketedProduct(schedulerBuckets);
        }

        private ReactionScheduler.Mode parseMode(String key) {
            if (mode == null || mode.isBlank()) {
                throw new IllegalArgumentException(key + " must not be blank");
            }
            try {
                return ReactionScheduler.Mode.valueOf(mode.toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                String options = Arrays.stream(ReactionScheduler.Mode.values())
                        .map(Enum::name)
                        .collect(Collectors.joining(", "));
                throw new IllegalArgumentException(key + " must be one of: " + options);
            }
        }

        public String getMode() {
            return mode;
        }

        public void setMode(String mode) {
            this.mode = mode;
        }

        public List<ReactionWeightBucketSettings> getBuckets() {
            return buckets;
        }

        public void setBuckets(List<ReactionWeightBucketSettings> buckets) {
            this.buckets = buckets;
        }
    }

    public static final class ReactionWeightBucketSettings {
        private Double maxProductExclusive;
        private double weight;

        public Double getMaxProductExclusive() {
            return maxProductExclusive;
        }

        public void setMaxProductExclusive(Double maxProductExclusive) {
            this.maxProductExclusive = maxProductExclusive;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }
    }

    public static final class RunControl {
        private Long iterations = null;
        private String maxRuntime;
        private Long maxRuntimeSeconds;
        private Long randomSeed = 13L;
        private Long effectiveRandomSeed;

        public long getIterations() {
            return iterations == null ? -1L : iterations;
        }

        public void setIterations(Long iterations) {
            this.iterations = iterations;
        }

        public String getMaxRuntime() {
            return maxRuntime;
        }

        public void setMaxRuntime(String maxRuntime) {
            this.maxRuntime = maxRuntime;
        }

        public Long getMaxRuntimeSeconds() {
            return maxRuntimeSeconds;
        }

        public void setMaxRuntimeSeconds(Long maxRuntimeSeconds) {
            this.maxRuntimeSeconds = maxRuntimeSeconds;
        }

        public Long getRandomSeed() {
            return randomSeed;
        }

        public void setRandomSeed(Long randomSeed) {
            this.randomSeed = randomSeed;
        }

        public long getEffectiveRandomSeed() {
            if (effectiveRandomSeed == null) {
                effectiveRandomSeed = randomSeed != null ? randomSeed : createTimeBasedSeed();
            }
            return effectiveRandomSeed;
        }

        public Duration resolveMaxRuntime() {
            if (maxRuntimeSeconds != null) {
                if (maxRuntimeSeconds <= 0L) {
                    throw new IllegalArgumentException("run.maxRuntimeSeconds must be positive");
                }
                return Duration.ofSeconds(maxRuntimeSeconds);
            }
            if (!hasText(maxRuntime)) {
                return null;
            }
            return parseDuration(maxRuntime);
        }

        private static Duration parseDuration(String value) {
            String trimmed = value.trim();
            if (trimmed.isEmpty()) {
                throw new IllegalArgumentException("run.maxRuntime must not be blank");
            }
            try {
                Duration duration = Duration.parse(trimmed);
                if (duration.isZero() || duration.isNegative()) {
                    throw new IllegalArgumentException("run.maxRuntime must be positive");
                }
                return duration;
            } catch (Exception ignored) {
                // Also accept compact values like 24h, 90m, 3600s, or plain seconds.
            }

            String lower = trimmed.toLowerCase(Locale.ROOT);
            long multiplier = 1L;
            String numberPart = lower;
            if (lower.endsWith("ms")) {
                multiplier = -1L;
                numberPart = lower.substring(0, lower.length() - 2);
            } else if (lower.endsWith("s")) {
                numberPart = lower.substring(0, lower.length() - 1);
            } else if (lower.endsWith("m")) {
                multiplier = 60L;
                numberPart = lower.substring(0, lower.length() - 1);
            } else if (lower.endsWith("h")) {
                multiplier = 3600L;
                numberPart = lower.substring(0, lower.length() - 1);
            } else if (lower.endsWith("d")) {
                multiplier = 86_400L;
                numberPart = lower.substring(0, lower.length() - 1);
            }

            try {
                long amount = Long.parseLong(numberPart.trim());
                Duration duration = multiplier < 0L
                        ? Duration.ofMillis(amount)
                        : Duration.ofSeconds(Math.multiplyExact(amount, multiplier));
                if (duration.isZero() || duration.isNegative()) {
                    throw new IllegalArgumentException("run.maxRuntime must be positive");
                }
                return duration;
            } catch (RuntimeException e) {
                throw new IllegalArgumentException("run.maxRuntime must be an ISO-8601 duration or a value like 24h, 90m, 3600s", e);
            }
        }

        private static long createTimeBasedSeed() {
            long now = System.nanoTime();
            long wall = System.currentTimeMillis();
            long thread = Thread.currentThread().getId();
            return now ^ Long.rotateLeft(wall, 21) ^ Long.rotateLeft(thread, 42);
        }
    }

    public static final class OutputPaths {
        private String hitsTsv;
        private double minReportedSimilarity = 0.0;

        public String getHitsTsv() {
            return hitsTsv;
        }

        public void setHitsTsv(String hitsTsv) {
            this.hitsTsv = hitsTsv;
        }

        public double getMinReportedSimilarity() {
            return minReportedSimilarity;
        }

        public void setMinReportedSimilarity(double minReportedSimilarity) {
            this.minReportedSimilarity = minReportedSimilarity;
        }
    }

    public static final class OptimizationRequestSettings {
        private int beamSize = 10;
        private int neighborPoolSize = 8;
        private int sampledNeighbors = 4;
        private int perPositionCap = 2;
        private int maxRounds = 6;
        private int patience = 3;
        private double minPhesaSimilarity = 0.55;
        private Double minScoreThreshold = 0.0;
        private double improvementTolerance = 1e-4;
        private boolean reportAllCandidates = true;
        private Long randomSeed = null;
        private String logLevel = LocalOptimizationLogLevel.NONE.name();

        static OptimizationRequestSettings fullDefaults() {
            return new OptimizationRequestSettings();
        }

        static OptimizationRequestSettings microDefaults() {
            OptimizationRequestSettings request = new OptimizationRequestSettings();
            request.beamSize = 5;
            request.neighborPoolSize = 4;
            request.sampledNeighbors = 4;
            request.perPositionCap = 2;
            request.maxRounds = 2;
            request.patience = 1;
            request.minPhesaSimilarity = 0.0;
            request.minScoreThreshold = 0.0;
            request.improvementTolerance = 1e-4;
            request.reportAllCandidates = true;
            request.randomSeed = null;
            request.logLevel = LocalOptimizationLogLevel.NONE.name();
            return request;
        }

        public void validate(String namePrefix) {
            if (beamSize <= 0) {
                throw new IllegalArgumentException(namePrefix + ".beamSize must be positive");
            }
            if (neighborPoolSize <= 0) {
                throw new IllegalArgumentException(namePrefix + ".neighborPoolSize must be positive");
            }
            if (sampledNeighbors <= 0) {
                throw new IllegalArgumentException(namePrefix + ".sampledNeighbors must be positive");
            }
            if (perPositionCap <= 0) {
                throw new IllegalArgumentException(namePrefix + ".perPositionCap must be positive");
            }
            if (maxRounds <= 0) {
                throw new IllegalArgumentException(namePrefix + ".maxRounds must be positive");
            }
            if (patience < 0) {
                throw new IllegalArgumentException(namePrefix + ".patience must be >= 0");
            }
            requireFinite(minPhesaSimilarity, namePrefix + ".minPhesaSimilarity");
            if (minScoreThreshold != null) {
                requireFinite(minScoreThreshold, namePrefix + ".minScoreThreshold");
            }
            requireFinite(improvementTolerance, namePrefix + ".improvementTolerance");
            parseLogLevel(logLevel, namePrefix + ".logLevel");
        }

        public LocalOptimizationRequest toLocalOptimizationRequest(long defaultSeed) {
            long effectiveSeed = randomSeed != null ? randomSeed : defaultSeed;
            LocalOptimizationRequest.Builder builder = LocalOptimizationRequest.builder()
                    .beamSize(beamSize)
                    .neighborPoolSize(neighborPoolSize)
                    .sampledNeighbors(sampledNeighbors)
                    .perPositionCap(perPositionCap)
                    .maxRounds(maxRounds)
                    .patience(patience)
                    .minPhesaSimilarity(minPhesaSimilarity)
                    .improvementTolerance(improvementTolerance)
                    .reportAllCandidates(reportAllCandidates)
                    .randomSeed(effectiveSeed)
                    .logLevel(parseLogLevel(logLevel, "request.logLevel"));
            if (minScoreThreshold != null) {
                builder.minScoreThreshold(minScoreThreshold);
            }
            return builder.build();
        }

        private static LocalOptimizationLogLevel parseLogLevel(String value, String key) {
            if (value == null || value.isBlank()) {
                throw new IllegalArgumentException(key + " must not be blank");
            }
            try {
                return LocalOptimizationLogLevel.valueOf(value.toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                String options = Arrays.stream(LocalOptimizationLogLevel.values())
                        .map(Enum::name)
                        .collect(Collectors.joining(", "));
                throw new IllegalArgumentException(key + " must be one of: " + options);
            }
        }

        public int getBeamSize() {
            return beamSize;
        }

        public void setBeamSize(int beamSize) {
            this.beamSize = beamSize;
        }

        public int getNeighborPoolSize() {
            return neighborPoolSize;
        }

        public void setNeighborPoolSize(int neighborPoolSize) {
            this.neighborPoolSize = neighborPoolSize;
        }

        public int getSampledNeighbors() {
            return sampledNeighbors;
        }

        public void setSampledNeighbors(int sampledNeighbors) {
            this.sampledNeighbors = sampledNeighbors;
        }

        public int getPerPositionCap() {
            return perPositionCap;
        }

        public void setPerPositionCap(int perPositionCap) {
            this.perPositionCap = perPositionCap;
        }

        public int getMaxRounds() {
            return maxRounds;
        }

        public void setMaxRounds(int maxRounds) {
            this.maxRounds = maxRounds;
        }

        public int getPatience() {
            return patience;
        }

        public void setPatience(int patience) {
            this.patience = patience;
        }

        public double getMinPhesaSimilarity() {
            return minPhesaSimilarity;
        }

        public void setMinPhesaSimilarity(double minPhesaSimilarity) {
            this.minPhesaSimilarity = minPhesaSimilarity;
        }

        public Double getMinScoreThreshold() {
            return minScoreThreshold;
        }

        public void setMinScoreThreshold(Double minScoreThreshold) {
            this.minScoreThreshold = minScoreThreshold;
        }

        public double getImprovementTolerance() {
            return improvementTolerance;
        }

        public void setImprovementTolerance(double improvementTolerance) {
            this.improvementTolerance = improvementTolerance;
        }

        public boolean isReportAllCandidates() {
            return reportAllCandidates;
        }

        public void setReportAllCandidates(boolean reportAllCandidates) {
            this.reportAllCandidates = reportAllCandidates;
        }

        public Long getRandomSeed() {
            return randomSeed;
        }

        public void setRandomSeed(Long randomSeed) {
            this.randomSeed = randomSeed;
        }

        public String getLogLevel() {
            return logLevel;
        }

        public void setLogLevel(String logLevel) {
            this.logLevel = logLevel;
        }
    }

    public InputPaths getInputs() {
        return inputs;
    }

    public void setInputs(InputPaths inputs) {
        this.inputs = inputs;
    }

    public QueryInput getQuery() {
        return query;
    }

    public void setQuery(QueryInput query) {
        this.query = query;
    }

    public String getQueryFile() {
        return queryFile;
    }

    public void setQueryFile(String queryFile) {
        this.queryFile = queryFile;
    }

    public CandidateSampling getSampling() {
        return sampling;
    }

    public void setSampling(CandidateSampling sampling) {
        this.sampling = sampling;
    }

    public FullOptimization getFullOptimization() {
        return fullOptimization;
    }

    public void setFullOptimization(FullOptimization fullOptimization) {
        this.fullOptimization = fullOptimization;
    }

    public MicroOptimization getMicroOptimization() {
        return microOptimization;
    }

    public void setMicroOptimization(MicroOptimization microOptimization) {
        this.microOptimization = microOptimization;
    }

    public Orchestration getOrchestration() {
        return orchestration;
    }

    public void setOrchestration(Orchestration orchestration) {
        this.orchestration = orchestration;
    }

    public RunControl getRun() {
        return run;
    }

    public void setRun(RunControl run) {
        this.run = run;
    }

    public OutputPaths getOutput() {
        return output;
    }

    public void setOutput(OutputPaths output) {
        this.output = output;
    }
}
