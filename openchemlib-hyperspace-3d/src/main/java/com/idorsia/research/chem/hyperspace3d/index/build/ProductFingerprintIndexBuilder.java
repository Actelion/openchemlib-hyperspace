package com.idorsia.research.chem.hyperspace3d.index.build;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.FeaturizationResult;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexShard;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexWriter;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintRecord;
import com.idorsia.research.chem.hyperspace3d.model.ProductEmbeddingBatchEncoder;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/** Builds a persistent, sampled, query-independent product embedding index. */
public final class ProductFingerprintIndexBuilder {
    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT)
            .enable(SerializationFeature.ORDER_MAP_ENTRIES_BY_KEYS);
    private static final String STATE_FILE = "build-state.json";
    private static final String PENDING_STATE_FILE = "build-state.pending.json";
    private static final String MANIFEST_FILE = "manifest.json";
    private static final ThreadLocal<IDCodeParser> PARSER = ThreadLocal.withInitial(IDCodeParser::new);

    private final RawSynthonSpace full;
    private final RawSynthonSpace downsampled;
    private final ProductEmbeddingBatchEncoder encoder;
    private final ProductFingerprintIndexBuildConfig config;
    private final Provenance provenance;
    private final OCLDeepSpaceFeaturizer featurizer = new OCLDeepSpaceFeaturizer();
    private final DeepSpaceTensorBatchBuilder tensorBuilder = new DeepSpaceTensorBatchBuilder(featurizer);

    public ProductFingerprintIndexBuilder(RawSynthonSpace full, RawSynthonSpace downsampled,
            ProductEmbeddingBatchEncoder encoder, ProductFingerprintIndexBuildConfig config,
            Provenance provenance) {
        this.full = full;
        this.downsampled = downsampled;
        this.encoder = encoder;
        this.config = config;
        this.provenance = provenance;
        config.validate();
    }

    public ProductFingerprintIndexBuildResult build(Path outputDirectory) {
        long started = System.nanoTime();
        try {
            Files.createDirectories(outputDirectory);
            Path manifestPath = outputDirectory.resolve(MANIFEST_FILE);
            if (Files.isRegularFile(manifestPath)) return completedResult(manifestPath);

            ProductTupleSampler sampler = new ProductTupleSampler(full, downsampled,
                    config.sampling.seed, config.sampling.reactionWeighting.toWeighting());
            if (config.sampling.minimumAcceptedPerReaction > 0
                    && config.output.targetRecordCount < (long) sampler.reactionIds().size()
                            * config.sampling.minimumAcceptedPerReaction) {
                throw new ProductFingerprintIndexBuildException(
                        "targetRecordCount is smaller than the configured reaction coverage floor");
            }
            BuildState state = loadOrCreateState(outputDirectory, sampler);
            boolean resumed = state.accepted > 0;
            sampler.restore(state.sampler);
            verifyCommittedShards(outputDirectory, state.shards);
            long committedRecords = state.shards.stream()
                    .mapToLong(ProductFingerprintIndexShard::recordCount).sum();
            if (committedRecords != state.accepted || state.metrics.accepted != state.accepted) {
                throw new ProductFingerprintIndexBuildException(
                        "checkpoint record counts are inconsistent");
            }
            int queueSize = Math.max(config.runtime.cpuWorkers,
                    config.runtime.encoderBatchSize * config.runtime.queueCapacity);
            ExecutorService workers = new ThreadPoolExecutor(config.runtime.cpuWorkers,
                    config.runtime.cpuWorkers, 0L, TimeUnit.MILLISECONDS,
                    new ArrayBlockingQueue<>(queueSize),
                    new ThreadPoolExecutor.CallerRunsPolicy());
            try {
                buildRemaining(outputDirectory, sampler, state, workers, started);
            } finally {
                workers.shutdownNow();
            }
            ProductFingerprintIndexManifest manifest = createManifest(state, started);
            writeJsonAtomic(outputDirectory.resolve(MANIFEST_FILE + ".partial"), manifestPath, manifest);
            state.status = "COMPLETE";
            writeStateAtomic(outputDirectory, state);
            return new ProductFingerprintIndexBuildResult(manifestPath, state.accepted,
                    state.shards.size(), resumed, state.metrics);
        } catch (ProductFingerprintIndexBuildException error) {
            throw error;
        } catch (IOException | RuntimeException error) {
            throw new ProductFingerprintIndexBuildException("product-index build failed", error);
        }
    }

    private void buildRemaining(Path outputDirectory, ProductTupleSampler sampler,
            BuildState state, ExecutorService workers, long started) throws IOException {
        long maxProposals;
        try {
            maxProposals = Math.addExact(Math.multiplyExact(config.output.targetRecordCount,
                    config.sampling.maxAttemptsMultiplier),
                    Math.multiplyExact((long) sampler.reactionIds().size(),
                            config.sampling.coverageAttemptsPerReaction));
        } catch (ArithmeticException overflow) {
            maxProposals = Long.MAX_VALUE;
        }
        while (state.accepted < config.output.targetRecordCount) {
            int shardIndex = state.shards.size();
            long needed = Math.min(config.output.recordsPerShard,
                    config.output.targetRecordCount - state.accepted);
            Path partial = outputDirectory.resolve(shardName(shardIndex) + ".partial");
            Files.deleteIfExists(partial);
            long recordsThisShard = 0;
            try (ProductFingerprintIndexWriter writer = new ProductFingerprintIndexWriter(partial)) {
                while (recordsThisShard < needed) {
                    int capacity = (int) Math.min(config.runtime.encoderBatchSize,
                            needed - recordsThisShard);
                    boolean wasCoverage = state.coverageIndex < sampler.reactionIds().size();
                    List<Prepared> prepared;
                    if (wasCoverage) {
                        prepared = prepareCoverageBatch(sampler, state, capacity);
                    } else {
                        prepared = prepareWeightedBatch(sampler, state, workers, capacity);
                    }
                    if (state.metrics.proposals > maxProposals) {
                        throw new ProductFingerprintIndexBuildException(
                                "maximum proposal count reached before target; build remains resumable");
                    }
                    List<Prepared> accepted = account(prepared, state);
                    if (accepted.isEmpty()) {
                        if (prepared.isEmpty()) {
                            if (wasCoverage) continue;
                            throw new ProductFingerprintIndexBuildException(
                                    "all sampled product tuples are exhausted before reaching target");
                        }
                        continue;
                    }
                    long packStarted = System.nanoTime();
                    var batch = tensorBuilder.buildFromFeatures(
                            accepted.stream().map(Prepared::features).toList());
                    state.metrics.tensorPackingNanos += System.nanoTime() - packStarted;
                    long encodeStarted = System.nanoTime();
                    float[][] embeddings = encoder.encode(batch);
                    state.metrics.encodingNanos += System.nanoTime() - encodeStarted;
                    if (embeddings.length != accepted.size()) {
                        throw new ProductFingerprintIndexBuildException("encoder batch-size mismatch");
                    }
                    long writeStarted = System.nanoTime();
                    for (int i = 0; i < embeddings.length; i++) {
                        writer.write(new ProductFingerprintRecord(accepted.get(i).sample.tuple(), embeddings[i]));
                        state.metrics.accepted(accepted.get(i).sample.tuple().reactionId());
                        state.accepted++;
                        recordsThisShard++;
                    }
                    state.metrics.writingNanos += System.nanoTime() - writeStarted;
                    reportProgress(state, started);
                }
                writer.flush();
            }
            commitShard(outputDirectory, partial, recordsThisShard, sampler, state);
        }
    }

    private List<Prepared> prepareCoverageBatch(ProductTupleSampler sampler,
            BuildState state, int capacity) {
        List<Prepared> result = new ArrayList<>(capacity);
        List<String> ids = sampler.reactionIds();
        while (result.size() < capacity && state.coverageIndex < ids.size()) {
            String reaction = ids.get(state.coverageIndex);
            long accepted = state.metrics.acceptedByReaction.getOrDefault(reaction, 0L);
            if (accepted >= config.sampling.minimumAcceptedPerReaction
                    || state.coverageAttempts >= config.sampling.coverageAttemptsPerReaction
                    || sampler.exhausted(reaction)) {
                if (accepted < config.sampling.minimumAcceptedPerReaction) {
                    state.metrics.rejected("REACTION_COVERAGE_NOT_REACHED");
                }
                state.coverageIndex++;
                state.coverageAttempts = 0;
                continue;
            }
            ProductTupleSampler.Sample sample = sampler.next(reaction);
            if (sample == null) continue;
            state.coverageAttempts++;
            state.metrics.proposals++;
            result.add(prepare(sample));
            // Coverage decisions depend on acceptance, so preserve strict order.
            if (result.get(result.size() - 1).accepted()) break;
        }
        return result;
    }

    private List<Prepared> prepareWeightedBatch(ProductTupleSampler sampler,
            BuildState state, ExecutorService workers, int capacity) {
        int window = capacity;
        List<Callable<Prepared>> tasks = new ArrayList<>(window);
        for (int i = 0; i < window; i++) {
            ProductTupleSampler.Sample sample = sampler.nextWeighted();
            if (sample == null) break;
            state.metrics.proposals++;
            tasks.add(() -> prepare(sample));
        }
        if (tasks.isEmpty()) return List.of();
        try {
            List<Future<Prepared>> futures = workers.invokeAll(tasks);
            List<Prepared> result = new ArrayList<>(futures.size());
            for (Future<Prepared> future : futures) result.add(future.get());
            return result;
        } catch (InterruptedException error) {
            Thread.currentThread().interrupt();
            throw new ProductFingerprintIndexBuildException("CPU preparation interrupted", error);
        } catch (ExecutionException error) {
            throw new ProductFingerprintIndexBuildException("CPU preparation failed", error.getCause());
        }
    }

    private Prepared prepare(ProductTupleSampler.Sample sample) {
        long assemblyStarted = System.nanoTime();
        StereoMolecule assembled;
        try {
            List<StereoMolecule> parts = new ArrayList<>(sample.synthons().size());
            for (RawSynthon synthon : sample.synthons()) {
                StereoMolecule molecule = new StereoMolecule();
                PARSER.get().parse(molecule, synthon.getIdcode());
                molecule.ensureHelperArrays(Molecule.cHelperCIP);
                parts.add(molecule);
            }
            assembled = SynthonAssembler.assembleSynthons_faster(parts);
        } catch (RuntimeException error) {
            return Prepared.rejected(sample, "ASSEMBLY_FAILED:" + error.getClass().getSimpleName(),
                    System.nanoTime() - assemblyStarted, 0);
        }
        long assemblyNanos = System.nanoTime() - assemblyStarted;
        int atoms = assembled.getAtoms();
        if (atoms < config.filters.minHeavyAtoms) {
            return Prepared.rejected(sample, "FILTER_TOO_FEW_HEAVY_ATOMS", assemblyNanos, 0);
        }
        if (atoms > config.filters.maxHeavyAtoms) {
            return Prepared.rejected(sample, "FILTER_TOO_MANY_HEAVY_ATOMS", assemblyNanos, 0);
        }
        if (config.filters.maxRotatableBonds > 0
                && countRotatableBonds(assembled) > config.filters.maxRotatableBonds) {
            return Prepared.rejected(sample, "FILTER_TOO_MANY_ROTATABLE_BONDS", assemblyNanos, 0);
        }
        long featureStarted = System.nanoTime();
        FeaturizationResult features = featurizer.featurize(assembled);
        long featureNanos = System.nanoTime() - featureStarted;
        return features.accepted()
                ? new Prepared(sample, features, null, assemblyNanos, featureNanos)
                : Prepared.rejected(sample, features.rejectionReason(), assemblyNanos, featureNanos);
    }

    private List<Prepared> account(List<Prepared> prepared, BuildState state) {
        List<Prepared> accepted = new ArrayList<>();
        for (Prepared item : prepared) {
            state.metrics.assemblyNanos += item.assemblyNanos;
            state.metrics.featurizationNanos += item.featurizationNanos;
            if (item.accepted()) accepted.add(item);
            else state.metrics.rejected(item.rejection);
        }
        return accepted;
    }

    private void commitShard(Path outputDirectory, Path partial, long records,
            ProductTupleSampler sampler, BuildState state) throws IOException {
        long hashStarted = System.nanoTime();
        String hash = sha256(partial);
        state.metrics.hashingNanos += System.nanoTime() - hashStarted;
        String name = shardName(state.shards.size());
        Path destination = outputDirectory.resolve(name);
        state.shards.add(new ProductFingerprintIndexShard(name, records, hash));
        state.sampler = sampler.snapshot();
        Path pending = outputDirectory.resolve(PENDING_STATE_FILE);
        writeJsonReplacing(pending, state);
        moveAtomic(partial, destination, false);
        moveAtomic(pending, outputDirectory.resolve(STATE_FILE), true);
    }

    private BuildState loadOrCreateState(Path outputDirectory, ProductTupleSampler sampler)
            throws IOException {
        recoverPendingCommit(outputDirectory);
        Path statePath = outputDirectory.resolve(STATE_FILE);
        if (Files.isRegularFile(statePath)) {
            if (!config.output.resume) throw new ProductFingerprintIndexBuildException(
                    "output contains an incomplete build and resume is disabled");
            BuildState state = MAPPER.readValue(statePath.toFile(), BuildState.class);
            validateState(state);
            return state;
        }
        try (var files = Files.list(outputDirectory)) {
            List<Path> unexpected = files.filter(path -> !path.getFileName().toString().endsWith(".partial"))
                    .toList();
            if (!unexpected.isEmpty()) throw new ProductFingerprintIndexBuildException(
                    "output directory is not empty or resumable: " + unexpected.get(0));
        }
        try (var files = Files.list(outputDirectory)) {
            for (Path partial : files.filter(path -> path.getFileName().toString().endsWith(".partial")).toList()) {
                Files.deleteIfExists(partial);
            }
        }
        BuildState state = new BuildState();
        state.configHash = provenance.configurationHash;
        state.rawspaceSha256 = provenance.rawspaceSha256;
        state.downsampledRawspaceSha256 = provenance.downsampledRawspaceSha256;
        state.modelBundleHash = provenance.modelBundleHash;
        state.sampler = sampler.snapshot();
        writeStateAtomic(outputDirectory, state);
        return state;
    }

    private void recoverPendingCommit(Path outputDirectory) throws IOException {
        Path pending = outputDirectory.resolve(PENDING_STATE_FILE);
        if (!Files.isRegularFile(pending)) return;
        BuildState state = MAPPER.readValue(pending.toFile(), BuildState.class);
        validateState(state);
        if (state.shards.isEmpty()) throw new ProductFingerprintIndexBuildException(
                "pending checkpoint contains no shard");
        ProductFingerprintIndexShard shard = state.shards.get(state.shards.size() - 1);
        Path destination = outputDirectory.resolve(shard.path());
        Path partial = outputDirectory.resolve(shard.path() + ".partial");
        if (!Files.isRegularFile(destination)) {
            if (!Files.isRegularFile(partial) || !sha256(partial).equals(shard.sha256())) {
                throw new ProductFingerprintIndexBuildException("pending shard is missing or corrupt");
            }
            moveAtomic(partial, destination, false);
        }
        if (!sha256(destination).equals(shard.sha256())) {
            throw new ProductFingerprintIndexBuildException("pending finalized shard checksum mismatch");
        }
        moveAtomic(pending, outputDirectory.resolve(STATE_FILE), true);
    }

    private void validateState(BuildState state) {
        if (state.formatVersion != 1
                || !provenance.configurationHash.equals(state.configHash)
                || !provenance.rawspaceSha256.equals(state.rawspaceSha256)
                || !provenance.downsampledRawspaceSha256.equals(state.downsampledRawspaceSha256)
                || !provenance.modelBundleHash.equals(state.modelBundleHash)) {
            throw new ProductFingerprintIndexBuildException(
                    "checkpoint configuration, inputs, or model do not match this run");
        }
    }

    private void verifyCommittedShards(Path outputDirectory,
            List<ProductFingerprintIndexShard> shards) throws IOException {
        long count = 0;
        for (ProductFingerprintIndexShard shard : shards) {
            Path path = outputDirectory.resolve(shard.path());
            if (!Files.isRegularFile(path) || !sha256(path).equals(shard.sha256())) {
                throw new ProductFingerprintIndexBuildException(
                        "committed shard is missing or corrupt: " + shard.path());
            }
            count += shard.recordCount();
        }
        if (count == 0 && !shards.isEmpty()) throw new ProductFingerprintIndexBuildException(
                "invalid zero-sized committed shards");
    }

    private ProductFingerprintIndexManifest createManifest(BuildState state, long started) {
        ProductFingerprintIndexManifest manifest = new ProductFingerprintIndexManifest();
        manifest.rawspaceIdentity = provenance.rawspaceIdentity;
        manifest.rawspaceSha256 = provenance.rawspaceSha256;
        manifest.downsampledRawspaceIdentity = provenance.downsampledRawspaceIdentity;
        manifest.downsampledRawspaceSha256 = provenance.downsampledRawspaceSha256;
        manifest.modelBundleHash = provenance.modelBundleHash;
        manifest.samplingAlgorithm = ProductTupleSampler.ALGORITHM;
        manifest.samplingSeed = config.sampling.seed;
        manifest.reactionWeightConfiguration = MAPPER.convertValue(
                config.sampling.reactionWeighting, new TypeReference<>() {});
        manifest.molecularFilters = MAPPER.convertValue(config.filters, new TypeReference<>() {});
        manifest.recordCount = state.accepted;
        manifest.shards = List.copyOf(state.shards);
        manifest.buildConfigurationHash = provenance.configurationHash;
        manifest.buildStatistics = state.metrics.manifestView(System.nanoTime() - started);
        manifest.validate();
        return manifest;
    }

    private ProductFingerprintIndexBuildResult completedResult(Path manifestPath) throws IOException {
        ProductFingerprintIndexManifest manifest = MAPPER.readValue(manifestPath.toFile(),
                ProductFingerprintIndexManifest.class);
        manifest.validate();
        verifyCommittedShards(manifestPath.getParent(), manifest.shards);
        if (!provenance.configurationHash.equals(manifest.buildConfigurationHash)
                || !provenance.rawspaceSha256.equals(manifest.rawspaceSha256)
                || !provenance.downsampledRawspaceSha256.equals(manifest.downsampledRawspaceSha256)
                || !provenance.modelBundleHash.equals(manifest.modelBundleHash)) {
            throw new ProductFingerprintIndexBuildException(
                    "completed output does not match this configuration, input, and model");
        }
        ProductFingerprintIndexBuildMetrics metrics = new ProductFingerprintIndexBuildMetrics();
        metrics.accepted = manifest.recordCount;
        return new ProductFingerprintIndexBuildResult(manifestPath, manifest.recordCount,
                manifest.shards.size(), true, metrics);
    }

    private void writeStateAtomic(Path outputDirectory, BuildState state) throws IOException {
        Path temporary = outputDirectory.resolve(STATE_FILE + ".partial");
        writeJsonAtomic(temporary, outputDirectory.resolve(STATE_FILE), state);
    }

    private static void writeJsonAtomic(Path temporary, Path destination, Object value)
            throws IOException {
        MAPPER.writeValue(temporary.toFile(), value);
        moveAtomic(temporary, destination, true);
    }

    private static void writeJsonReplacing(Path destination, Object value) throws IOException {
        Path temporary = destination.resolveSibling(destination.getFileName() + ".partial");
        MAPPER.writeValue(temporary.toFile(), value);
        moveAtomic(temporary, destination, true);
    }

    private static void moveAtomic(Path source, Path destination, boolean replace) throws IOException {
        try {
            if (replace) Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE,
                    StandardCopyOption.REPLACE_EXISTING);
            else Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE);
        } catch (AtomicMoveNotSupportedException error) {
            throw new IOException("atomic moves are required in the output directory", error);
        }
    }

    private static int countRotatableBonds(StereoMolecule molecule) {
        molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
        boolean[] rotatable = new boolean[molecule.getBonds()];
        TorsionDB.findRotatableBonds(molecule, true, rotatable);
        int count = 0;
        for (boolean value : rotatable) if (value) count++;
        return count;
    }

    private void reportProgress(BuildState state, long started) {
        if (config.runtime.progressIntervalSeconds == 0) return;
        long now = System.nanoTime();
        if (now - state.lastProgressNanos < Duration.ofSeconds(
                config.runtime.progressIntervalSeconds).toNanos()) return;
        state.lastProgressNanos = now;
        double seconds = (now - started) / 1_000_000_000d;
        System.out.printf("Hyperspace3D index: %,d/%,d accepted, %,d proposals, %.1f accepted/s%n",
                state.accepted, config.output.targetRecordCount, state.metrics.proposals,
                state.accepted / Math.max(0.001, seconds));
    }

    private static String shardName(int index) {
        return String.format("shard-%05d.h3di", index);
    }

    public static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int read; (read = input.read(buffer)) >= 0;) digest.update(buffer, 0, read);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }

    public static String configurationHash(ProductFingerprintIndexBuildConfig config) {
        try {
            Map<String, Object> semantic = new LinkedHashMap<>();
            semantic.put("formatVersion", config.formatVersion);
            semantic.put("targetRecordCount", config.output.targetRecordCount);
            semantic.put("recordsPerShard", config.output.recordsPerShard);
            semantic.put("sampling", config.sampling);
            semantic.put("filters", config.filters);
            semantic.put("device", config.runtime.device);
            semantic.put("cudaDeviceId", config.runtime.cudaDeviceId);
            semantic.put("encoderBatchSize", config.runtime.encoderBatchSize);
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            return HexFormat.of().formatHex(digest.digest(MAPPER.writeValueAsBytes(semantic)));
        } catch (IOException | NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }

    public static String identity(RawSynthonSpace space) {
        return space.getName() + ":" + space.getVersion();
    }

    public record Provenance(String rawspaceIdentity, String rawspaceSha256,
            String downsampledRawspaceIdentity, String downsampledRawspaceSha256,
            String modelBundleHash, String configurationHash) {}

    public static final class BuildState {
        public int formatVersion = 1;
        public String status = "IN_PROGRESS";
        public String configHash;
        public String rawspaceSha256;
        public String downsampledRawspaceSha256;
        public String modelBundleHash;
        public long accepted;
        public int coverageIndex;
        public int coverageAttempts;
        public ProductTupleSampler.Snapshot sampler;
        public List<ProductFingerprintIndexShard> shards = new ArrayList<>();
        public ProductFingerprintIndexBuildMetrics metrics = new ProductFingerprintIndexBuildMetrics();
        public transient long lastProgressNanos;
    }

    private record Prepared(ProductTupleSampler.Sample sample, FeaturizationResult features,
            String rejection, long assemblyNanos, long featurizationNanos) {
        private boolean accepted() { return rejection == null; }
        private static Prepared rejected(ProductTupleSampler.Sample sample, String reason,
                long assemblyNanos, long featurizationNanos) {
            return new Prepared(sample, null, reason, assemblyNanos, featurizationNanos);
        }
    }
}
