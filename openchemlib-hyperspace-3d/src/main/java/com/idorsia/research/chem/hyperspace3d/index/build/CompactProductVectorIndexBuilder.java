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
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexShard;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexWriter;
import com.idorsia.research.chem.hyperspace3d.model.ProductEmbeddingBatchEncoder;
import java.io.IOException;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/** Direct graph -> 128D -> 16D builder. It never materializes an intermediate 128D index. */
public final class CompactProductVectorIndexBuilder {
    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT).enable(SerializationFeature.ORDER_MAP_ENTRIES_BY_KEYS);
    private static final ThreadLocal<IDCodeParser> PARSER = ThreadLocal.withInitial(IDCodeParser::new);
    private final RawSynthonSpace full;
    private final RawSynthonSpace downsampled;
    private final ProductEmbeddingBatchEncoder encoder;
    private final CompactProductVectorIndexBuildConfig config;
    private final Provenance provenance;
    private final DeepSpaceTensorBatchBuilder tensors = new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer());

    public CompactProductVectorIndexBuilder(RawSynthonSpace full, RawSynthonSpace downsampled,
            ProductEmbeddingBatchEncoder encoder, CompactProductVectorIndexBuildConfig config,
            Provenance provenance) {
        this.full = full; this.downsampled = downsampled; this.encoder = encoder;
        this.config = config; this.provenance = provenance; config.validate();
    }

    public ProductFingerprintIndexBuildResult build(Path directory) {
        long started = System.nanoTime();
        try {
            Files.createDirectories(directory);
            try (var contents = Files.list(directory)) {
                if (contents.findAny().isPresent()) throw new ProductFingerprintIndexBuildException(
                        "compact index output directory must be empty");
            }
            ProductTupleSampler sampler = new ProductTupleSampler(full, downsampled,
                    config.sampling.seed, config.sampling.reactionWeighting.toWeighting());
            List<String> reactionDictionary = sampler.reactionIds();
            long minimum = (long) reactionDictionary.size() * config.sampling.minimumAcceptedPerReaction;
            if (config.output.targetRecordCount < minimum) throw new ProductFingerprintIndexBuildException(
                    "targetRecordCount is smaller than the reaction coverage floor");
            var metrics = new ProductFingerprintIndexBuildMetrics();
            List<ProductVectorIndexShard> shards = new ArrayList<>();
            long accepted = 0;
            try (var workers = Executors.newFixedThreadPool(config.runtime.cpuWorkers);
                 var sink = new ShardSink(directory, reactionDictionary, shards, metrics)) {
                for (int reactionIndex = 0; reactionIndex < reactionDictionary.size()
                        && accepted < config.output.targetRecordCount; reactionIndex++) {
                    String reaction = reactionDictionary.get(reactionIndex); int reactionAccepted = 0; int attempts = 0;
                    while (reactionAccepted < config.sampling.minimumAcceptedPerReaction
                            && attempts++ < config.sampling.coverageAttemptsPerReaction) {
                        var sample = sampler.next(reaction); if (sample == null) break;
                        metrics.proposals++; Prepared prepared = prepare(sample);
                        if (prepared.features == null) { accountRejected(metrics, prepared); continue; }
                        accepted += writeBatch(sink, List.of(prepared), 1, metrics);
                        reactionAccepted++;
                    }
                    if (reactionAccepted < config.sampling.minimumAcceptedPerReaction)
                        throw new ProductFingerprintIndexBuildException("reaction coverage not reached: " + reaction);
                }
                while (accepted < config.output.targetRecordCount) {
                    int batchLimit = (int) Math.min(config.runtime.encoderBatchSize,
                            config.output.targetRecordCount - accepted);
                    List<ProductTupleSampler.Sample> samples = new ArrayList<>(batchLimit);
                    for (int i = 0; i < batchLimit; i++) {
                        var sample = sampler.nextWeighted(); if (sample == null) break;
                        samples.add(sample); metrics.proposals++;
                    }
                    if (samples.isEmpty()) throw new ProductFingerprintIndexBuildException("sampled tuple space exhausted");
                    List<Future<Prepared>> futures = new ArrayList<>();
                    for (var sample : samples) futures.add(workers.submit(() -> prepare(sample)));
                    List<Prepared> prepared = new ArrayList<>();
                    for (Future<Prepared> future : futures) {
                        Prepared item = future.get();
                        if (item.features == null) accountRejected(metrics, item); else prepared.add(item);
                    }
                    if (prepared.isEmpty()) continue;
                    accepted += writeBatch(sink, prepared,
                            config.output.targetRecordCount - accepted, metrics);
                }
            }
            ProductVectorIndexManifest manifest = manifest(reactionDictionary, shards, accepted, metrics,
                    System.nanoTime() - started);
            Path partial = directory.resolve("manifest.json.partial");
            MAPPER.writeValue(partial.toFile(), manifest);
            move(partial, directory.resolve("manifest.json"));
            return new ProductFingerprintIndexBuildResult(directory.resolve("manifest.json"), accepted,
                    shards.size(), false, metrics);
        } catch (ProductFingerprintIndexBuildException error) { throw error; }
        catch (Exception error) { throw new ProductFingerprintIndexBuildException("compact index build failed", error); }
    }

    private long writeBatch(ShardSink sink, List<Prepared> prepared, long maximum,
            ProductFingerprintIndexBuildMetrics metrics) throws IOException {
        int take = (int) Math.min(maximum, prepared.size());
        long pack = System.nanoTime();
        var batch = tensors.buildFromFeatures(prepared.subList(0, take).stream().map(Prepared::features).toList());
        metrics.tensorPackingNanos += System.nanoTime() - pack;
        long encode = System.nanoTime(); float[][] vectors = encoder.encode(batch);
        metrics.encodingNanos += System.nanoTime() - encode;
        if (vectors.length != take) throw new ProductFingerprintIndexBuildException("encoder batch-size mismatch");
        long write = System.nanoTime();
        for (int i = 0; i < take; i++) {
            Prepared item = prepared.get(i);
            metrics.assemblyNanos += item.assemblyNanos;
            metrics.featurizationNanos += item.featurizationNanos;
            sink.write(item.sample.tuple(), vectors[i]);
            metrics.accepted(item.sample.tuple().reactionId());
        }
        metrics.writingNanos += System.nanoTime() - write;
        return take;
    }

    private Prepared prepare(ProductTupleSampler.Sample sample) {
        long assemblyStarted = System.nanoTime();
        StereoMolecule assembled;
        try {
            List<StereoMolecule> parts = new ArrayList<>();
            for (RawSynthon synthon : sample.synthons()) {
                StereoMolecule part = new StereoMolecule(); PARSER.get().parse(part, synthon.getIdcode());
                part.ensureHelperArrays(Molecule.cHelperCIP); parts.add(part);
            }
            assembled = SynthonAssembler.assembleSynthons_faster(parts);
        } catch (RuntimeException error) {
            return new Prepared(sample, null, "ASSEMBLY_FAILED", System.nanoTime() - assemblyStarted, 0);
        }
        long assemblyNanos = System.nanoTime() - assemblyStarted;
        if (assembled.getAtoms() < config.filters.minHeavyAtoms)
            return new Prepared(sample, null, "FILTER_TOO_FEW_HEAVY_ATOMS", assemblyNanos, 0);
        if (assembled.getAtoms() > config.filters.maxHeavyAtoms)
            return new Prepared(sample, null, "FILTER_TOO_MANY_HEAVY_ATOMS", assemblyNanos, 0);
        if (config.filters.maxRotatableBonds > 0 && rotatable(assembled) > config.filters.maxRotatableBonds)
            return new Prepared(sample, null, "FILTER_TOO_MANY_ROTATABLE_BONDS", assemblyNanos, 0);
        long featureStarted = System.nanoTime(); FeaturizationResult features = new OCLDeepSpaceFeaturizer().featurize(assembled);
        long featureNanos = System.nanoTime() - featureStarted;
        return features.accepted() ? new Prepared(sample, features, null, assemblyNanos, featureNanos)
                : new Prepared(sample, null, features.rejectionReason(), assemblyNanos, featureNanos);
    }

    private static void accountRejected(ProductFingerprintIndexBuildMetrics metrics, Prepared item) {
        metrics.assemblyNanos += item.assemblyNanos; metrics.featurizationNanos += item.featurizationNanos;
        metrics.rejected(item.rejection);
    }
    private static int rotatable(StereoMolecule molecule) {
        molecule.ensureHelperArrays(Molecule.cHelperNeighbours); boolean[] flags = new boolean[molecule.getBonds()];
        TorsionDB.findRotatableBonds(molecule, true, flags); int count = 0; for (boolean flag : flags) if (flag) count++;
        return count;
    }
    private ProductVectorIndexManifest manifest(List<String> dictionary,
            List<ProductVectorIndexShard> shards, long count, ProductFingerprintIndexBuildMetrics metrics,
            long elapsed) {
        ProductVectorIndexManifest result = new ProductVectorIndexManifest();
        result.rawspaceIdentity = provenance.rawspaceIdentity; result.rawspaceSha256 = provenance.rawspaceSha256;
        result.downsampledRawspaceIdentity = provenance.downsampledRawspaceIdentity;
        result.downsampledRawspaceSha256 = provenance.downsampledRawspaceSha256;
        result.sourceModelBundleHash = provenance.sourceModelBundleHash;
        result.projectionModelBundleHash = provenance.projectionModelBundleHash;
        result.modelVersion = provenance.modelVersion; result.samplingAlgorithm = ProductTupleSampler.ALGORITHM;
        result.samplingSeed = config.sampling.seed;
        result.reactionWeightConfiguration = MAPPER.convertValue(config.sampling.reactionWeighting, new TypeReference<>() {});
        result.molecularFilters = MAPPER.convertValue(config.filters, new TypeReference<>() {});
        result.recordCount = count; result.reactionDictionary = List.copyOf(dictionary);
        result.shards = List.copyOf(shards); result.buildConfigurationHash = provenance.configurationHash;
        result.buildStatistics = metrics.manifestView(elapsed); result.validate(); return result;
    }
    private static void move(Path source, Path destination) throws IOException {
        try { Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE); }
        catch (AtomicMoveNotSupportedException error) { throw new IOException("atomic moves are required", error); }
    }
    private final class ShardSink implements AutoCloseable {
        private final Path directory;
        private final List<String> dictionary;
        private final List<ProductVectorIndexShard> shards;
        private final ProductFingerprintIndexBuildMetrics metrics;
        private ProductVectorIndexWriter writer;
        private Path vectors;
        private Path rows;
        private Path tuples;
        private long records;

        private ShardSink(Path directory, List<String> dictionary,
                List<ProductVectorIndexShard> shards, ProductFingerprintIndexBuildMetrics metrics) {
            this.directory = directory; this.dictionary = dictionary;
            this.shards = shards; this.metrics = metrics;
        }
        private void write(com.idorsia.research.chem.hyperspace3d.index.ProductTuple tuple,
                float[] vector) throws IOException {
            if (writer == null) open();
            writer.write(tuple, vector); records++;
            if (records == config.output.recordsPerShard) commit();
        }
        private void open() throws IOException {
            String stem = String.format("shard-%05d", shards.size());
            vectors = directory.resolve(stem + ".vec.partial");
            rows = directory.resolve(stem + ".rows.partial");
            tuples = directory.resolve(stem + ".tuples.partial");
            writer = new ProductVectorIndexWriter(vectors, rows, tuples, dictionary);
        }
        private void commit() throws IOException {
            if (writer == null) return;
            writer.close();
            long hashStarted = System.nanoTime();
            String vh = ProductFingerprintIndexBuilder.sha256(vectors);
            String rh = ProductFingerprintIndexBuilder.sha256(rows);
            String th = ProductFingerprintIndexBuilder.sha256(tuples);
            metrics.hashingNanos += System.nanoTime() - hashStarted;
            String stem = String.format("shard-%05d", shards.size());
            String vn = stem + ".vec", rn = stem + ".rows", tn = stem + ".tuples";
            move(vectors, directory.resolve(vn)); move(rows, directory.resolve(rn));
            move(tuples, directory.resolve(tn));
            shards.add(new ProductVectorIndexShard(vn, rn, tn, records, vh, rh, th));
            writer = null; records = 0;
        }
        @Override public void close() throws IOException { commit(); }
    }

    public static String configurationHash(CompactProductVectorIndexBuildConfig config) {
        try {
            Map<String,Object> value = new LinkedHashMap<>();
            value.put("formatVersion", config.formatVersion); value.put("output", config.output);
            value.put("sampling", config.sampling); value.put("filters", config.filters);
            value.put("runtimeDevice", config.runtime.device); value.put("batch", config.runtime.encoderBatchSize);
            Path temp = Files.createTempFile("h3d-compact-config-", ".json");
            try { MAPPER.writeValue(temp.toFile(), value); return ProductFingerprintIndexBuilder.sha256(temp); }
            finally { Files.deleteIfExists(temp); }
        } catch (IOException error) { throw new IllegalStateException(error); }
    }
    public record Provenance(String rawspaceIdentity, String rawspaceSha256,
            String downsampledRawspaceIdentity, String downsampledRawspaceSha256,
            String sourceModelBundleHash, String projectionModelBundleHash,
            String modelVersion, String configurationHash) {}
    private record Prepared(ProductTupleSampler.Sample sample, FeaturizationResult features,
            String rejection, long assemblyNanos, long featurizationNanos) {}
}
