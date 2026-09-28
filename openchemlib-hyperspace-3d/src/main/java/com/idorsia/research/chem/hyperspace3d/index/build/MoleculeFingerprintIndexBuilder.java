package com.idorsia.research.chem.hyperspace3d.index.build;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.FeaturizationResult;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexShard;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexWriter;
import com.idorsia.research.chem.hyperspace3d.model.MoleculeFingerprintBatch;
import com.idorsia.research.chem.hyperspace3d.model.MoleculeFingerprintBatchEncoder;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.Duration;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

/**
 * Streams a flat molecule table through parallel OCL preparation and a dual-output
 * ONNX encoder. Source-order shards are independently committed and resumable.
 */
public final class MoleculeFingerprintIndexBuilder {
    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT)
            .enable(SerializationFeature.ORDER_MAP_ENTRIES_BY_KEYS);
    private static final String STATE_FILE = "build-state.json";
    private static final String MANIFEST_FILE = "manifest.json";
    private static final ThreadLocal<SmilesParser> PARSER =
            ThreadLocal.withInitial(SmilesParser::new);

    private final MoleculeFingerprintBatchEncoder encoder;
    private final MoleculeFingerprintIndexBuildConfig config;
    private final MoleculeFingerprintIndexBuildConfig.ResolvedPaths paths;
    private final OCLDeepSpaceFeaturizer featurizer = new OCLDeepSpaceFeaturizer();
    private final DeepSpaceTensorBatchBuilder tensorBuilder =
            new DeepSpaceTensorBatchBuilder(featurizer);

    public MoleculeFingerprintIndexBuilder(MoleculeFingerprintBatchEncoder encoder,
            MoleculeFingerprintIndexBuildConfig config,
            MoleculeFingerprintIndexBuildConfig.ResolvedPaths paths) {
        this.encoder = Objects.requireNonNull(encoder);
        this.config = Objects.requireNonNull(config);
        this.paths = Objects.requireNonNull(paths);
        config.validate();
    }

    public MoleculeFingerprintIndexBuildResult build() {
        long started = System.nanoTime();
        Path output = paths.outputDirectory();
        try {
            requireInputFiles();
            Files.createDirectories(output);
            Path manifestPath = output.resolve(MANIFEST_FILE);
            if (Files.isRegularFile(manifestPath)) {
                return completedResult(manifestPath);
            }
            cleanupPartialDirectories(output);
            BuildState state = loadOrCreateState(output);
            recoverCommittedShards(output, state);
            boolean resumed = state.nextSourceRow > 0;
            int queueSize = Math.max(config.runtime.cpuWorkers,
                    config.runtime.encoderBatchSize * config.runtime.queueCapacity);
            ExecutorService workers = new ThreadPoolExecutor(
                    config.runtime.cpuWorkers, config.runtime.cpuWorkers,
                    0L, TimeUnit.MILLISECONDS, new ArrayBlockingQueue<>(queueSize),
                    new ThreadPoolExecutor.CallerRunsPolicy());
            try (BufferedReader reader = openTable(paths.library())) {
                Header header = readHeader(reader);
                skipRows(reader, state.nextSourceRow);
                buildRemaining(reader, header, output, state, workers, started);
            } finally {
                workers.shutdownNow();
            }
            MoleculeFingerprintIndexManifest manifest = manifest(state, started);
            manifest.validate();
            writeJsonAtomic(output.resolve(MANIFEST_FILE + ".partial"),
                    output.resolve(MANIFEST_FILE), manifest);
            state.status = "COMPLETE";
            writeState(output, state);
            return new MoleculeFingerprintIndexBuildResult(
                    output.resolve(MANIFEST_FILE), state.nextSourceRow,
                    state.metrics.accepted, state.shards.size(), resumed, state.metrics);
        } catch (IOException | RuntimeException error) {
            throw new IllegalStateException("molecule fingerprint index build failed", error);
        }
    }

    private void buildRemaining(BufferedReader reader, Header header, Path output,
            BuildState state, ExecutorService workers, long started) throws IOException {
        boolean eof = false;
        while (!eof) {
            int shardIndex = state.shards.size();
            String name = shardName(shardIndex);
            Path partial = output.resolve(name + ".partial");
            Path destination = output.resolve(name);
            deleteTree(partial);
            Files.createDirectories(partial);
            ShardRun run;
            try (MoleculeFingerprintIndexWriter writer =
                         new MoleculeFingerprintIndexWriter(partial)) {
                run = processShard(reader, header, state.nextSourceRow, writer, workers,
                        state.metrics, started);
            }
            if (run.sourceRows == 0) {
                deleteTree(partial);
                break;
            }
            MoleculeFingerprintIndexShard shard = new MoleculeFingerprintIndexShard();
            shard.shardIndex = shardIndex;
            shard.directory = name;
            shard.firstSourceRow = state.nextSourceRow;
            shard.sourceRowCount = run.sourceRows;
            shard.recordCount = run.accepted;
            shard.rejectedCount = run.sourceRows - run.accepted;
            shard.rejections = run.rejections;
            MAPPER.writeValue(partial.resolve("shard.json").toFile(), shard);
            Files.createFile(partial.resolve(".complete"));
            moveAtomic(partial, destination, false);
            state.shards.add(shard);
            state.nextSourceRow += shard.sourceRowCount;
            writeState(output, state);
            eof = run.eof;
            reportProgress(state, started);
        }
    }

    private ShardRun processShard(BufferedReader reader, Header header,
            long firstSourceRow, MoleculeFingerprintIndexWriter writer,
            ExecutorService workers, MoleculeFingerprintIndexBuildMetrics metrics,
            long started) throws IOException {
        Deque<PendingBatch> pending = new ArrayDeque<>();
        long queuedRows = 0;
        boolean eof = false;
        for (int i = 0; i < config.runtime.queueCapacity && !eof
                && queuedRows < config.output.sourceRowsPerShard; i++) {
            PendingBatch batch = readAndSubmit(reader, header,
                    firstSourceRow + queuedRows,
                    config.output.sourceRowsPerShard - queuedRows, workers);
            if (batch == null) eof = true;
            else {
                pending.addLast(batch);
                queuedRows += batch.sourceRows;
            }
        }

        long processed = 0;
        long accepted = 0;
        Map<String, Long> shardRejections = new LinkedHashMap<>();
        while (!pending.isEmpty()) {
            PendingBatch current = pending.removeFirst();
            long waitStarted = System.nanoTime();
            List<Prepared> prepared = await(current.futures);
            metrics.encoderWaitNanos += System.nanoTime() - waitStarted;

            if (!eof && queuedRows < config.output.sourceRowsPerShard) {
                PendingBatch next = readAndSubmit(reader, header,
                        firstSourceRow + queuedRows,
                        config.output.sourceRowsPerShard - queuedRows, workers);
                if (next == null) eof = true;
                else {
                    pending.addLast(next);
                    queuedRows += next.sourceRows;
                }
            }

            List<Prepared> acceptedRows = new ArrayList<>(prepared.size());
            List<FeaturizationResult> features = new ArrayList<>(prepared.size());
            for (Prepared row : prepared) {
                processed++;
                metrics.sourceRows++;
                metrics.parsingAndFeaturizationNanos += row.preparationNanos;
                if (row.accepted()) {
                    acceptedRows.add(row);
                    features.add(row.features);
                    metrics.accepted++;
                    accepted++;
                } else {
                    metrics.rejected(row.rejection);
                    shardRejections.merge(row.rejection, 1L, Long::sum);
                }
            }
            if (!features.isEmpty()) {
                long packingStarted = System.nanoTime();
                DeepSpaceTensorBatch tensors = tensorBuilder.buildFromFeatures(features);
                metrics.tensorPackingNanos += System.nanoTime() - packingStarted;
                long encodingStarted = System.nanoTime();
                MoleculeFingerprintBatch embeddings = encoder.encode(tensors);
                metrics.encodingNanos += System.nanoTime() - encodingStarted;
                if (embeddings.size() != acceptedRows.size()) {
                    throw new IllegalStateException("encoder changed the batch row count");
                }
                long writingStarted = System.nanoTime();
                for (int row = 0; row < acceptedRows.size(); row++) {
                    Prepared molecule = acceptedRows.get(row);
                    writer.write(molecule.sourceRow, molecule.heavyAtoms,
                            molecule.id, molecule.smiles,
                            embeddings.base128()[row], embeddings.compact16()[row]);
                }
                metrics.writingNanos += System.nanoTime() - writingStarted;
            }
            reportProgress(metrics, started);
        }
        return new ShardRun(processed, accepted, eof, shardRejections);
    }

    private PendingBatch readAndSubmit(BufferedReader reader, Header header,
            long firstSourceRow, long remaining, ExecutorService workers)
            throws IOException {
        int count = (int) Math.min(config.runtime.encoderBatchSize, remaining);
        List<Future<Prepared>> futures = new ArrayList<>(count);
        for (int index = 0; index < count; index++) {
            String line = reader.readLine();
            if (line == null) break;
            long sourceRow = firstSourceRow + index;
            RawRecord raw = rawRecord(line, sourceRow, header);
            futures.add(workers.submit(() -> prepare(raw)));
        }
        return futures.isEmpty() ? null : new PendingBatch(futures, futures.size());
    }

    private Prepared prepare(RawRecord raw) {
        long started = System.nanoTime();
        if (raw.rejection != null) {
            return Prepared.rejected(raw, raw.rejection, System.nanoTime() - started);
        }
        try {
            StereoMolecule molecule = new StereoMolecule();
            PARSER.get().parse(molecule, raw.smiles);
            FeaturizationResult features = featurizer.featurize(molecule);
            if (!features.accepted()) {
                return Prepared.rejected(raw, features.rejectionReason(),
                        System.nanoTime() - started);
            }
            int atoms = 0;
            for (boolean present : features.atomMask()) if (present) atoms++;
            return new Prepared(raw.sourceRow, raw.id, raw.smiles, features,
                    atoms, null, System.nanoTime() - started);
        } catch (Exception error) {
            return Prepared.rejected(raw, "SMILES_PARSE_ERROR",
                    System.nanoTime() - started);
        }
    }

    private static RawRecord rawRecord(String line, long sourceRow, Header header) {
        if (line.isBlank()) return new RawRecord(sourceRow, "", "", "BLANK_ROW");
        String[] fields = line.split("\\t", -1);
        int maximum = Math.max(header.smilesIndex, header.idIndex);
        if (fields.length <= maximum) {
            return new RawRecord(sourceRow, "", "", "MISSING_COLUMN");
        }
        String smiles = fields[header.smilesIndex];
        String id = fields[header.idIndex];
        if (smiles.isBlank()) return new RawRecord(sourceRow, id, smiles, "MISSING_SMILES");
        if (id.isBlank()) return new RawRecord(sourceRow, id, smiles, "MISSING_ID");
        return new RawRecord(sourceRow, id, smiles, null);
    }

    private static List<Prepared> await(List<Future<Prepared>> futures) {
        List<Prepared> result = new ArrayList<>(futures.size());
        for (Future<Prepared> future : futures) {
            try {
                result.add(future.get());
            } catch (InterruptedException error) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException("molecule preparation interrupted", error);
            } catch (ExecutionException error) {
                throw new IllegalStateException("molecule preparation failed", error.getCause());
            }
        }
        return result;
    }

    private BufferedReader openTable(Path input) throws IOException {
        InputStream stream = new BufferedInputStream(Files.newInputStream(input), 1 << 20);
        String name = input.getFileName().toString().toLowerCase();
        try {
            if (name.endsWith(".bz2")) {
                stream = new BZip2CompressorInputStream(stream, true);
            } else if (name.endsWith(".gz")) {
                stream = new GZIPInputStream(stream, 1 << 20);
            }
            return new BufferedReader(
                    new InputStreamReader(stream, StandardCharsets.UTF_8), 1 << 20);
        } catch (IOException error) {
            stream.close();
            throw error;
        }
    }

    private Header readHeader(BufferedReader reader) throws IOException {
        String line = reader.readLine();
        if (line == null) throw new IOException("input table is empty");
        String[] names = line.split("\\t", -1);
        int smiles = findColumn(names, config.inputs.smilesColumn);
        int id = findColumn(names, config.inputs.idColumn);
        if (smiles < 0 || id < 0) {
            throw new IOException("input header does not contain configured SMILES and ID columns");
        }
        return new Header(smiles, id);
    }

    private static int findColumn(String[] names, String wanted) {
        for (int index = 0; index < names.length; index++) {
            if (names[index].equals(wanted)) return index;
        }
        return -1;
    }

    private static void skipRows(BufferedReader reader, long count) throws IOException {
        for (long row = 0; row < count; row++) {
            if (reader.readLine() == null) {
                throw new IOException("input ended before the resume source row");
            }
        }
    }

    private BuildState loadOrCreateState(Path output) throws IOException {
        Path statePath = output.resolve(STATE_FILE);
        if (Files.isRegularFile(statePath)) {
            if (!config.output.resume) {
                throw new IOException("output contains an incomplete build and resume is disabled");
            }
            BuildState state = MAPPER.readValue(statePath.toFile(), BuildState.class);
            validateState(state);
            return state;
        }
        try (var entries = Files.list(output)) {
            List<Path> unexpected = entries.toList();
            if (!unexpected.isEmpty()) {
                throw new IOException("output directory is not empty or resumable: "
                        + unexpected.get(0));
            }
        }
        BuildState state = new BuildState();
        state.input = paths.library().toString();
        state.smilesColumn = config.inputs.smilesColumn;
        state.idColumn = config.inputs.idColumn;
        state.modelBundle = paths.modelBundle().toString();
        state.compactBundle = paths.compactBundle().toString();
        state.sourceRowsPerShard = config.output.sourceRowsPerShard;
        writeState(output, state);
        return state;
    }

    private void validateState(BuildState state) throws IOException {
        if (state.formatVersion != 1
                || !paths.library().toString().equals(state.input)
                || !config.inputs.smilesColumn.equals(state.smilesColumn)
                || !config.inputs.idColumn.equals(state.idColumn)
                || !paths.modelBundle().toString().equals(state.modelBundle)
                || !paths.compactBundle().toString().equals(state.compactBundle)
                || config.output.sourceRowsPerShard != state.sourceRowsPerShard) {
            throw new IOException("resume state does not match the structural build settings");
        }
    }

    private void recoverCommittedShards(Path output, BuildState state) throws IOException {
        boolean changed = false;
        while (true) {
            int index = state.shards.size();
            Path directory = output.resolve(shardName(index));
            if (!Files.isDirectory(directory)
                    || !Files.isRegularFile(directory.resolve(".complete"))) break;
            MoleculeFingerprintIndexShard shard = MAPPER.readValue(
                    directory.resolve("shard.json").toFile(),
                    MoleculeFingerprintIndexShard.class);
            if (shard.shardIndex != index || shard.firstSourceRow != state.nextSourceRow
                    || !shardName(index).equals(shard.directory)) {
                throw new IOException("committed shard metadata is not contiguous");
            }
            state.shards.add(shard);
            applyRecoveredShard(state, shard);
            changed = true;
        }
        if (changed) writeState(output, state);
    }

    private static void applyRecoveredShard(BuildState state,
            MoleculeFingerprintIndexShard shard) {
        state.nextSourceRow += shard.sourceRowCount;
        state.metrics.sourceRows += shard.sourceRowCount;
        state.metrics.accepted += shard.recordCount;
        state.metrics.rejected += shard.rejectedCount;
        shard.rejections.forEach((reason, count) ->
                state.metrics.rejections.merge(reason, count, Long::sum));
    }

    private MoleculeFingerprintIndexManifest manifest(BuildState state, long started) {
        MoleculeFingerprintIndexManifest result = new MoleculeFingerprintIndexManifest();
        result.input = paths.library().toString();
        result.smilesColumn = config.inputs.smilesColumn;
        result.idColumn = config.inputs.idColumn;
        result.modelBundle = paths.modelBundle().toString();
        result.compactBundle = paths.compactBundle().toString();
        result.sourceRowsPerShard = config.output.sourceRowsPerShard;
        result.sourceRowCount = state.nextSourceRow;
        result.recordCount = state.metrics.accepted;
        result.rejectedCount = state.metrics.rejected;
        result.shards = List.copyOf(state.shards);
        result.runtime = config.runtimeDescription();
        result.buildStatistics = state.metrics.manifestView(System.nanoTime() - started);
        return result;
    }

    private MoleculeFingerprintIndexBuildResult completedResult(Path manifestPath)
            throws IOException {
        MoleculeFingerprintIndexManifest manifest = MAPPER.readValue(
                manifestPath.toFile(), MoleculeFingerprintIndexManifest.class);
        manifest.validate();
        if (!paths.library().toString().equals(manifest.input)
                || !config.inputs.smilesColumn.equals(manifest.smilesColumn)
                || !config.inputs.idColumn.equals(manifest.idColumn)
                || !paths.modelBundle().toString().equals(manifest.modelBundle)
                || !paths.compactBundle().toString().equals(manifest.compactBundle)
                || config.output.sourceRowsPerShard != manifest.sourceRowsPerShard) {
            throw new IOException("completed index does not match this configuration");
        }
        for (MoleculeFingerprintIndexShard shard : manifest.shards) {
            if (!Files.isRegularFile(paths.outputDirectory().resolve(shard.directory)
                    .resolve(".complete"))) {
                throw new IOException("completed index has a missing shard");
            }
        }
        MoleculeFingerprintIndexBuildMetrics metrics =
                new MoleculeFingerprintIndexBuildMetrics();
        metrics.sourceRows = manifest.sourceRowCount;
        metrics.accepted = manifest.recordCount;
        metrics.rejected = manifest.rejectedCount;
        return new MoleculeFingerprintIndexBuildResult(manifestPath,
                manifest.sourceRowCount, manifest.recordCount,
                manifest.shards.size(), true, metrics);
    }

    private void requireInputFiles() throws IOException {
        if (!Files.isRegularFile(paths.library())) throw new IOException(
                "input molecule table is missing: " + paths.library());
        if (!Files.isDirectory(paths.modelBundle())) throw new IOException(
                "model bundle is missing: " + paths.modelBundle());
        if (!Files.isDirectory(paths.compactBundle())) throw new IOException(
                "compact model bundle is missing: " + paths.compactBundle());
    }

    private static void cleanupPartialDirectories(Path output) throws IOException {
        try (var entries = Files.list(output)) {
            for (Path path : entries.filter(p ->
                    p.getFileName().toString().endsWith(".partial")).toList()) {
                deleteTree(path);
            }
        }
    }

    private static void deleteTree(Path path) throws IOException {
        if (!Files.exists(path)) return;
        try (var entries = Files.walk(path)) {
            for (Path entry : entries.sorted(java.util.Comparator.reverseOrder()).toList()) {
                Files.delete(entry);
            }
        }
    }

    private static void writeState(Path output, BuildState state) throws IOException {
        writeJsonAtomic(output.resolve(STATE_FILE + ".partial"),
                output.resolve(STATE_FILE), state);
    }

    private static void writeJsonAtomic(Path temporary, Path destination, Object value)
            throws IOException {
        MAPPER.writeValue(temporary.toFile(), value);
        moveAtomic(temporary, destination, true);
    }

    private static void moveAtomic(Path source, Path destination, boolean replace)
            throws IOException {
        try {
            if (replace) {
                Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE,
                        StandardCopyOption.REPLACE_EXISTING);
            } else {
                Files.move(source, destination, StandardCopyOption.ATOMIC_MOVE);
            }
        } catch (AtomicMoveNotSupportedException error) {
            throw new IOException("atomic moves are required in the output directory", error);
        }
    }

    private void reportProgress(BuildState state, long started) {
        reportProgress(state.metrics, started);
    }

    private void reportProgress(MoleculeFingerprintIndexBuildMetrics metrics, long started) {
        if (config.runtime.progressIntervalSeconds == 0) return;
        long now = System.nanoTime();
        if (now - metrics.lastProgressNanos < Duration.ofSeconds(
                config.runtime.progressIntervalSeconds).toNanos()) return;
        metrics.lastProgressNanos = now;
        double seconds = (now - started) / 1_000_000_000d;
        System.out.printf(
                "Molecule fingerprints: %,d source rows, %,d accepted, %.1f rows/s, %.1f accepted/s%n",
                metrics.sourceRows, metrics.accepted,
                metrics.sourceRows / Math.max(0.001, seconds),
                metrics.accepted / Math.max(0.001, seconds));
    }

    private static String shardName(int index) {
        return String.format("shard-%05d", index);
    }

    public static final class BuildState {
        public int formatVersion = 1;
        public String status = "IN_PROGRESS";
        public String input;
        public String smilesColumn;
        public String idColumn;
        public String modelBundle;
        public String compactBundle;
        public long sourceRowsPerShard;
        public long nextSourceRow;
        public List<MoleculeFingerprintIndexShard> shards = new ArrayList<>();
        public MoleculeFingerprintIndexBuildMetrics metrics =
                new MoleculeFingerprintIndexBuildMetrics();
    }

    private record Header(int smilesIndex, int idIndex) {}
    private record RawRecord(long sourceRow, String id, String smiles, String rejection) {}
    private record PendingBatch(List<Future<Prepared>> futures, int sourceRows) {}
    private record ShardRun(long sourceRows, long accepted, boolean eof,
                            Map<String, Long> rejections) {}
    private record Prepared(long sourceRow, String id, String smiles,
                            FeaturizationResult features, int heavyAtoms,
                            String rejection, long preparationNanos) {
        boolean accepted() { return rejection == null; }
        static Prepared rejected(RawRecord raw, String reason, long nanos) {
            return new Prepared(raw.sourceRow, raw.id, raw.smiles,
                    null, 0, reason, nanos);
        }
    }
}
