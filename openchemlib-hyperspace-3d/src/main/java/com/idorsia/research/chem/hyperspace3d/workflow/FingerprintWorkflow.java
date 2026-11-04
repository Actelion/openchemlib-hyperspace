package com.idorsia.research.chem.hyperspace3d.workflow;

import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.DigestInputStream;
import java.util.*;
import java.util.zip.*;
import static com.idorsia.research.chem.hyperspace3d.workflow.WorkflowFiles.*;
import static com.idorsia.research.chem.hyperspace3d.workflow.FingerprintWorkflowManifest.*;

/** One process owns each partition; only verified immutable outputs are shared between tasks. */
public final class FingerprintWorkflow {
    public final FingerprintWorkflowConfig config;
    public final Path root;

    public FingerprintWorkflow(Path configPath) throws IOException {
        config = JSON.readValue(configPath.toFile(), FingerprintWorkflowConfig.class);
        config.resolve(configPath);
        root = Path.of(config.sharedDirectory);
    }

    public FingerprintWorkflowManifest prepare() throws Exception {
        Files.createDirectories(root);
        Files.createDirectories(Path.of(config.scratchRoot));
        try (var lock = lock(root.resolve("prepare.lock"))) {
            FingerprintWorkflowManifest m = new FingerprintWorkflowManifest();
            m.configHash = identity(config);
            var base = DeepSpaceModelBundle.load(Path.of(config.modelBundle));
            var compact = CompactSkelSpheresModelBundle.load(Path.of(config.compactBundle), base.manifest());
            m.modelBundleHash = base.bundleHash();
            m.compactBundleHash = compact.bundleHash();
            m.applicationFiles = applicationHashes();
            m.environmentSetupSha256 = environmentHash();
            Path binding = root.resolve("preparation.json");
            String bindingId = identity(List.of(m.configHash, m.modelBundleHash, m.compactBundleHash, m.applicationFiles, m.environmentSetupSha256));
            if (Files.exists(binding) && !JSON.readTree(binding.toFile()).path("identity").asText().equals(bindingId))
                throw new IOException("workflow directory belongs to another configuration/model; use a new directory");
            json(binding, Map.of("identity", bindingId));
            Path scratch = Files.createTempDirectory(Path.of(config.scratchRoot), "hyperspace-prepare-");
            try {
                var digest = digest();
                try (InputStream raw = new DigestInputStream(Files.newInputStream(Path.of(config.library)), digest);
                     BufferedReader reader = table(raw, config.library)) {
                    String header = reader.readLine();
                    validateHeader(header, config.smilesColumn, config.idColumn);
                    String line = reader.readLine();
                    while (line != null) {
                        int index = m.partitions.size();
                        Path local = scratch.resolve("chunk.tsv.gz");
                        long count = 0;
                        try (var out = new BufferedWriter(new OutputStreamWriter(
                                new GZIPOutputStream(Files.newOutputStream(local)), StandardCharsets.UTF_8))) {
                            out.write(header); out.write('\n');
                            while (line != null && count < config.rowsPerPartition) {
                                out.write(line); out.write('\n'); count++;
                                line = reader.readLine();
                            }
                        }
                        String relative = String.format(Locale.ROOT, "inputs/part-%06d.tsv.gz", index);
                        String checksum = hash(local);
                        Path destination = root.resolve(relative);
                        publishFile(local, destination, checksum);
                        m.partitions.add(new Partition(index, relative, m.sourceRowCount, count, checksum));
                        m.sourceRowCount += count;
                        System.out.printf("Prepared partition %d: %,d source rows%n", index, count);
                    }
                }
                m.inputSha256 = HexFormat.of().formatHex(digest.digest());
                if (m.sourceRowCount == 0) throw new IOException("input contains no source rows");
                m.workflowId = identity(m);
                Path manifest = root.resolve("workflow.json");
                if (Files.exists(manifest) && !JSON.readTree(manifest.toFile()).path("workflowId").asText().equals(m.workflowId))
                    throw new IOException("source changed since preparation; use a new workflow directory");
                json(manifest, m);
                return m;
            } finally { deleteTree(scratch); }
        }
    }

    public FingerprintWorkflowManifest load() throws IOException {
        var m = JSON.readValue(root.resolve("workflow.json").toFile(), FingerprintWorkflowManifest.class);
        String storedId = m.workflowId;
        m.workflowId = null;
        if (m.formatVersion != 1 || !Objects.equals(storedId, identity(m))
                || !Objects.equals(m.configHash, identity(config)))
            throw new IOException("workflow manifest/configuration mismatch");
        m.workflowId = storedId;
        long next = 0;
        for (int i = 0; i < m.partitions.size(); i++) {
            var p = m.partitions.get(i);
            if (p.index() != i || p.firstSourceRow() != next || p.sourceRowCount() < 1
                    || p.sourceRowCount() > config.rowsPerPartition)
                throw new IOException("invalid partition interval");
            child(root, p.input());
            next = Math.addExact(next, p.sourceRowCount());
        }
        if (next == 0 || next != m.sourceRowCount) throw new IOException("invalid source row totals");
        return m;
    }

    public void worker(int task, int tasks) throws Exception {
        var m = load();
        if (tasks < 1 || tasks > config.maxArrayTasks || task < 0 || task >= tasks)
            throw new IllegalArgumentException("invalid array task index/count");
        Files.createDirectories(Path.of(config.scratchRoot));
        Path scratch = Files.createTempDirectory(Path.of(config.scratchRoot), "hyperspace-worker-");
        try {
            stageRuntime(scratch, m);
            for (int i = task; i < m.partitions.size(); i += tasks) {
                Partition p = m.partitions.get(i);
                Path destination = partitionDirectory(i);
                try (var lock = lock(root.resolve(String.format(Locale.ROOT, "locks/part-%06d.lock", i)))) {
                    if (Files.exists(destination)) {
                        verifyPartition(m, p, destination);
                        System.out.println("Already complete: partition " + i);
                        continue;
                    }
                    Path work = scratch.resolve("partition");
                    deleteTree(work); Files.createDirectories(work);
                    Path input = work.resolve("input.tsv.gz");
                    Files.copy(child(root, p.input()), input);
                    if (!hash(input).equals(p.sha256())) throw new IOException("input partition checksum mismatch");
                    Path output = work.resolve("index");
                    runBuilder(scratch, input, output);
                    validateIndex(output, m, p.sourceRowCount());
                    json(output.resolve("completion.json"), new Completion(m.workflowId, i, hashes(output)));
                    Files.createDirectories(destination.getParent());
                    Path upload = destination.resolveSibling(destination.getFileName() + ".upload");
                    deleteTree(upload);
                    copyTree(output, upload);
                    verifyPartition(m, p, upload);
                    Files.move(upload, destination, StandardCopyOption.ATOMIC_MOVE);
                    deleteTree(work);
                    System.out.println("Published partition " + i);
                }
            }
        } finally { deleteTree(scratch); }
    }

    public void smoke() throws Exception {
        var m = load();
        try (var lock = lock(root.resolve("smoke.lock"))) {
            Files.deleteIfExists(root.resolve("smoke-success.json"));
            Files.createDirectories(Path.of(config.scratchRoot));
            Path scratch = Files.createTempDirectory(Path.of(config.scratchRoot), "hyperspace-smoke-");
            try {
                stageRuntime(scratch, m);
                Path input = scratch.resolve("smoke.tsv");
                Files.writeString(input, config.smilesColumn + "\t" + config.idColumn
                        + "\nCCCCCC\tsmoke-1\nc1ccccc1\tsmoke-2\nCC(=O)NCC\tsmoke-3\n");
                runBuilder(scratch, input, scratch.resolve("index"));
                var result = validateIndex(scratch.resolve("index"), m, 3);
                if (result.recordCount != 3) throw new IOException("smoke fixture did not encode all molecules");
                json(root.resolve("smoke-success.json"), Map.of("workflowId", m.workflowId,
                        "device", config.runtime.device, "host", java.net.InetAddress.getLocalHost().getHostName(),
                        "javaVersion", System.getProperty("java.version"), "completedAt", java.time.Instant.now().toString()));
            } finally { deleteTree(scratch); }
        }
    }

    public int arrayTasks() throws IOException {
        var m = load();
        Path stamp = root.resolve("smoke-success.json");
        if (!Files.isRegularFile(stamp)
                || !JSON.readTree(stamp.toFile()).path("workflowId").asText().equals(m.workflowId))
            throw new IOException("run a successful matching smoke job before production submission");
        return Math.min(config.maxArrayTasks, m.partitions.size());
    }

    public Path partitionDirectory(int i) { return root.resolve(String.format(Locale.ROOT, "partitions/part-%06d", i)); }

    public MoleculeFingerprintIndexManifest verifyPartition(FingerprintWorkflowManifest m, Partition p, Path directory) throws IOException {
        Completion completion = JSON.readValue(directory.resolve("completion.json").toFile(), Completion.class);
        if (!m.workflowId.equals(completion.workflowId()) || completion.partition() != p.index())
            throw new IOException("partition belongs to another workflow");
        verify(directory, completion.files());
        return validateIndex(directory, m, p.sourceRowCount());
    }

    static MoleculeFingerprintIndexManifest validateIndex(Path directory, FingerprintWorkflowManifest m, long sources) throws IOException {
        var index = MoleculeFingerprintIndexReader.loadManifest(directory);
        if (index.sourceRowCount != sources || !m.modelBundleHash.equals(index.modelBundleHash)
                || !m.compactBundleHash.equals(index.compactBundleHash))
            throw new IOException("partition counts/model identity mismatch");
        for (int i = 0; i < index.shards.size(); i++) {
            child(directory, index.shards.get(i).directory);
            try (var reader = new MoleculeFingerprintIndexReader(directory, index, i)) {
                var shard = index.shards.get(i);
                long previous = shard.firstSourceRow - 1;
                // Sequential fixed-row validation avoids billions of random string-table seeks on shared storage.
                try (var rows = new DataInputStream(new BufferedInputStream(Files.newInputStream(
                        child(directory, shard.directory).resolve("rows.bin")), 1 << 20))) {
                    rows.skipNBytes(24);
                    for (long row = 0; row < reader.recordCount(); row++) {
                        long source = Long.reverseBytes(rows.readLong());
                        rows.skipNBytes(16);
                        if (source <= previous || source >= shard.firstSourceRow + shard.sourceRowCount)
                            throw new IOException("invalid source row in partition");
                        previous = source;
                    }
                }
            }
        }
        return index;
    }

    private Map<String, String> applicationHashes() throws IOException {
        Map<String, String> values = new TreeMap<>();
        values.put("application.jar", hash(Path.of(config.applicationJar)));
        try (var entries = Files.list(Path.of(config.libraryDirectory))) {
            for (Path p : entries.filter(p -> p.getFileName().toString().endsWith(".jar")).sorted().toList())
                values.put("lib/" + p.getFileName(), hash(p));
        }
        if (values.size() < 2) throw new IOException("distribution dependency directory has no jars");
        long ort = values.keySet().stream().filter(k -> k.startsWith("lib/onnxruntime") && k.endsWith(".jar")).count();
        if (ort != 1) throw new IOException("distribution must contain exactly one ONNX Runtime jar; use a clean build");
        return values;
    }

    private void stageRuntime(Path scratch, FingerprintWorkflowManifest m) throws IOException {
        if (!environmentHash().equals(m.environmentSetupSha256)) throw new IOException("environment setup script changed");
        if (!applicationHashes().equals(m.applicationFiles)) throw new IOException("application distribution changed");
        Path runtime = scratch.resolve("runtime");
        Files.createDirectories(runtime.resolve("lib"));
        Files.copy(Path.of(config.applicationJar), runtime.resolve("application.jar"));
        for (String key : m.applicationFiles.keySet()) if (key.startsWith("lib/"))
            Files.copy(Path.of(config.libraryDirectory).resolve(key.substring(4)), child(runtime, key));
        verify(runtime, m.applicationFiles);
        copyTree(Path.of(config.modelBundle), scratch.resolve("model"));
        copyTree(Path.of(config.compactBundle), scratch.resolve("compact"));
        var base = DeepSpaceModelBundle.load(scratch.resolve("model"));
        var compact = CompactSkelSpheresModelBundle.load(scratch.resolve("compact"), base.manifest());
        if (!base.bundleHash().equals(m.modelBundleHash) || !compact.bundleHash().equals(m.compactBundleHash))
            throw new IOException("model bundle changed");
    }

    private String environmentHash() throws IOException {
        return config.environmentSetup == null ? "none" : hash(Path.of(config.environmentSetup));
    }

    private void runBuilder(Path scratch, Path input, Path output) throws Exception {
        Path request = scratch.resolve("build.json");
        json(request, config.build(input, output, scratch.resolve("model"), scratch.resolve("compact")));
        String cp = scratch.resolve("runtime/application.jar") + File.pathSeparator + scratch.resolve("runtime/lib/*");
        Process process = new ProcessBuilder(config.javaExecutable, "-Xmx" + config.heap, "-cp", cp,
                "com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintIndexBuilderCLI", "--config", request.toString())
                .directory(scratch.toFile()).redirectInput(new File("/dev/null")).redirectErrorStream(true).start();
        try {
            try (var log = process.getInputStream()) { log.transferTo(System.out); }
            if (process.waitFor() != 0) throw new IOException("fingerprint subprocess failed; see job log");
        } finally {
            if (process.isAlive()) { process.destroyForcibly(); process.waitFor(); }
        }
    }

    static BufferedReader table(InputStream raw, String name) throws IOException {
        InputStream decoded = name.endsWith(".gz") ? new GZIPInputStream(raw)
                : name.endsWith(".bz2") ? new BZip2CompressorInputStream(raw, true) : raw;
        return new BufferedReader(new InputStreamReader(decoded, StandardCharsets.UTF_8), 1 << 20);
    }
    static void validateHeader(String header, String smiles, String id) throws IOException {
        if (header == null) throw new IOException("empty input");
        List<String> columns = Arrays.asList(header.split("\t", -1));
        if (Collections.frequency(columns, smiles) != 1 || Collections.frequency(columns, id) != 1)
            throw new IOException("TSV must have unique configured SMILES and ID columns");
    }
    static void publishFile(Path local, Path destination, String checksum) throws IOException {
        Files.createDirectories(destination.getParent());
        if (Files.exists(destination)) {
            if (!hash(destination).equals(checksum)) throw new IOException("existing chunk differs; use a new workflow directory");
            return;
        }
        Path temporary = destination.resolveSibling(destination.getFileName() + ".upload");
        Files.copy(local, temporary, StandardCopyOption.REPLACE_EXISTING);
        if (!hash(temporary).equals(checksum)) throw new IOException("chunk upload checksum mismatch");
        Files.move(temporary, destination, StandardCopyOption.ATOMIC_MOVE);
    }
}
