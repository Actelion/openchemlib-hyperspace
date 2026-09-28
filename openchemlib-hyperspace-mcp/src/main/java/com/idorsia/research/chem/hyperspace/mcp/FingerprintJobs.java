package com.idorsia.research.chem.hyperspace.mcp;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/** MCP stays chemistry-free: it inspects bounded text and launches the thin 3D distribution. */
final class FingerprintJobs {
    static Map<String, Object> inspect(Path input) throws IOException {
        if (!Files.isRegularFile(input)) throw new IOException("library is not a regular file");
        try (InputStream raw = Files.newInputStream(input);
             InputStream decoded = input.toString().endsWith(".gz") ? new GZIPInputStream(raw)
                     : input.toString().endsWith(".bz2") ? new BZip2CompressorInputStream(raw, true) : raw;
             Reader reader = new InputStreamReader(decoded, StandardCharsets.UTF_8)) {
            char[] chars = new char[65536];
            int count = 0, n;
            while (count < chars.length && (n = reader.read(chars, count, chars.length - count)) != -1) count += n;
            String text = new String(chars, 0, count);
            String[] lines = text.split("\\R", -1);
            List<String> warnings = new ArrayList<>();
            if (text.isEmpty()) warnings.add("Input is empty");
            if (!lines[0].contains("\t")) warnings.add("Expected a headered tab-delimited table");
            if (lines[0].length() == chars.length) warnings.add("Header exceeds inspection limit");
            return Map.of("input", input.toAbsolutePath().toString(), "columns", List.of(lines[0].split("\t", -1)),
                    "sampleLines", Arrays.stream(lines).skip(1).limit(Math.min(5, Math.max(0, lines.length - 2))).toList(),
                    "bounded", true, "warnings", warnings);
        }
    }

    static Map<String, Object> normalize(Map<String, Object> input, ServerConfig c) throws Exception {
        if (c.fingerprintJar() == null || c.fingerprintLibDirectory() == null || c.modelBundle() == null || c.compactBundle() == null)
            throw new IllegalArgumentException("Configure fingerprintJar, fingerprintLibDirectory, modelBundle and compactBundle first");
        for (Path p : List.of(c.fingerprintLibDirectory(), c.modelBundle(), c.compactBundle()))
            if (!Files.isDirectory(p)) throw new IllegalArgumentException("Required directory missing: " + p);
        Path library = Path.of(JsonFiles.string(input, "input")).toRealPath();
        if (!Files.isRegularFile(library)) throw new IllegalArgumentException("input must be a file");
        Path output = input.containsKey("outputDirectory") ? Path.of(JsonFiles.string(input, "outputDirectory"))
                : c.workspace().resolve("libraries").resolve(UUID.randomUUID().toString());
        output = output.toAbsolutePath().normalize();
        if (output.getParent() == null || output.equals(library)) throw new IllegalArgumentException("invalid output directory");
        Map<String, Object> r = new LinkedHashMap<>();
        r.put("jobKind", "fingerprint"); r.put("input", library.toString()); r.put("outputDirectory", output.toString());
        r.put("smilesColumn", input.getOrDefault("smilesColumn", "smiles"));
        r.put("idColumn", input.getOrDefault("idColumn", "id"));
        if (JsonFiles.string(r, "smilesColumn").equals(JsonFiles.string(r, "idColumn")))
            throw new IllegalArgumentException("SMILES and ID columns must differ");
        String device = Objects.toString(input.getOrDefault("device", "CPU")).toUpperCase(Locale.ROOT);
        if (!Set.of("CPU", "CUDA").contains(device)) throw new IllegalArgumentException("device must be CPU or CUDA");
        r.put("device", device);
        r.put("threads", JsonFiles.integer(input, "threads", c.threads(), 1, 4096));
        r.put("encoderBatchSize", JsonFiles.integer(input, "encoderBatchSize", device.equals("CUDA") ? 512 : 64, 1, 65536));
        r.put("cudaDeviceId", JsonFiles.integer(input, "cudaDeviceId", 0, 0, 1024));
        r.put("sourceRowsPerShard", JsonFiles.integer(input, "sourceRowsPerShard", 1_000_000, 1, 8_000_000));
        r.put("heap", ServerConfig.validHeap((String) input.getOrDefault("heap", c.heap())));
        return r;
    }

    static Map<String, Object> buildConfig(ServerConfig c, Map<String, Object> r) {
        return Map.of("inputs", Map.of("library", r.get("input"), "smilesColumn", r.get("smilesColumn"),
                        "idColumn", r.get("idColumn"), "modelBundle", c.modelBundle().toString(), "compactBundle", c.compactBundle().toString()),
                "output", Map.of("directory", r.get("outputDirectory"), "sourceRowsPerShard", r.get("sourceRowsPerShard"), "resume", true),
                "runtime", Map.of("device", r.get("device"), "cpuWorkers", r.get("threads"),
                        "encoderBatchSize", r.get("encoderBatchSize"), "cudaDeviceId", r.get("cudaDeviceId")));
    }
    static List<String> command(ServerConfig c, Map<String, Object> r, Path request) {
        return List.of(c.javaExecutable(), "-Xmx" + r.get("heap"), "-cp",
                c.fingerprintJar() + File.pathSeparator + c.fingerprintLibDirectory().resolve("*"),
                "com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintIndexBuilderCLI", "--config", request.toString());
    }

    @SuppressWarnings("unchecked")
    static void validateIndex(Path output, Map<String, Object> m) throws Exception {
        if (!"hyperspace-molecule-fingerprint-index".equals(m.get("artifactType"))
                || ((Number) m.getOrDefault("recordCount", 0)).longValue() <= 0)
            throw new IOException("No valid molecule fingerprint index was produced");
        long records = 0, sources = 0;
        List<Map<String, Object>> shards = (List<Map<String, Object>>) m.get("shards");
        if (shards == null) throw new IOException("missing index shards");
        for (Map<String, Object> shard : shards) {
            Path p = output.resolve(JsonFiles.string(shard, "directory")).normalize();
            if (!p.startsWith(output) || p.equals(output)) throw new IOException("invalid shard path");
            for (String name : List.of(".complete", "shard.json", "vectors-128.fp16", "vectors-16.fp16", "rows.bin", "strings.bin"))
                if (!Files.isRegularFile(p.resolve(name))) throw new IOException("Missing fingerprint artifact: " + p.resolve(name));
            long count = ((Number) shard.get("recordCount")).longValue();
            if (count < 0 || Files.size(p.resolve("vectors-128.fp16")) != 24 + Math.multiplyExact(count, 256)
                    || Files.size(p.resolve("vectors-16.fp16")) != 24 + Math.multiplyExact(count, 32)
                    || Files.size(p.resolve("rows.bin")) != 24 + Math.multiplyExact(count, 24))
                throw new IOException("Invalid vector/row file size");
            records += count; sources += ((Number) shard.get("sourceRowCount")).longValue();
        }
        if (records != ((Number) m.get("recordCount")).longValue()
                || sources != ((Number) m.get("sourceRowCount")).longValue()) throw new IOException("inconsistent index totals");
    }
}
