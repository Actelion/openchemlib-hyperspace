package com.idorsia.research.chem.hyperspace3d.workflow;

import com.idorsia.research.chem.hyperspace3d.index.*;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.*;
import java.util.*;
import static com.idorsia.research.chem.hyperspace3d.workflow.WorkflowFiles.*;
import static com.idorsia.research.chem.hyperspace3d.workflow.FingerprintWorkflowManifest.*;

/** Retains vector bytes unchanged; adjusts only the fixed-width source-row metadata. */
public final class FingerprintFinalizer {
    private FingerprintFinalizer() {}

    public static Path finish(FingerprintWorkflow workflow) throws Exception {
        var m = workflow.load();
        Path output = workflow.root.resolve("cache");
        try (var lock = lock(workflow.root.resolve("finalize.lock"))) {
            Files.createDirectories(output);
            Path owner = output.resolve("workflow-id.json");
            if (Files.exists(owner)) {
                if (!JSON.readTree(owner.toFile()).path("workflowId").asText().equals(m.workflowId))
                    throw new IOException("cache belongs to another workflow");
            } else {
                try (var entries = Files.list(output)) {
                    if (entries.findAny().isPresent()) throw new IOException("unrecognized nonempty cache directory");
                }
                json(owner, Map.of("workflowId", m.workflowId));
            }
            var result = new MoleculeFingerprintIndexManifest();
            result.input = workflow.config.library;
            result.smilesColumn = workflow.config.smilesColumn;
            result.idColumn = workflow.config.idColumn;
            result.modelBundle = workflow.config.modelBundle;
            result.compactBundle = workflow.config.compactBundle;
            result.modelBundleHash = m.modelBundleHash;
            result.compactBundleHash = m.compactBundleHash;
            result.sourceRowsPerShard = workflow.config.rowsPerPartition;
            result.shards = new ArrayList<>();
            Map<String, Long> rejections = new TreeMap<>();
            for (Partition p : m.partitions) {
                Path source = workflow.partitionDirectory(p.index());
                var index = workflow.verifyPartition(m, p, source);
                for (var original : index.shards) {
                    int number = result.shards.size();
                    var shard = JSON.convertValue(original, MoleculeFingerprintIndexShard.class);
                    shard.shardIndex = number;
                    shard.directory = String.format(Locale.ROOT, "shard-%05d", number);
                    shard.firstSourceRow = Math.addExact(p.firstSourceRow(), original.firstSourceRow);
                    Path target = output.resolve(shard.directory);
                    if (Files.exists(target)) {
                        var c = JSON.readValue(target.resolve("completion.json").toFile(), Completion.class);
                        if (!m.workflowId.equals(c.workflowId()) || c.partition() != number)
                            throw new IOException("final shard identity mismatch");
                        verify(target, c.files());
                        if (!identity(JSON.readValue(target.resolve("shard.json").toFile(), MoleculeFingerprintIndexShard.class)).equals(identity(shard)))
                            throw new IOException("final shard metadata mismatch");
                    } else {
                        Path partial = output.resolve(shard.directory + ".partial");
                        deleteTree(partial); Files.createDirectories(partial);
                        Path from = child(source, original.directory);
                        var sourceCompletion = JSON.readValue(source.resolve("completion.json").toFile(), Completion.class);
                        Map<String, String> checksums = new TreeMap<>();
                        for (String name : List.of("vectors-128.fp16", "vectors-16.fp16", "strings.bin")) {
                            boolean linked = linkOrCopy(from.resolve(name), partial.resolve(name));
                            String expected = sourceCompletion.files().get(original.directory + "/" + name);
                            // The source was verified above. A hard link needs no second payload scan.
                            if (!linked && !hash(partial.resolve(name)).equals(expected))
                                throw new IOException("final payload copy checksum mismatch");
                            checksums.put(name, expected);
                        }
                        offsetRows(from.resolve("rows.bin"), partial.resolve("rows.bin"), p.firstSourceRow());
                        json(partial.resolve("shard.json"), shard);
                        Files.createFile(partial.resolve(".complete"));
                        for (String name : List.of("rows.bin", "shard.json", ".complete"))
                            checksums.put(name, hash(partial.resolve(name)));
                        json(partial.resolve("completion.json"), new Completion(m.workflowId, number, checksums));
                        Files.move(partial, target, StandardCopyOption.ATOMIC_MOVE);
                    }
                    result.shards.add(shard);
                    result.sourceRowCount += shard.sourceRowCount;
                    result.recordCount += shard.recordCount;
                    result.rejectedCount += shard.rejectedCount;
                    shard.rejections.forEach((k, v) -> rejections.merge(k, v, Long::sum));
                }
                System.out.println("Finalized partition " + p.index());
            }
            if (result.recordCount == 0) throw new IOException("no accepted molecules; no searchable cache was published");
            result.runtime = Map.of("workflowId", m.workflowId, "partitionCount", m.partitions.size());
            result.buildStatistics = Map.of("rejections", rejections, "sourceRows", result.sourceRowCount,
                    "accepted", result.recordCount, "rejected", result.rejectedCount);
            result.validate();
            for (int i = 0; i < result.shards.size(); i++)
                try (var reader = new MoleculeFingerprintIndexReader(output, result, i)) {
                    if (reader.recordCount() > 0) {
                        reader.readMetadata(0);
                        reader.readMetadata(reader.recordCount() - 1);
                    }
                }
            json(output.resolve("manifest.json"), result);
            return output.resolve("manifest.json");
        }
    }

    static boolean linkOrCopy(Path from, Path to) throws IOException {
        try { Files.createLink(to, from); return true; }
        catch (UnsupportedOperationException | IOException e) {
            System.err.println("Hard link unavailable; copying immutable payload (additional storage needed): " + from);
            Files.copy(from, to);
            return false;
        }
    }

    static void offsetRows(Path from, Path to, long offset) throws IOException {
        try (var in = new BufferedInputStream(Files.newInputStream(from));
             var out = new BufferedOutputStream(Files.newOutputStream(to))) {
            byte[] header = in.readNBytes(24);
            if (header.length != 24) throw new IOException("truncated row header");
            out.write(header);
            byte[] rows = new byte[24 * 32768];
            ByteBuffer values = ByteBuffer.wrap(rows).order(ByteOrder.LITTLE_ENDIAN);
            while (true) {
                int n = in.readNBytes(rows, 0, rows.length);
                if (n == 0) break;
                if (n % 24 != 0) throw new IOException("truncated row record");
                for (int at = 0; at < n; at += 24)
                    values.putLong(at, Math.addExact(values.getLong(at), offset));
                out.write(rows, 0, n);
            }
        }
    }
}
