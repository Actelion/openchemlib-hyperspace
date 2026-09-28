package com.idorsia.research.chem.hyperspace3d;

import com.idorsia.research.chem.hyperspace3d.workflow.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;
import static com.idorsia.research.chem.hyperspace3d.workflow.WorkflowFiles.*;

class FingerprintWorkflowTest {
    @TempDir Path temp;
    static final String INPUT = "smiles\tid\nCCO\tduplicate\nc1ccccc1\tduplicate\n\nwrong\nCCN\tlast\n";

    Path config(String extension) throws Exception {
        Path input = temp.resolve("input.tsv" + extension);
        try (OutputStream raw = Files.newOutputStream(input);
             OutputStream out = extension.equals(".gz") ? new GZIPOutputStream(raw)
                     : extension.equals(".bz2") ? new BZip2CompressorOutputStream(raw) : raw) {
            out.write(INPUT.getBytes(java.nio.charset.StandardCharsets.UTF_8));
        }
        Path lib = Files.createDirectories(temp.resolve("lib"));
        Files.writeString(lib.resolve("onnxruntime-test.jar"), "test distribution");
        var c = new FingerprintWorkflowConfig();
        c.library = input.toString(); c.sharedDirectory = temp.resolve("shared").toString();
        c.scratchRoot = temp.resolve("scratch").toString();
        c.applicationJar = Files.writeString(temp.resolve("app.jar"), "application").toString();
        c.libraryDirectory = lib.toString();
        c.modelBundle = Path.of(System.getProperty("hyperspace3d.modelBundle", "model-bundles/deepspace7-v1")).toAbsolutePath().toString();
        c.compactBundle = Path.of(System.getProperty("hyperspace3d.compactBundle", "model-bundles/deepspace7-skelspheres16")).toAbsolutePath().toString();
        c.rowsPerPartition = 2; c.runtime.device = "CPU";
        Path path = temp.resolve("config.json"); json(path, c); return path;
    }

    @Test void prepareResumeAndFinalizationPreserveRowsAndVectors() throws Exception {
        var w = new FingerprintWorkflow(config(".gz"));
        var m = w.prepare();
        assertEquals(5, m.sourceRowCount); assertEquals(3, m.partitions.size());
        assertEquals(m.workflowId, w.prepare().workflowId);
        assertThrows(IOException.class, w::arrayTasks);
        json(w.root.resolve("smoke-success.json"), Map.of("workflowId", m.workflowId));
        assertEquals(3, w.arrayTasks());
        partition(w, m, 0, List.of(0L, 1L));
        partition(w, m, 1, List.of());
        // A partial upload must never count as a completed partition.
        Files.createDirectories(w.root.resolve("partitions/part-000002.upload"));
        assertThrows(IOException.class, () -> FingerprintFinalizer.finish(w));
        assertFalse(Files.exists(w.root.resolve("cache/manifest.json")));
        assertTrue(Files.exists(w.root.resolve("cache/shard-00000/.complete")));
        partition(w, m, 2, List.of(0L));
        Path manifest = FingerprintFinalizer.finish(w);
        assertEquals(manifest, FingerprintFinalizer.finish(w));
        Path cache = manifest.getParent();
        var result = MoleculeFingerprintIndexReader.loadManifest(cache);
        assertEquals(3, result.recordCount); assertEquals(2, result.rejectedCount);
        List<Long> rows = new ArrayList<>(); List<String> ids = new ArrayList<>();
        for (int i = 0; i < result.shards.size(); i++) {
            try (var reader = new MoleculeFingerprintIndexReader(cache, result, i)) {
                for (var record : reader.readBatch(10)) {
                    rows.add(record.sourceRow()); ids.add(record.moleculeId());
                    assertEquals(1f, record.base128()[0]); assertEquals(1f, record.compact16()[0]);
                }
            }
            assertArrayEquals(Files.readAllBytes(w.partitionDirectory(i).resolve("shard-00000/vectors-128.fp16")),
                    Files.readAllBytes(cache.resolve(result.shards.get(i).directory).resolve("vectors-128.fp16")));
        }
        assertEquals(List.of(0L, 1L, 4L), rows);
        assertEquals(List.of("duplicate", "duplicate", "last"), ids);
        // Model identity, not an old absolute pathname, permits a relocated bundle.
        Path moved = temp.resolve("moved-model"); copyTree(Path.of(w.config.modelBundle), moved);
        var base = DeepSpaceModelBundle.load(moved);
        var compact = CompactSkelSpheresModelBundle.load(Path.of(w.config.compactBundle), base.manifest());
        try (var source = new JavaMoleculeFingerprintDataSource(cache)) {
            source.validateCompatibility(base, compact, true);
        }
        result.modelBundleHash = "wrong"; json(manifest, result);
        try (var source = new JavaMoleculeFingerprintDataSource(cache)) {
            assertThrows(IllegalArgumentException.class, () -> source.validateCompatibility(base, compact, true));
        }
    }

    @Test void rejectsChangedInputsModelsAndCorruptPartitions() throws Exception {
        var w = new FingerprintWorkflow(config("")); var m = w.prepare();
        partition(w, m, 0, List.of(0L));
        try (var lock = lock(w.root.resolve("prepare.lock"))) {
            assertThrows(Exception.class, w::prepare);
        }
        Files.writeString(Path.of(w.config.library), INPUT.replace("CCO", "CCC"));
        assertThrows(IOException.class, w::prepare);
        Files.writeString(w.partitionDirectory(0).resolve("shard-00000/strings.bin"), "broken");
        assertThrows(IOException.class, () -> w.verifyPartition(m, m.partitions.get(0), w.partitionDirectory(0)));
        Files.writeString(Path.of(w.config.applicationJar), "changed");
        assertThrows(IOException.class, w::prepare);
    }

    @Test void bz2AndInterruptedPreparationReuseVerifiedChunks() throws Exception {
        var w = new FingerprintWorkflow(config(".bz2")); var first = w.prepare();
        Files.delete(w.root.resolve("workflow.json"));
        assertEquals(first.workflowId, w.prepare().workflowId);
        var edited = JSON.readTree(w.root.resolve("workflow.json").toFile());
        ((com.fasterxml.jackson.databind.node.ObjectNode) edited).put("sourceRowCount", 99);
        json(w.root.resolve("workflow.json"), edited);
        assertThrows(IOException.class, w::load);
    }

    @Test void zeroAcceptedCacheIsNotPublished() throws Exception {
        var w = new FingerprintWorkflow(config("")); var m = w.prepare();
        for (int i = 0; i < m.partitions.size(); i++) partition(w, m, i, List.of());
        assertThrows(IOException.class, () -> FingerprintFinalizer.finish(w));
        assertFalse(Files.exists(w.root.resolve("cache/manifest.json")));
    }

    static void partition(FingerprintWorkflow w, FingerprintWorkflowManifest m, int i, List<Long> rows) throws Exception {
        Path output = w.partitionDirectory(i); Path shardDir = output.resolve("shard-00000");
        float[] base = new float[128], compact = new float[16]; base[0] = 1; compact[0] = 1;
        try (var writer = new MoleculeFingerprintIndexWriter(shardDir)) {
            for (long row : rows) writer.write(row, 3, i == 2 ? "last" : "duplicate", "CCO", base, compact);
        }
        var shard = new MoleculeFingerprintIndexShard();
        shard.directory = "shard-00000"; shard.sourceRowCount = m.partitions.get(i).sourceRowCount();
        shard.recordCount = rows.size(); shard.rejectedCount = shard.sourceRowCount - rows.size();
        shard.rejections.put("TEST_REJECTION", shard.rejectedCount);
        json(shardDir.resolve("shard.json"), shard); Files.createFile(shardDir.resolve(".complete"));
        var index = new MoleculeFingerprintIndexManifest();
        index.input = "partition"; index.smilesColumn = "smiles"; index.idColumn = "id";
        index.modelBundle = w.config.modelBundle; index.compactBundle = w.config.compactBundle;
        index.modelBundleHash = m.modelBundleHash; index.compactBundleHash = m.compactBundleHash;
        index.sourceRowsPerShard = w.config.rowsPerPartition; index.sourceRowCount = shard.sourceRowCount;
        index.recordCount = shard.recordCount; index.rejectedCount = shard.rejectedCount; index.shards = List.of(shard);
        json(output.resolve("manifest.json"), index);
        json(output.resolve("completion.json"), new FingerprintWorkflowManifest.Completion(m.workflowId, i, hashes(output)));
    }
}
