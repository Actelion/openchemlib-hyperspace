package com.idorsia.research.chem.hyperspace3d;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintRecord;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.MoleculeFingerprintBatch;
import com.idorsia.research.chem.hyperspace3d.model.MoleculeFingerprintBatchEncoder;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class MoleculeFingerprintIndexBuilderTest {
    private static final String TABLE = """
            id	smiles	extra
            mol-a	CCCCCCC	a
            mol-b	CC	b
            mol-c	CCCCCC	c
            mol-d		d
            mol-e	c1ccccc1	e
            """;

    @TempDir Path temporary;

    @Test void buildsDeterministicFlatShardsFromPlainGzipAndBzip2() throws Exception {
        for (String extension : List.of(".tsv", ".tsv.gz", ".tsv.bz2")) {
            Path input = temporary.resolve("molecules" + extension);
            writeTable(input);
            Path output = temporary.resolve("index-" + extension.replace('.', '-'));
            MoleculeFingerprintIndexBuildConfig config = config(input, output);
            var result = new MoleculeFingerprintIndexBuilder(
                    fakeEncoder(), config, paths(config, input, output)).build();

            assertEquals(5, result.sourceRowCount());
            assertEquals(3, result.recordCount());
            assertEquals(3, result.shardCount());
            MoleculeFingerprintIndexManifest manifest =
                    MoleculeFingerprintIndexReader.loadManifest(output);
            assertEquals(2, manifest.rejectedCount);
            assertEquals(1, manifest.shards.get(0).recordCount);
            assertEquals(1, manifest.shards.get(1).recordCount);
            assertEquals(1, manifest.shards.get(2).recordCount);
            assertEquals(2, manifest.buildStatistics.get("rejected"));

            List<MoleculeFingerprintRecord> records = readAll(output, manifest);
            assertEquals(List.of(0L, 2L, 4L),
                    records.stream().map(MoleculeFingerprintRecord::sourceRow).toList());
            assertEquals(List.of("mol-a", "mol-c", "mol-e"),
                    records.stream().map(MoleculeFingerprintRecord::moleculeId).toList());
            assertEquals(List.of("CCCCCCC", "CCCCCC", "c1ccccc1"),
                    records.stream().map(MoleculeFingerprintRecord::smiles).toList());
            assertEquals(7f, records.get(0).base128()[0]);
            assertEquals(0.7f, records.get(0).compact16()[0], 0.001f);
            try (MoleculeFingerprintIndexReader compactOnly =
                         new MoleculeFingerprintIndexReader(output, manifest, 0,
                                 MoleculeFingerprintIndexReader.VectorColumns.COMPACT_16)) {
                assertEquals(0.7f, compactOnly.compact16(0)[0], 0.001f);
                assertThrows(IllegalStateException.class,
                        () -> compactOnly.base128(0));
                assertThrows(IllegalStateException.class,
                        () -> compactOnly.read(0));
            }

            var repeated = new MoleculeFingerprintIndexBuilder(
                    fakeEncoder(), config, paths(config, input, output)).build();
            assertTrue(repeated.resumed());
            assertEquals(3, repeated.recordCount());
        }
    }

    @Test void resumesAfterAnEncoderFailureWithoutReencodingCommittedShard()
            throws Exception {
        Path input = temporary.resolve("resume.tsv");
        writeTable(input);
        Path resumedOutput = temporary.resolve("resume-index");
        MoleculeFingerprintIndexBuildConfig config = config(input, resumedOutput);
        AtomicInteger calls = new AtomicInteger();
        MoleculeFingerprintBatchEncoder failOnSecondBatch = batch -> {
            if (calls.incrementAndGet() == 2) throw new IllegalStateException("simulated stop");
            return encode(batch);
        };
        assertThrows(IllegalStateException.class, () ->
                new MoleculeFingerprintIndexBuilder(failOnSecondBatch, config,
                        paths(config, input, resumedOutput)).build());
        assertTrue(Files.isRegularFile(
                resumedOutput.resolve("shard-00000/.complete")));
        assertTrue(Files.isDirectory(
                resumedOutput.resolve("shard-00001.partial")));

        var resumed = new MoleculeFingerprintIndexBuilder(fakeEncoder(), config,
                paths(config, input, resumedOutput)).build();
        assertTrue(resumed.resumed());
        assertEquals(3, resumed.recordCount());
        assertFalse(Files.exists(resumedOutput.resolve("shard-00001.partial")));

        Path cleanOutput = temporary.resolve("clean-index");
        MoleculeFingerprintIndexBuildConfig cleanConfig = config(input, cleanOutput);
        new MoleculeFingerprintIndexBuilder(fakeEncoder(), cleanConfig,
                paths(cleanConfig, input, cleanOutput)).build();
        for (int shard = 0; shard < 3; shard++) {
            for (String name : List.of("vectors-128.fp16", "vectors-16.fp16",
                    "rows.bin", "strings.bin")) {
                assertArrayEquals(
                        Files.readAllBytes(resumedOutput.resolve(
                                "shard-%05d/%s".formatted(shard, name))),
                        Files.readAllBytes(cleanOutput.resolve(
                                "shard-%05d/%s".formatted(shard, name))));
            }
        }
    }

    @Test void recoversShardCommittedBeforeStateUpdate() throws Exception {
        Path input = temporary.resolve("recovery.tsv");
        writeTable(input);
        Path output = temporary.resolve("recovery-index");
        MoleculeFingerprintIndexBuildConfig config = config(input, output);
        var paths = paths(config, input, output);
        new MoleculeFingerprintIndexBuilder(fakeEncoder(), config, paths).build();

        ObjectMapper mapper = new ObjectMapper();
        Path statePath = output.resolve("build-state.json");
        MoleculeFingerprintIndexBuilder.BuildState state = mapper.readValue(
                statePath.toFile(), MoleculeFingerprintIndexBuilder.BuildState.class);
        var removed = state.shards.remove(state.shards.size() - 1);
        state.nextSourceRow -= removed.sourceRowCount;
        state.metrics.sourceRows -= removed.sourceRowCount;
        state.metrics.accepted -= removed.recordCount;
        state.metrics.rejected -= removed.rejectedCount;
        removed.rejections.forEach((reason, count) ->
                state.metrics.rejections.computeIfPresent(reason, (key, old) -> old - count));
        mapper.writeValue(statePath.toFile(), state);
        Files.delete(output.resolve("manifest.json"));

        var recovered = new MoleculeFingerprintIndexBuilder(
                fakeEncoder(), config, paths).build();
        assertTrue(recovered.resumed());
        assertEquals(5, recovered.sourceRowCount());
        assertEquals(3, recovered.recordCount());
        assertEquals(3, recovered.shardCount());
    }

    @Test void rejectsChangedStructuralResumeSettings() throws Exception {
        Path input = temporary.resolve("mismatch.tsv");
        writeTable(input);
        Path output = temporary.resolve("mismatch-index");
        MoleculeFingerprintIndexBuildConfig config = config(input, output);
        AtomicInteger calls = new AtomicInteger();
        assertThrows(IllegalStateException.class, () ->
                new MoleculeFingerprintIndexBuilder(batch -> {
                    if (calls.incrementAndGet() == 2) throw new IllegalStateException("stop");
                    return encode(batch);
                }, config, paths(config, input, output)).build());

        MoleculeFingerprintIndexBuildConfig changed = config(input, output);
        changed.output.sourceRowsPerShard = 3;
        assertThrows(IllegalStateException.class, () ->
                new MoleculeFingerprintIndexBuilder(fakeEncoder(), changed,
                        paths(changed, input, output)).build());
    }

    private static List<MoleculeFingerprintRecord> readAll(Path directory,
            MoleculeFingerprintIndexManifest manifest) throws Exception {
        List<MoleculeFingerprintRecord> result = new ArrayList<>();
        for (int shard = 0; shard < manifest.shards.size(); shard++) {
            try (MoleculeFingerprintIndexReader reader =
                         new MoleculeFingerprintIndexReader(directory, manifest, shard)) {
                result.addAll(reader.readBatch(100));
            }
        }
        return result;
    }

    private MoleculeFingerprintIndexBuildConfig config(Path input, Path output)
            throws IOException {
        Path model = temporary.resolve("model");
        Path compact = temporary.resolve("compact");
        Files.createDirectories(model);
        Files.createDirectories(compact);
        MoleculeFingerprintIndexBuildConfig config =
                new MoleculeFingerprintIndexBuildConfig();
        config.inputs.library = input.toString();
        config.inputs.modelBundle = model.toString();
        config.inputs.compactBundle = compact.toString();
        config.output.directory = output.toString();
        config.output.sourceRowsPerShard = 2;
        config.runtime.device = "CPU";
        config.runtime.cpuWorkers = 3;
        config.runtime.encoderBatchSize = 2;
        config.runtime.queueCapacity = 3;
        config.runtime.progressIntervalSeconds = 0;
        config.validate();
        return config;
    }

    private static MoleculeFingerprintIndexBuildConfig.ResolvedPaths paths(
            MoleculeFingerprintIndexBuildConfig config, Path input, Path output) {
        return new MoleculeFingerprintIndexBuildConfig.ResolvedPaths(
                input.toAbsolutePath().normalize(),
                Path.of(config.inputs.modelBundle).toAbsolutePath().normalize(),
                Path.of(config.inputs.compactBundle).toAbsolutePath().normalize(),
                output.toAbsolutePath().normalize());
    }

    private static MoleculeFingerprintBatchEncoder fakeEncoder() {
        return MoleculeFingerprintIndexBuilderTest::encode;
    }

    private static MoleculeFingerprintBatch encode(
            com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch batch) {
        float[][] base = new float[batch.batchSize()][128];
        float[][] compact = new float[batch.batchSize()][16];
        for (int row = 0; row < batch.batchSize(); row++) {
            int atoms = 0;
            for (int atom = 0; atom < 32; atom++) {
                if (batch.atomMask()[row * 32 + atom]) atoms++;
            }
            base[row][0] = atoms;
            base[row][1] = row;
            compact[row][0] = atoms / 10f;
            compact[row][1] = row;
        }
        return new MoleculeFingerprintBatch(base, compact);
    }

    private static void writeTable(Path path) throws IOException {
        try (OutputStream output = compressedOutput(path)) {
            output.write(TABLE.getBytes(StandardCharsets.UTF_8));
        }
    }

    private static OutputStream compressedOutput(Path path) throws IOException {
        OutputStream raw = Files.newOutputStream(path);
        if (path.toString().endsWith(".gz")) return new GZIPOutputStream(raw);
        if (path.toString().endsWith(".bz2")) return new BZip2CompressorOutputStream(raw);
        return raw;
    }
}
