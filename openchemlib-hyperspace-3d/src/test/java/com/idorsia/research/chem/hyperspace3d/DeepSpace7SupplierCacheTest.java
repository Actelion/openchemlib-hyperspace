package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.cli.MoleculeSimilaritySearchCLI;
import com.idorsia.research.chem.hyperspace3d.cli.MoleculeSimilaritySearchConfig;
import com.idorsia.research.chem.hyperspace3d.index.DeepSpace7SupplierCacheDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSources;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresManifest;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;
import com.idorsia.research.chem.hyperspace3d.model.EmbeddingComparator;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculePheSAScreener;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresScreener;
import com.idorsia.research.chem.hyperspace3d.screening.ScreeningObjective;
import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class DeepSpace7SupplierCacheTest {
    private static final String FOUNDATION = "713ec867036277d1197484e6124b876329d23a48cc7792d4fad57e6cb2c9d3e0";
    private static final String PREDICTOR = "29b2d7770ae3740f7b58492e66d88c7b2743ffa9735f7a74543ad06c7182f2d2";
    private static final String PROJECTION = "f58db63b723d00d77cfd4ebfd77bcfd30cdfcc6833cb1c0b56895a509b8392f3";
    @TempDir Path temporary;

    @Test void readsRawFp16ColumnsAndDelayedGzipMetadata() throws Exception {
        fixture();
        try (var source = MoleculeFingerprintDataSources.open(temporary)) {
            assertEquals(DeepSpace7SupplierCacheDataSource.ARTIFACT_TYPE, source.artifactType());
            assertEquals(4, source.recordCount());
            assertEquals(2, source.shardCount());
            try (var reader = source.openShard(0, MoleculeFingerprintColumn.BASE_128)) {
                var batch = reader.readBatch(1);
                assertEquals(1, batch.size());
                assertEquals(128, batch.dimension());
                assertEquals(0.1015625f, batch.values()[0], 0.002f);
                assertEquals(0.8984375f, reader.readVector(1)[0], 0.002f);
            }
            var metadata = source.resolveMetadata(List.of(
                    new com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference(1, 1),
                    new com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference(0, 0)));
            assertEquals("mol-d", metadata.get(new com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference(1, 1)).moleculeId());
            assertEquals("C", metadata.get(new com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference(0, 0)).canonicalSmiles());
        }
    }

    @Test void rejectsIncorrectVectorSizeWithoutHashingPayloads() throws Exception {
        fixture();
        Files.write(temporary.resolve("shard_000001/fingerprints_128d_fp16.bin"), new byte[1]);
        assertThrows(IllegalArgumentException.class,
                () -> new DeepSpace7SupplierCacheDataSource(temporary));
    }

    @Test void exhaustive3dAndCompactCascadeMatchExpectedGlobalRanking() throws Exception {
        fixture();
        DeepSpaceModelManifest manifest = modelManifest();
        EmbeddingComparator comparator = new EmbeddingComparator() {
            @Override public float[][] compare(float[] query, float[][] candidates) {
                float[][] scores = new float[candidates.length][6];
                for (int row = 0; row < candidates.length; row++) scores[row][3] = candidates[row][0];
                return scores;
            }
        };
        try (var source = MoleculeFingerprintDataSources.open(temporary)) {
            var phesa = new MoleculePheSAScreener(comparator, manifest,
                    ScreeningObjective.direct("phesa_total"));
            var exhaustive = phesa.screen(source, new float[128], 2, 1);
            assertEquals(List.of("mol-b", "mol-c"), exhaustive.hits().stream()
                    .map(hit -> hit.molecule().moleculeId()).toList());

            float[] compactQuery = new float[16];
            compactQuery[0] = 1;
            var compact = new MoleculeSkelSpheresScreener(
                    new CompactSkelSpheresScorer(compactManifest())).scanUnresolved(
                            source, compactQuery, 3, 2);
            var cascade = phesa.rerank(source, new float[128], compact.hits(), 2, 2);
            assertEquals(List.of("mol-b", "mol-c"), cascade.hits().stream()
                    .map(hit -> hit.molecule().moleculeId()).toList());
            assertEquals(List.of(2, 3), cascade.hits().stream()
                    .map(hit -> hit.compactRank()).toList());
        }
    }


    @Test void unifiedCliRuns3dAgainstPythonCacheAndCommittedOnnx() throws Exception {
        fixture();
        var config = config();
        config.inputs.index = temporary.toString();
        config.inputs.modelBundle = Path.of("model-bundles/deepspace7-v1")
                .toAbsolutePath().normalize().toString();
        config.query.structure = "CCOc1ccccc1";
        config.screening.mode = "phesa_3d";
        config.screening.resultTopK = 2;
        config.output.reportTopK = 2;
        config.runtime.scanBatchSize = 2;
        config.runtime.comparatorBatchSize = 2;
        config.output.hitsTsv = temporary.resolve("hits.tsv").toString();
        config.output.hitsSdf = temporary.resolve("hits.sdf").toString();
        config.output.summaryMarkdown = temporary.resolve("summary.md").toString();
        config.output.runManifest = temporary.resolve("run.json").toString();
        Path configPath = temporary.resolve("search.json");
        new ObjectMapper().writeValue(configPath.toFile(), config);

        MoleculeSimilaritySearchCLI.main(new String[]{"--config", configPath.toString()});

        assertEquals(3, Files.readAllLines(temporary.resolve("hits.tsv")).size());
        assertEquals("hyperspace-molecule-similarity-search-run",
                new ObjectMapper().readTree(temporary.resolve("run.json").toFile())
                        .path("artifactType").asText());
    }

    @Test void rejectsMismatchedPrimaryModelCheckpoint() throws Exception {
        fixture();
        DeepSpaceModelManifest manifest = modelManifest();
        manifest.predictorCheckpointSha256 = "0".repeat(64);
        var model = new DeepSpaceModelBundle(Path.of("."), manifest, "test");
        try (var source = MoleculeFingerprintDataSources.open(temporary)) {
            assertThrows(IllegalArgumentException.class,
                    () -> source.validateCompatibility(model, null, false));
        }
    }
    @Test void unifiedConfigValidatesModeSpecificInputs() {
        var config = config();
        config.screening.mode = "phesa_3d";
        config.inputs.compactBundle = null;
        config.validate();
        config.screening.mode = "skelspheres_then_phesa";
        assertThrows(IllegalArgumentException.class, config::validate);
        config.inputs.compactBundle = "compact";
        config.screening.objective.target = "ffp_similarity";
        assertThrows(IllegalArgumentException.class, config::validate);
    }

    private void fixture() throws Exception {
        float[] base = {0.1f, 0.9f, 0.8f, 0.2f};
        float[] compact = {1.0f, 0.9f, 0.8f, 0.7f};
        writeShard(0, 0, new float[]{base[0], base[1]},
                new float[]{compact[0], compact[1]}, List.of("mol-a", "mol-b"));
        writeShard(1, 2, new float[]{base[2], base[3]},
                new float[]{compact[2], compact[3]}, List.of("mol-c", "mol-d"));
        Map<String, Object> root = new LinkedHashMap<>();
        root.put("artifact_type", DeepSpace7SupplierCacheDataSource.ARTIFACT_TYPE);
        root.put("complete_source", true);
        root.put("fingerprints", Map.of(
                "128d", Map.of("bytes_per_molecule", 256, "dtype", "float16"),
                "16d", Map.of("bytes_per_molecule", 32, "dtype", "float16")));
        root.put("provenance", provenance());
        root.put("source_rows", 4);
        root.put("valid_molecules", 4);
        root.put("rejected_molecules", 0);
        root.put("shards", 2);
        new ObjectMapper().writeValue(temporary.resolve("manifest.json").toFile(), root);
    }

    private void writeShard(int index, long sourceStart, float[] baseFirst,
            float[] compactFirst, List<String> ids) throws Exception {
        Path directory = temporary.resolve(String.format("shard_%06d", index));
        Files.createDirectories(directory);
        Files.writeString(directory.resolve(".complete"), "complete\n");
        writeVectors(directory.resolve("fingerprints_128d_fp16.bin"), baseFirst, 128);
        writeVectors(directory.resolve("fingerprints_16d_fp16.bin"), compactFirst, 16);
        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                new GZIPOutputStream(Files.newOutputStream(directory.resolve("records.tsv.gz")))))) {
            writer.write("local_index\tsource_row\tmolecule_id\tsmiles\tcanonical_smiles\tnum_atoms\n");
            for (int row = 0; row < ids.size(); row++) {
                writer.write(row + "\t" + (sourceStart + row) + "\t" + ids.get(row)
                        + "\tC\tC\t1\n");
            }
        }
        Map<String, Object> manifest = new LinkedHashMap<>();
        manifest.put("artifact_type", "deepspace7_streamed_supplier_fingerprint_shard_v1");
        manifest.put("shard_index", index);
        manifest.put("source_row_start", sourceStart);
        manifest.put("source_row_stop", sourceStart + ids.size());
        manifest.put("source_rows", ids.size());
        manifest.put("valid_molecules", ids.size());
        manifest.put("rejected_molecules", 0);
        manifest.put("provenance", provenance());
        manifest.put("fingerprints", Map.of(
                "128d", vectorManifest("fingerprints_128d_fp16.bin", ids.size(), 128),
                "16d", vectorManifest("fingerprints_16d_fp16.bin", ids.size(), 16)));
        manifest.put("records", Map.of("path", "records.tsv.gz", "sha256", "d".repeat(64)));
        new ObjectMapper().writeValue(directory.resolve("manifest.json").toFile(), manifest);
    }

    private static Map<String, Object> vectorManifest(String path, int records, int dimension) {
        return Map.of("dtype", "float16", "path", path,
                "shape", List.of(records, dimension), "sha256", "e".repeat(64));
    }

    private static void writeVectors(Path path, float[] firstValues, int dimension)
            throws Exception {
        ByteBuffer data = ByteBuffer.allocate(firstValues.length * dimension * 2)
                .order(ByteOrder.LITTLE_ENDIAN);
        for (float first : firstValues) {
            data.putShort(Float.floatToFloat16(first));
            for (int column = 1; column < dimension; column++) data.putShort((short) 0);
        }
        Files.write(path, data.array());
    }

    private static Map<String, Object> provenance() {
        return Map.of("integrated_checkpoint_sha256", FOUNDATION,
                "predictor_checkpoint_sha256", PREDICTOR,
                "projection_checkpoint_sha256", PROJECTION);
    }

    private static DeepSpaceModelManifest modelManifest() {
        var value = new DeepSpaceModelManifest();
        value.availableTargets = List.of("ffp_similarity", "skelspheres_similarity",
                "flexophore_similarity", "phesa_total", "phesa_pharmacophore", "phesa_shape");
        value.foundationCheckpointSha256 = FOUNDATION;
        value.predictorCheckpointSha256 = PREDICTOR;
        return value;
    }

    private static CompactSkelSpheresManifest compactManifest() {
        var value = new CompactSkelSpheresManifest();
        value.artifactType = "deepspace7-skelspheres16-onnx-bundle";
        value.formatVersion = 1;
        value.modelVersion = "test";
        value.architecture = "mlp_128_128_16_l2";
        value.target = "skelspheres_similarity";
        value.inputEmbeddingDim = 128;
        value.hiddenDim = 128;
        value.outputEmbeddingDim = 16;
        value.l2Normalized = true;
        value.canonicalSeed = 17;
        value.bestEpoch = 1;
        value.recommendedStorageDtype = "float16";
        value.calibration = new CompactSkelSpheresManifest.Calibration();
        value.calibration.type = "sigmoid_affine_cosine";
        value.calibration.scale = 1;
        value.calibration.bias = 0;
        value.sourceFoundationCheckpointSha256 = FOUNDATION;
        value.sourcePredictorCheckpointSha256 = PREDICTOR;
        value.projectionCheckpointSha256 = PROJECTION;
        value.projectionSha256 = "f".repeat(64);
        return value;
    }

    private static MoleculeSimilaritySearchConfig config() {
        var value = new MoleculeSimilaritySearchConfig();
        value.inputs.index = "index";
        value.inputs.modelBundle = "model";
        value.query.structure = "C";
        value.query.identifier = "query";
        value.output.hitsTsv = "hits.tsv";
        value.output.hitsSdf = "hits.sdf";
        value.output.summaryMarkdown = "summary.md";
        value.output.runManifest = "run.json";
        return value;
    }
}
