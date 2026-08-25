package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace3d.cli.MoleculeSkelSpheresSearchConfig;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexShard;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexWriter;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresManifest;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeExactSkelSpheresReranker;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresHit;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresScreener;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MoleculeSkelSpheresSearchTest {
    @TempDir Path temporary;

    @Test void compactOnlyScanRanksAcrossShardsAndResolvesMetadata() throws Exception {
        writeShard(0, new Row("mol-a", "c1ccccc1", vector(1f, 0f)),
                new Row("mol-b", "CCCCCC", vector(0.8f, 0.6f)));
        writeShard(1, new Row("mol-c", "c1ccncc1", vector(0.9f, 0.4359f)));
        MoleculeFingerprintIndexManifest manifest = manifest();
        new com.fasterxml.jackson.databind.ObjectMapper().writeValue(
                temporary.resolve("manifest.json").toFile(), manifest);
        Files.delete(temporary.resolve("shard-00000/vectors-128.fp16"));
        Files.delete(temporary.resolve("shard-00001/vectors-128.fp16"));

        var loaded = MoleculeFingerprintIndexReader.loadManifest(temporary);
        var result = new MoleculeSkelSpheresScreener(
                new CompactSkelSpheresScorer(compactManifest())).screen(
                        temporary, loaded, vector(1f, 0f), 2, 1);
        assertEquals(3, result.recordsScanned());
        assertEquals(List.of("mol-a", "mol-c"), result.hits().stream()
                .map(hit -> hit.molecule().moleculeId()).toList());
        assertEquals(0, result.hits().get(0).reference().shardIndex());
        assertEquals(1, result.hits().get(1).reference().shardIndex());
        assertEquals("c1ccncc1", result.hits().get(1).molecule().smiles());
    }

    @Test void exactRerankerUsesOclBinarySkelSpheres() throws Exception {
        StereoMolecule query = new StereoMolecule();
        new SmilesParser().parse(query, "c1ccccc1");
        var benzene = hit("benzene", "c1ccccc1", 0.5);
        var cyclohexane = hit("cyclohexane", "C1CCCCC1", 0.6);
        var result = new MoleculeExactSkelSpheresReranker().rerank(
                query, List.of(cyclohexane, benzene));
        assertEquals("benzene", result.get(0).molecule().moleculeId());
        assertEquals(1.0, result.get(0).exactSimilarity(), 1.0e-12);
        assertTrue(Double.isFinite(result.get(1).exactSimilarity()));
    }

    @Test void configRejectsInconsistentReportLimit() {
        var config = new MoleculeSkelSpheresSearchConfig();
        config.inputs.index = "index";
        config.inputs.modelBundle = "model";
        config.inputs.compactBundle = "compact";
        config.query.structure = "CC";
        config.query.identifier = "query";
        config.output.hitsTsv = "hits.tsv";
        config.output.hitsSdf = "hits.sdf";
        config.output.summaryMarkdown = "summary.md";
        config.output.runManifest = "run.json";
        config.output.learnedTopK = 10;
        config.output.reportTopK = 11;
        assertThrows(IllegalArgumentException.class, config::validate);
    }

    private void writeShard(int shard, Row... rows) throws Exception {
        Path directory = temporary.resolve(String.format("shard-%05d", shard));
        try (var writer = new MoleculeFingerprintIndexWriter(directory)) {
            for (int i = 0; i < rows.length; i++) {
                Row row = rows[i];
                writer.write(shard * 10L + i, 6, row.id, row.smiles,
                        new float[128], row.compact);
            }
        }
        Files.writeString(directory.resolve(".complete"), "complete\n");
    }

    private MoleculeFingerprintIndexManifest manifest() {
        var value = new MoleculeFingerprintIndexManifest();
        value.input = temporary.resolve("input.tsv").toString();
        value.smilesColumn = "smiles";
        value.idColumn = "id";
        value.modelBundle = temporary.resolve("model").toString();
        value.compactBundle = temporary.resolve("compact").toString();
        value.sourceRowsPerShard = 2;
        value.sourceRowCount = 3;
        value.recordCount = 3;
        value.rejectedCount = 0;
        var first = shard(0, 0, 2);
        var second = shard(1, 2, 1);
        value.shards = List.of(first, second);
        value.runtime = Map.of();
        value.buildStatistics = Map.of();
        value.validate();
        return value;
    }

    private static MoleculeFingerprintIndexShard shard(int index, long first, long count) {
        var value = new MoleculeFingerprintIndexShard();
        value.shardIndex = index;
        value.directory = String.format("shard-%05d", index);
        value.firstSourceRow = first;
        value.sourceRowCount = count;
        value.recordCount = count;
        value.rejectedCount = 0;
        return value;
    }

    private static MoleculeSkelSpheresHit hit(String id, String smiles, double predicted) {
        return new MoleculeSkelSpheresHit(
                new com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference(0,
                        "benzene".equals(id) ? 0 : 1),
                new com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata(
                        0, 6, id, smiles), predicted, predicted, null);
    }

    private static float[] vector(float first, float second) {
        float[] result = new float[16];
        result[0] = first;
        result[1] = second;
        return result;
    }

    private static CompactSkelSpheresManifest compactManifest() {
        CompactSkelSpheresManifest value = new CompactSkelSpheresManifest();
        String hash = "0".repeat(64);
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
        value.sourceFoundationCheckpointSha256 = hash;
        value.sourcePredictorCheckpointSha256 = hash;
        value.projectionCheckpointSha256 = hash;
        value.projectionSha256 = hash;
        return value;
    }

    private record Row(String id, String smiles, float[] compact) {}
}
