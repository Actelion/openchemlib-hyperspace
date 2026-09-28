package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexShard;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexWriter;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresManifest;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.screening.CompactSkelSpheresScreener;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class ProductVectorIndexTest {
    @TempDir Path temporary;

    @Test void columnarRoundTripAndStreamingElites() throws Exception {
        Path vectors = temporary.resolve("shard-00000.vec");
        Path rows = temporary.resolve("shard-00000.rows");
        Path tuples = temporary.resolve("shard-00000.tuples");
        try (var writer = new ProductVectorIndexWriter(vectors, rows, tuples, List.of("r1", "r2"))) {
            writer.write(tuple("r1", "a"), vector(1f, 0f));
            writer.write(tuple("r1", "b"), vector(0.8f, 0.2f));
            writer.write(tuple("r2", "c"), vector(-1f, 0f));
            assertEquals(3, writer.recordCount());
        }
        ProductVectorIndexManifest manifest = manifest(vectors, rows, tuples);
        new com.fasterxml.jackson.databind.ObjectMapper().writeValue(
                temporary.resolve("manifest.json").toFile(), manifest);
        var loaded = ProductVectorIndexReader.loadManifest(temporary, true);
        try (var reader = new ProductVectorIndexReader(temporary, loaded, 0)) {
            var first = reader.readBatch(2);
            assertEquals(2, first.vectors().length);
            assertEquals("b", reader.resolveTuple(first.references().get(1)).synthonIds().get(0));
            assertEquals(0.8f, first.vectors()[1][0], 0.001f);
        }
        var result = new CompactSkelSpheresScreener(new CompactSkelSpheresScorer(compactManifest()))
                .screen(temporary, loaded, vector(1f, 0f), 2, 1, 2);
        assertEquals(3, result.recordsScanned());
        assertEquals(List.of("a", "b"), result.global().stream()
                .map(hit -> hit.tuple().synthonIds().get(0)).toList());
        assertEquals("a", result.byReaction().get("r1").get(0).tuple().synthonIds().get(0));
        assertEquals("c", result.byReaction().get("r2").get(0).tuple().synthonIds().get(0));
    }

    @Test void checksumsAndManifestContractAreEnforced() throws Exception {
        Path vectors = temporary.resolve("shard-00000.vec");
        Path rows = temporary.resolve("shard-00000.rows");
        Path tuples = temporary.resolve("shard-00000.tuples");
        try (var writer = new ProductVectorIndexWriter(vectors, rows, tuples, List.of("r1"))) {
            writer.write(tuple("r1", "a"), vector(1f, 0f));
        }
        ProductVectorIndexManifest manifest = manifest(vectors, rows, tuples);
        manifest.embeddingDimension = 128;
        assertThrows(IllegalArgumentException.class, manifest::validate);
        manifest.embeddingDimension = 16;
        new com.fasterxml.jackson.databind.ObjectMapper().writeValue(
                temporary.resolve("manifest.json").toFile(), manifest);
        Files.write(vectors, new byte[]{1}, java.nio.file.StandardOpenOption.APPEND);
        assertThrows(java.io.IOException.class,
                () -> ProductVectorIndexReader.loadManifest(temporary, true));
    }

    private ProductVectorIndexManifest manifest(Path vectors, Path rows, Path tuples) throws Exception {
        ProductVectorIndexManifest value = new ProductVectorIndexManifest();
        String hash = "0".repeat(64);
        value.rawspaceIdentity = "full:1"; value.rawspaceSha256 = hash;
        value.downsampledRawspaceIdentity = "small:1"; value.downsampledRawspaceSha256 = hash;
        value.sourceModelBundleHash = hash; value.projectionModelBundleHash = hash;
        value.modelVersion = "test"; value.samplingAlgorithm = "test"; value.samplingSeed = 17;
        value.reactionWeightConfiguration = Map.of("mode", "test");
        value.molecularFilters = Map.of("maxAtoms", 32); value.recordCount = 3;
        value.reactionDictionary = List.of("r1", "r2");
        value.shards = List.of(new ProductVectorIndexShard(vectors.getFileName().toString(),
                rows.getFileName().toString(), tuples.getFileName().toString(), 3,
                ProductFingerprintIndexBuilder.sha256(vectors),
                ProductFingerprintIndexBuilder.sha256(rows),
                ProductFingerprintIndexBuilder.sha256(tuples)));
        value.validate(); return value;
    }
    private static ProductTuple tuple(String reaction, String id) {
        return new ProductTuple(reaction, List.of(id), List.of(0));
    }
    private static float[] vector(float first, float second) {
        float[] value = new float[16]; value[0] = first; value[1] = second; return value;
    }
    private static CompactSkelSpheresManifest compactManifest() {
        CompactSkelSpheresManifest value = new CompactSkelSpheresManifest();
        String hash = "0".repeat(64);
        value.artifactType = "deepspace7-skelspheres16-onnx-bundle"; value.formatVersion = 1;
        value.modelVersion = "test"; value.architecture = "mlp_128_128_16_l2";
        value.target = "skelspheres_similarity"; value.inputEmbeddingDim = 128;
        value.hiddenDim = 128; value.outputEmbeddingDim = 16; value.l2Normalized = true;
        value.canonicalSeed = 17; value.bestEpoch = 1; value.recommendedStorageDtype = "float16";
        value.calibration = new CompactSkelSpheresManifest.Calibration();
        value.calibration.type = "sigmoid_affine_cosine"; value.calibration.scale = 1; value.calibration.bias = 0;
        value.sourceFoundationCheckpointSha256 = hash; value.sourcePredictorCheckpointSha256 = hash;
        value.projectionCheckpointSha256 = hash; value.projectionSha256 = hash; return value;
    }
}
