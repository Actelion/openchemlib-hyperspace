package com.idorsia.research.chem.hyperspace3d;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresManifest;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class CompactSkelSpheresModelTest {
    @Test void manifestRejectsWrongDimensionsAndSource() {
        CompactSkelSpheresManifest manifest = validManifest();
        manifest.validate();
        manifest.outputEmbeddingDim = 32;
        assertThrows(IllegalArgumentException.class, manifest::validate);
    }

    @Test void javaProjectionAndCalibrationMatchPythonFixture() throws Exception {
        String primaryPath = System.getProperty("hyperspace3d.modelBundle");
        String compactPath = System.getProperty("hyperspace3d.compactBundle");
        Assumptions.assumeTrue(primaryPath != null && compactPath != null,
                "set primary and compact bundle properties to run compact ONNX parity");
        DeepSpaceModelBundle primary = DeepSpaceModelBundle.load(Path.of(primaryPath));
        CompactSkelSpheresModelBundle compact = CompactSkelSpheresModelBundle.load(
                Path.of(compactPath), primary.manifest());
        JsonNode fixture;
        try (InputStream resource = getClass().getResourceAsStream(
                "/com/idorsia/research/chem/hyperspace3d/deepspace7-skelspheres16-golden.json.gz")) {
            assertNotNull(resource);
            fixture = new ObjectMapper().readTree(new GZIPInputStream(resource));
        }
        List<JsonNode> records = new ArrayList<>();
        fixture.get("records").forEach(records::add);
        float[][] base = new float[records.size()][128];
        for (int row = 0; row < records.size(); row++) {
            copy(records.get(row).get("baseEmbedding"), base[row]);
        }
        DeepSpaceOnnxEnvironment runtime =
                new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.CPU);
        try (CompactSkelSpheresProjector projector =
                     new CompactSkelSpheresProjector(runtime, compact)) {
            float[][] actual = projector.project(base);
            CompactSkelSpheresScorer scorer =
                    new CompactSkelSpheresScorer(compact.manifest());
            for (int row = 0; row < actual.length; row++) {
                assertVector(records.get(row).get("compactEmbedding"), actual[row], 1e-6f);
                var score = scorer.score(actual[0], actual[row]);
                assertEquals(records.get(row).get("dotAgainstFirst").doubleValue(),
                        score.dotProduct(), 1e-6);
                assertEquals(records.get(row).get("calibratedScoreAgainstFirst").doubleValue(),
                        score.calibratedSimilarity(), 1e-6);
            }
        }
    }

    private static CompactSkelSpheresManifest validManifest() {
        CompactSkelSpheresManifest manifest = new CompactSkelSpheresManifest();
        manifest.artifactType = "deepspace7-skelspheres16-onnx-bundle";
        manifest.formatVersion = 1;
        manifest.modelVersion = "deepspace7-v2-skelspheres16-seed17";
        manifest.architecture = "mlp_128_128_16_l2";
        manifest.target = "skelspheres_similarity";
        manifest.inputEmbeddingDim = 128;
        manifest.hiddenDim = 128;
        manifest.outputEmbeddingDim = 16;
        manifest.l2Normalized = true;
        manifest.canonicalSeed = 17;
        manifest.bestEpoch = 45;
        manifest.recommendedStorageDtype = "float16";
        manifest.calibration = new CompactSkelSpheresManifest.Calibration();
        manifest.calibration.type = "sigmoid_affine_cosine";
        manifest.calibration.scale = 1.0;
        manifest.calibration.bias = 0.0;
        manifest.sourceFoundationCheckpointSha256 = "a".repeat(64);
        manifest.sourcePredictorCheckpointSha256 = "b".repeat(64);
        manifest.projectionCheckpointSha256 = "c".repeat(64);
        manifest.projectionSha256 = "d".repeat(64);
        return manifest;
    }

    private static void copy(JsonNode source, float[] target) {
        for (int i = 0; i < target.length; i++) target[i] = source.get(i).floatValue();
    }

    private static void assertVector(JsonNode expected, float[] actual, float tolerance) {
        assertEquals(expected.size(), actual.length);
        for (int i = 0; i < actual.length; i++) {
            assertEquals(expected.get(i).floatValue(), actual[i], tolerance, "element " + i);
        }
    }
}
