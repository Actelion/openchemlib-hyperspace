package com.idorsia.research.chem.hyperspace3d;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class JavaOnnxParityTest {
    @Test void javaOrtMatchesPythonOrtGoldenOutputs() throws Exception {
        String bundlePath = System.getProperty("hyperspace3d.modelBundle");
        Assumptions.assumeTrue(bundlePath != null,
                "set -Dhyperspace3d.modelBundle to run ONNX integration parity");
        JsonNode fixture;
        try (InputStream resource = getClass().getResourceAsStream(
                "/com/idorsia/research/chem/hyperspace3d/deepspace7-v1-v3graph-golden.json.gz")) {
            fixture = new ObjectMapper().readTree(new GZIPInputStream(resource));
        }
        List<JsonNode> records = new ArrayList<>();
        fixture.get("records").forEach(record -> { if (!record.has("rejection")) records.add(record); });
        int batch = records.size();
        float[] atom = new float[batch * 32 * 56];
        float[] pair = new float[batch * 32 * 32 * 36];
        boolean[] mask = new boolean[batch * 32];
        for (int b = 0; b < batch; b++) {
            copyFloats(records.get(b).get("atomX"), atom, b * 32 * 56);
            copyFloats(records.get(b).get("pairX"), pair, b * 32 * 32 * 36);
            for (int i = 0; i < 32; i++) mask[b * 32 + i] = records.get(b).get("atomMask").get(i).asBoolean();
        }
        DeepSpaceModelBundle bundle = DeepSpaceModelBundle.load(Path.of(bundlePath));
        DeepSpaceOnnxEnvironment runtime =
                new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.CPU);
        try (DeepSpaceV1Encoder encoder = new DeepSpaceV1Encoder(runtime, bundle);
             DeepSpaceV1Comparator comparator = new DeepSpaceV1Comparator(runtime, bundle)) {
            float[][] embeddings = encoder.encode(new DeepSpaceTensorBatch(batch, atom, pair, mask));
            for (int b = 0; b < batch; b++) {
                assertVector(records.get(b).get("embedding"), embeddings[b], 1e-5f);
            }
            float[][] scores = comparator.compare(embeddings[0], embeddings);
            for (int b = 0; b < batch; b++) {
                assertVector(records.get(b).get("scoresAgainstFirst"), scores[b], 1e-5f);
            }
            DeepSpaceTensorBatch one = new DeepSpaceTensorBatch(1,
                    java.util.Arrays.copyOfRange(atom, 0, 32 * 56),
                    java.util.Arrays.copyOfRange(pair, 0, 32 * 32 * 36),
                    java.util.Arrays.copyOfRange(mask, 0, 32));
            assertArrayEquals(embeddings[0], encoder.encode(one)[0], 1e-5f);
            DeepSpaceTensorBatch permuted = reverseAtomOrder(one);
            assertArrayEquals(embeddings[0], encoder.encode(permuted)[0], 1e-5f);
        }
    }

    private static void copyFloats(JsonNode source, float[] target, int offset) {
        for (int i = 0; i < source.size(); i++) target[offset + i] = source.get(i).floatValue();
    }

    private static DeepSpaceTensorBatch reverseAtomOrder(DeepSpaceTensorBatch source) {
        float[] atom = new float[32 * 56];
        float[] pair = new float[32 * 32 * 36];
        boolean[] mask = new boolean[32];
        for (int a = 0; a < 32; a++) {
            int oldA = 31 - a;
            mask[a] = source.atomMask()[oldA];
            System.arraycopy(source.atomFeatures(), oldA * 56, atom, a * 56, 56);
            for (int b = 0; b < 32; b++) {
                int oldB = 31 - b;
                System.arraycopy(source.pairFeatures(),
                        ((oldA * 32) + oldB) * 36,
                        pair, ((a * 32) + b) * 36, 36);
            }
        }
        return new DeepSpaceTensorBatch(1, atom, pair, mask);
    }

    private static void assertVector(JsonNode expected, float[] actual, float tolerance) {
        assertEquals(expected.size(), actual.length);
        for (int i = 0; i < actual.length; i++) {
            assertEquals(expected.get(i).floatValue(), actual[i], tolerance, "element " + i);
        }
    }
}
