package com.idorsia.research.chem.hyperspace3d;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.feature.FeaturizationResult;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class GoldenFeatureParityTest {
    @Test void strictV2ChannelsMatchPythonFixture() throws Exception {
        InputStream resource = getClass().getResourceAsStream(
                "/com/idorsia/research/chem/hyperspace3d/deepspace7-v1-v3graph-golden.json.gz");
        assertNotNull(resource);
        JsonNode fixture = new ObjectMapper().readTree(new GZIPInputStream(resource));
        OCLDeepSpaceFeaturizer featurizer = new OCLDeepSpaceFeaturizer();
        for (JsonNode record : fixture.get("records")) {
            if (record.has("rejection")) continue;
            FeaturizationResult actual = featurizer.featurize(
                    OCLDeepSpaceFeaturizerTest.molecule(record.get("inputSmiles").asText()));
            assertTrue(actual.accepted(), record.get("id").asText());
            JsonNode expectedAtom = record.get("atomX");
            for (int atom = 0; atom < 32; atom++) {
                for (int feature = 0; feature < 50; feature++) {
                    int index = atom * 56 + feature;
                    assertEquals(expectedAtom.get(index).floatValue(), actual.atomFeatures()[index],
                            0f, record.get("id").asText() + " atom feature " + index);
                }
            }
            JsonNode expectedPair = record.get("pairX");
            for (int i = 0; i < actual.pairFeatures().length; i++) {
                assertEquals(expectedPair.get(i).floatValue(), actual.pairFeatures()[i],
                        0f, record.get("id").asText() + " pair feature " + i);
            }
        }
    }
}
