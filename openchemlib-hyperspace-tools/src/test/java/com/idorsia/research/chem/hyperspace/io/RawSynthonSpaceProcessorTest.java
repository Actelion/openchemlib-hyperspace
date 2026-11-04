package com.idorsia.research.chem.hyperspace.io;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.junit.jupiter.api.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

class RawSynthonSpaceProcessorTest {

    @Test
    void descriptorAugmentationAddsAttributes() throws Exception {
        Path sample = Paths.get("..", "openchemlib-hyperspace-core", "src", "main", "resources",
                "testdata", "idorsia_toy_space_a.txt").toAbsolutePath().normalize();
        RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
        RawSynthonSpace raw = importer.importEnamine(SynthonSpaceParser3.EnamineOptions.builder()
                .input(sample)
                .spaceName("toy")
                .mode("FragFp")
                .threads(1)
                .maxSynthonSets(3)
                .buildSynthonSpace(false)
                .build()).getRawSpace();

        RawSynthonSpace processed = RawSynthonSpaceProcessor.process(raw,
                List.of(new DescriptorAugmentationTask(List.of("FragFp"), 1, true)));

        String tags = processed.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS);
        assertTrue(tags != null && tags.contains("FragFp"), "descriptor tags should include FragFp");
        RawSynthonSpace.ReactionData data = processed.getReactions().get("benzoimidazole_b-8");
        assertNotNull(data, "expected reaction benzoimidazole_b-8");
        Map<String, Map<String, String>> attributes = data.getFragmentAttributes();
        assertTrue(attributes.values().stream().anyMatch(map ->
                map != null && map.containsKey("descriptor.FragFp")), "Fragment missing descriptor attribute");
    }
}
