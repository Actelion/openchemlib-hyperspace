package com.idorsia.research.chem.hyperspace.downsampling;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class SkelSpheresKCentersDownsamplerTest {

    @Test
    void downsamplingToyReactionIsDeterministic() throws Exception {
        SynthonSpace space = loadToySpace("testdata/idorsia_toy_space_a.txt",
                Collections.singleton("benzoimidazole_b-8"));

        SynthonDownsamplingRequest request = SynthonDownsamplingRequest.builder()
                .withMaxCenters(8)
                .withMinSimilarity(0.75)
                .withRandomSeed(42L)
                .enforceConnectorEquivalence(true)
                .build();

        SynthonDownsampler downsampler = new SkelSpheresKCentersDownsampler();
        SynthonDownsamplingResult result = new SynthonDownsamplingOrchestrator()
                .downsample(space, downsampler, request);

        assertEquals(3, result.getSetResults().size());
        long totalOriginal = result.getTotalOriginalSynthons();
        long totalRetained = result.getTotalRetainedSynthons();
        assertEquals(50, totalOriginal);
        assertEquals(19, totalRetained);
        Map<Integer, Integer> retainedPerFragment = new HashMap<>();
        result.getSetResults().forEach(r -> {
            assertTrue(r.getRetainedSize() <= request.getMaxCenters());
            retainedPerFragment.put(r.getFragType().frag, r.getRetainedSize());
        });
        assertEquals(8, retainedPerFragment.get(0));
        assertEquals(8, retainedPerFragment.get(1));
        assertEquals(3, retainedPerFragment.get(2));
    }

    private SynthonSpace loadToySpace(String resource, Set<String> allowedReactions) throws Exception {
        InputStream in = getClass().getClassLoader().getResourceAsStream(resource);
        assertNotNull(in, "Test data file missing: " + resource);
        Map<String, Map<Integer, List<Object>>> moleculesByReaction = new HashMap<>();
        Map<String, Map<String, String>> idcodeToIds = new HashMap<>();
        SmilesParser parser = new SmilesParser();

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
            reader.readLine(); // header
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\t");
                if (parts.length < 4) {
                    continue;
                }
                String rxnId = parts[3];
                if (!allowedReactions.contains(rxnId)) {
                    continue;
                }
                StereoMolecule mol = new StereoMolecule();
                parser.parse(mol, parts[0]);
                String idcode = mol.getIDCode();
                int synthonIdx = Integer.parseInt(parts[2]);
                String fragmentId = parts[1];
                moleculesByReaction.computeIfAbsent(rxnId, key -> new HashMap<>())
                        .computeIfAbsent(synthonIdx, key -> new ArrayList<>())
                        .add(idcode);
                idcodeToIds.computeIfAbsent(rxnId, key -> new HashMap<>()).put(idcode, fragmentId);
            }
        } catch (IOException ex) {
            fail("Unable to read toy space", ex);
        }

        SynthonSpace space = new SynthonSpace();
        space.setFP(new DescriptorHandlerLongFFP1024_plus("ffp"), 1024);
        for (String rxnId : moleculesByReaction.keySet()) {
            space.addReaction(rxnId, moleculesByReaction.get(rxnId), idcodeToIds.get(rxnId), null);
        }
        space.initAfterJavaDeserialization();
        space.reinitBitTree();
        return space;
    }
}
