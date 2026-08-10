package com.idorsia.research.chem.hyperspace.rawspace;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.downsampling.SkelSpheresKCentersDownsampler;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonSetDownsamplingResult;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsampler;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingOrchestrator;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingRequest;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingResult;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class RawSynthonSpaceTest {

    @Test
    void fullRawSpaceJsonRoundTripDoesNotEmitEmbeddedDownsampledSets() throws Exception {
        SynthonSpace space = loadToySpace("testdata/idorsia_toy_space_a.txt",
                Collections.singleton("benzoimidazole_b-8"));

        RawSynthonSpace.Builder builder = RawSynthonSpace.builder("toy")
                .withFullSynthonSpace(space);
        String sampleFragId = space.getSynthonSet("benzoimidazole_b-8", 0).get(0).fragment_id;
        builder.addFragmentAttribute("benzoimidazole_b-8", sampleFragId, "price", "42");
        builder.addReactionMetadata("benzoimidazole_b-8", "source.space", "toy_vendor");
        RawSynthonSpace raw = builder.build();

        Map<Integer, List<SynthonSpace.FragId>> fragments = raw.getFragmentSets("benzoimidazole_b-8");
        assertEquals(3, fragments.size());

        Path temp = Files.createTempFile("rawspace", ".rawspace");
        Path tempGz = Files.createTempFile("rawspace", ".rawspace.gz");
        try {
            RawSynthonSpaceIO.write(raw, temp);
            String json = Files.readString(temp, StandardCharsets.UTF_8);
            assertFalse(json.contains("downsampledFragmentSets"));
            assertFalse(json.contains("downsamplingAlgorithm"));
            assertFalse(json.contains("downsamplingRequest"));

            RawSynthonSpace loaded = RawSynthonSpaceIO.read(temp);
            assertEquals(raw.getFragmentSets("benzoimidazole_b-8").get(0).size(),
                    loaded.getFragmentSets("benzoimidazole_b-8").get(0).size());
            assertNotNull(loaded.findFragment("benzoimidazole_b-8",
                    raw.getFragmentSets("benzoimidazole_b-8").get(0).get(0).fragment_id));
            assertEquals("42", loaded.getReactions().get("benzoimidazole_b-8")
                    .getFragmentAttributes().get(sampleFragId).get("price"));
            assertEquals("toy_vendor", loaded.getReactions().get("benzoimidazole_b-8")
                    .getReactionMetadata().get("source.space"));

            RawSynthonSpaceIO.write(raw, tempGz);
            RawSynthonSpace loadedGz = RawSynthonSpaceIO.read(tempGz);
            assertEquals(loaded.getFragmentSets("benzoimidazole_b-8").get(1).size(),
                    loadedGz.getFragmentSets("benzoimidazole_b-8").get(1).size());
        } finally {
            Files.deleteIfExists(temp);
            Files.deleteIfExists(tempGz);
        }
    }

    @Test
    void downsampledRawSpaceStoresRepresentativesInNormalFragmentSets() throws Exception {
        SynthonSpace space = loadToySpace("testdata/idorsia_toy_space_a.txt",
                Collections.singleton("benzoimidazole_b-8"));

        SynthonDownsamplingRequest request = SynthonDownsamplingRequest.builder()
                .withMaxCenters(8)
                .withMinSimilarity(0.75)
                .withRandomSeed(42L)
                .enforceConnectorEquivalence(true)
                .withSizeCapScale(1.5)
                .withSizeCapOffset(2.0)
                .build();

        SynthonDownsampler downsampler = new SkelSpheresKCentersDownsampler();
        SynthonDownsamplingResult result = new SynthonDownsamplingOrchestrator()
                .downsample(space, downsampler, request);

        RawSynthonSpace.Builder builder = RawSynthonSpace.builder("toy_downsampled")
                .withDownsamplingMetadata(downsampler.getName(), request);
        for (SynthonSetDownsamplingResult setResult : result.getSetResults()) {
            builder.addFragments(setResult.getFragType().rxn_id,
                    setResult.getFragType().frag,
                    setResult.getRepresentatives());
        }
        RawSynthonSpace raw = builder.build();

        assertEquals(RawSynthonSpace.MetadataKeys.SPACE_ROLE_DOWNSAMPLED,
                raw.getMetadata().get(RawSynthonSpace.MetadataKeys.SPACE_ROLE));
        assertEquals(downsampler.getName(),
                raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_ALGORITHM));
        assertEquals("8", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_MAX_CENTERS));
        assertEquals("0.75", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_MIN_SIMILARITY));
        assertEquals("42", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SEED));
        assertEquals("true", raw.getMetadata().get(
                RawSynthonSpace.MetadataKeys.DOWNSAMPLING_ENFORCE_CONNECTOR_EQUIVALENCE));
        assertEquals("1.5", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SIZE_CAP_SCALE));
        assertEquals("2.0", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SIZE_CAP_OFFSET));

        Path temp = Files.createTempFile("rawspace-downsampled", ".rawspace");
        try {
            RawSynthonSpaceIO.write(raw, temp);
            String json = Files.readString(temp, StandardCharsets.UTF_8);
            assertFalse(json.contains("downsampledFragmentSets"));

            RawSynthonSpace loaded = RawSynthonSpaceIO.read(temp);
            for (SynthonSetDownsamplingResult setResult : result.getSetResults()) {
                List<SynthonSpace.FragId> loadedSet = loaded.getFragmentSets(setResult.getFragType().rxn_id)
                        .get(setResult.getFragType().frag);
                assertEquals(setResult.getRepresentatives().size(), loadedSet.size());
            }
        } finally {
            Files.deleteIfExists(temp);
        }
    }

    @Test
    void rawSpaceCanRehydrateSynthonSpace() throws Exception {
        SynthonSpace original = loadToySpace("testdata/idorsia_toy_space_a.txt",
                Collections.singleton("benzoimidazole_b-8"));
        RawSynthonSpace raw = RawSynthonSpace.builder("toy")
                .withFullSynthonSpace(original)
                .build();

        SynthonSpace rebuilt = RawSynthonSpaceAssembler.buildSynthonSpace(raw,
                RawSynthonSpaceAssembler.BuildOptions.builder()
                        .descriptorShortName("FragFp")
                        .build());
        assertEquals(original.getRxnIds().size(), rebuilt.getRxnIds().size());

        for (String rxn : original.getRxnIds()) {
            Map<Integer, SynthonSpace.FragType> types = original.getFragTypes(rxn);
            assertNotNull(types);
            for (Integer idx : types.keySet()) {
                List<SynthonSpace.FragId> originalSet = original.getSynthonSet(rxn, idx);
                List<SynthonSpace.FragId> rebuiltSet = rebuilt.getSynthonSet(rxn, idx);
                assertNotNull(rebuiltSet, "Missing synthon set for " + rxn + ":" + idx);
                assertEquals(originalSet.size(), rebuiltSet.size(),
                        "Mismatched synthon count for " + rxn + ":" + idx);
            }
        }
    }

    private SynthonSpace loadToySpace(String resource, Set<String> allowedReactions) throws Exception {
        InputStream in = getClass().getClassLoader().getResourceAsStream(resource);
        assertNotNull(in, "Test data file missing: " + resource);
        Map<String, Map<Integer, List<Object>>> moleculesByReaction = new HashMap<>();
        Map<String, Map<String, String>> idcodeToIds = new HashMap<>();
        SmilesParser parser = new SmilesParser();

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
            reader.readLine();
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
