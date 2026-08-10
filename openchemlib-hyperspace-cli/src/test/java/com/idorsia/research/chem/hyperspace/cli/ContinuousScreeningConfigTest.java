package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.screening.ReactionScheduler;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

class ContinuousScreeningConfigTest {

    @TempDir
    Path tempDir;

    @Test
    void resolvesAndLoadsEncodedPhesaQueryFile() throws Exception {
        Path phesaPath = tempDir.resolve("query.phesa");
        String encoded = encodeToyPhesa();
        Files.writeString(phesaPath, "\n" + encoded + "\n");
        Path configPath = writeConfig("""
                  "query": { "phesaFile": "query.phesa" },
                """);

        ContinuousScreeningConfig config = ContinuousScreeningConfig.fromJson(configPath);
        ContinuousScreeningConfig.QueryInput query = config.resolveQuery(configPath);
        PheSAMolecule descriptor = ContinuousScreeningCLI.buildQueryDescriptor(query);

        assertEquals(phesaPath.toAbsolutePath().normalize().toString(), query.getPhesaFile());
        assertNotNull(descriptor);
        assertFalse(new DescriptorHandlerShape().calculationFailed(descriptor));
    }

    @Test
    void parsesBucketedReactionWeighting() throws Exception {
        Path configPath = writeConfig("""
                  "query": { "smiles": "CCO" },
                  "orchestration": {
                    "workerThreads": 2,
                    "queueCapacity": 10,
                    "reactionWeighting": {
                      "mode": "BUCKETED_PRODUCT",
                      "buckets": [
                        { "maxProductExclusive": 2000, "weight": 0.01 },
                        { "maxProductExclusive": 50000, "weight": 0.1 },
                        { "weight": 1.0 }
                      ]
                    }
                  },
                """);

        ContinuousScreeningConfig config = ContinuousScreeningConfig.fromJson(configPath);
        ReactionScheduler.Weighting weighting = config.toReactionWeighting();
        Map<String, List<Integer>> sizes = new LinkedHashMap<>();
        sizes.put("small", List.of(20, 99));
        sizes.put("medium", List.of(200, 100));
        sizes.put("large", List.of(500, 100));

        ReactionScheduler scheduler = new ReactionScheduler(sizes, weighting);

        assertEquals(0.01, scheduler.getWeightsByReaction().get("small"), 1e-12);
        assertEquals(0.1, scheduler.getWeightsByReaction().get("medium"), 1e-12);
        assertEquals(1.0, scheduler.getWeightsByReaction().get("large"), 1e-12);
    }


    @Test
    void parsesRuntimeLimitAndNullRandomSeed() throws Exception {
        Path configPath = writeConfig("""
                  "query": { "smiles": "CCO" },
                  "run": { "maxRuntime": "24h", "randomSeed": null },
                """);

        ContinuousScreeningConfig config = ContinuousScreeningConfig.fromJson(configPath);

        assertEquals(Duration.ofHours(24), config.getRun().resolveMaxRuntime());
        assertNull(config.getRun().getRandomSeed());
        assertEquals(config.getRun().getEffectiveRandomSeed(), config.getRun().getEffectiveRandomSeed());
    }

    private Path writeConfig(String queryAndOptionalOrchestrationJson) throws Exception {
        Path configPath = tempDir.resolve("screening.json");
        String orchestrationFallback = queryAndOptionalOrchestrationJson.contains("\"orchestration\"")
                ? ""
                : "  \"orchestration\": { \"workerThreads\": 2, \"queueCapacity\": 10 },\n";
        Files.writeString(configPath, """
                {
                  "inputs": { "rawFull": "full.rawspace.gz", "rawDownsampled": "down.rawspace.gz" },
                %s
                  "sampling": { "attemptsPerReaction": 1, "minAtoms": 1, "maxRotatableBonds": 10 },
                %s  "output": { "hitsTsv": "hits.tsv" }
                }
                """.formatted(queryAndOptionalOrchestrationJson, orchestrationFallback));
        return configPath;
    }

    private static String encodeToyPhesa() throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, "CCO");
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        ConformerSet conformers = new ConformerSetGenerator(1).generateConformerSet(molecule);
        DescriptorHandlerShape handler = new DescriptorHandlerShape();
        PheSAMolecule descriptor = handler.createDescriptor(conformers);
        if (handler.calculationFailed(descriptor)) {
            throw new IllegalStateException("Unable to create test PheSA descriptor");
        }
        return handler.encode(descriptor);
    }
}
