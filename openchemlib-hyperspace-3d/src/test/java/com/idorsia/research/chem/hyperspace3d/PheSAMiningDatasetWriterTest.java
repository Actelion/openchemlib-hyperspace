package com.idorsia.research.chem.hyperspace3d;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAMiningDatasetWriter;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAMiningManifest;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAQueryPairRecord;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.junit.jupiter.api.Assertions.*;

class PheSAMiningDatasetWriterTest {
    @TempDir Path temporary;
    @Test void finalizesACompleteHashBoundCache() throws Exception {
        float[] query = new float[128]; query[0] = 1;
        float[] candidate = new float[128]; candidate[1] = 1;
        var molecules = List.of(
                new PheSAMiningDatasetWriter.MoleculeEntry(0, "query", "A", "query", 0,
                        "CC", "query", -1, -1, 2, 0, 1, 0, query),
                new PheSAMiningDatasetWriter.MoleculeEntry(1, "candidate", "MINING_TRAIN",
                        "candidate", 1, "CCC", "candidate", 0, 5, 3, 0, 1, 0, candidate));
        var pair = new PheSAQueryPairRecord(1, .6f, .7f, .5f, 1, 2, 3,
                .6f, .7f, .5f, 1, 0, 0);
        var provenance = new PheSAMiningDatasetWriter.Provenance("a".repeat(64),
                "b".repeat(64), "d".repeat(64), "c".repeat(64), "OCL-test", .5, 8, 0, 17);
        Path output = temporary.resolve("cache");
        new PheSAMiningDatasetWriter().write(output, provenance, molecules,
                List.of(new PheSAMiningDatasetWriter.QueryEntry(0, "query", "A", 100, 0,
                        List.of(pair))));
        assertTrue(Files.isRegularFile(output.resolve(".complete")));
        assertEquals(512, Files.size(output.resolve("descriptors.f16")));
        var manifest = new ObjectMapper().readValue(output.resolve("manifest.json").toFile(),
                PheSAMiningManifest.class);
        manifest.validate(); assertEquals(1, manifest.pairCount);
        assertThrows(java.io.IOException.class, () -> new PheSAMiningDatasetWriter().write(
                output, provenance, molecules, List.of()));
    }
}
