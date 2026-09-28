package com.idorsia.research.chem.hyperspace3d;

import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.screening.*;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class ProductIndexScreenerTest {
    @TempDir Path temporary;

    @Test void roundTripsAndCollectsGlobalAndReactionElites() throws Exception {
        Path shard = temporary.resolve("products.h3di");
        ProductTuple a = tuple("r1", "a");
        ProductTuple b = tuple("r1", "b");
        ProductTuple c = tuple("r2", "c");
        try (ProductFingerprintIndexWriter writer = new ProductFingerprintIndexWriter(shard)) {
            writer.write(record(a, 0.2f));
            writer.write(record(b, 0.9f));
            writer.write(record(c, 0.7f));
            writer.write(record(b, 0.9f));
            assertEquals(4, writer.recordCount());
        }
        var comparator = (com.idorsia.research.chem.hyperspace3d.model.EmbeddingComparator)
                (query, candidates) -> {
                    float[][] scores = new float[candidates.length][6];
                    for (int i = 0; i < candidates.length; i++) scores[i][3] = candidates[i][0];
                    return scores;
                };
        try (ProductFingerprintIndexReader reader = new ProductFingerprintIndexReader(shard)) {
            FingerprintIndexScreener.Result result = new FingerprintIndexScreener(
                    comparator, DeepSpaceModelManifestTest.valid()).screen(reader,
                    new float[128], ScreeningObjective.direct("phesa_total"), 2, 2, 1);
            assertEquals(3, result.uniqueTupleCount());
            assertEquals(b, result.global().get(0).tuple());
            assertEquals(c, result.byReaction().get("r2").get(0).tuple());
            assertEquals(1, result.byReaction().get("r1").size());
        }
    }

    private static ProductTuple tuple(String reaction, String id) {
        return new ProductTuple(reaction, List.of(id), List.of(0));
    }
    private static ProductFingerprintRecord record(ProductTuple tuple, float score) {
        float[] embedding = new float[128]; embedding[0] = score;
        return new ProductFingerprintRecord(tuple, embedding);
    }
}
