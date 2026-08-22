package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildException;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.ProductEmbeddingBatchEncoder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class ProductFingerprintIndexResumeTest {
    @TempDir Path temporary;

    @Test void resumesFromLastCommittedShardWithoutDuplicateTuples() throws Exception {
        RawSynthonSpace space = chemicalSpace();
        ProductFingerprintIndexBuildConfig config = config();
        var provenance = new ProductFingerprintIndexBuilder.Provenance(
                "chemical:1.0", "raw", "chemical:1.0", "reduced", "bundle",
                ProductFingerprintIndexBuilder.configurationHash(config));
        AtomicInteger batches = new AtomicInteger();
        ProductEmbeddingBatchEncoder failing = batch -> {
            if (batches.incrementAndGet() == 3) throw new IllegalStateException("simulated interruption");
            return embeddings(batch.batchSize());
        };
        Path output = temporary.resolve("resume-index");
        assertThrows(ProductFingerprintIndexBuildException.class, () ->
                new ProductFingerprintIndexBuilder(space, space, failing, config, provenance).build(output));
        assertTrue(Files.isRegularFile(output.resolve("shard-00000.h3di")));
        assertFalse(Files.exists(output.resolve("manifest.json")));

        ProductEmbeddingBatchEncoder working = batch -> embeddings(batch.batchSize());
        var result = new ProductFingerprintIndexBuilder(space, space, working,
                config, provenance).build(output);
        assertTrue(result.resumed());
        assertEquals(4, result.recordCount());
        Set<ProductTuple> tuples = new HashSet<>();
        for (int shard = 0; shard < 2; shard++) {
            try (ProductFingerprintIndexReader reader = new ProductFingerprintIndexReader(
                    output.resolve(String.format("shard-%05d.h3di", shard)))) {
                reader.readBatch(10).forEach(record -> assertTrue(tuples.add(record.tuple())));
            }
        }
        assertEquals(4, tuples.size());
    }

    private static float[][] embeddings(int batch) {
        float[][] result = new float[batch][128];
        for (float[] embedding : result) embedding[0] = 1f;
        return result;
    }

    private static ProductFingerprintIndexBuildConfig config() {
        ProductFingerprintIndexBuildConfig config = new ProductFingerprintIndexBuildConfig();
        config.inputs.rawFull = "full";
        config.inputs.rawDownsampled = "reduced";
        config.inputs.modelBundle = "bundle";
        config.output.directory = "output";
        config.output.targetRecordCount = 4;
        config.output.recordsPerShard = 2;
        config.runtime.cpuWorkers = 2;
        config.runtime.encoderBatchSize = 2;
        config.runtime.progressIntervalSeconds = 0;
        config.filters.maxRotatableBonds = 0;
        config.validate();
        return config;
    }

    private static RawSynthonSpace chemicalSpace() throws Exception {
        return RawSynthonSpace.builder("chemical")
                .addRawFragments("r", 0, List.of(chemical(0, "a", "[U]CCC"),
                        chemical(0, "b", "[U]CCCC")))
                .addRawFragments("r", 1, List.of(chemical(1, "x", "[U]CCC"),
                        chemical(1, "y", "[U]CCO"))).build();
    }

    private static RawSynthon chemical(int position, String id, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule("r", position, id, molecule);
    }
}
