package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.build.CompactProductVectorIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.CompactProductVectorIndexBuilder;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class CompactProductVectorIndexBuilderTest {
    @TempDir Path temporary;

    @Test void buildsDirectlyIntoRequestedColumnarShards() throws Exception {
        RawSynthonSpace space = space();
        CompactProductVectorIndexBuildConfig config = new CompactProductVectorIndexBuildConfig();
        config.inputs.rawFull = "full"; config.inputs.rawDownsampled = "downsampled";
        config.inputs.modelBundle = "model"; config.inputs.compactBundle = "compact";
        config.output.directory = "output"; config.output.targetRecordCount = 4;
        config.output.recordsPerShard = 2; config.runtime.encoderBatchSize = 3;
        config.runtime.cpuWorkers = 2; config.filters.maxRotatableBonds = 0; config.validate();
        String hash = "0".repeat(64);
        var provenance = new CompactProductVectorIndexBuilder.Provenance(
                "space:1.0", hash, "space:1.0", hash, hash, hash, "test-model",
                CompactProductVectorIndexBuilder.configurationHash(config));
        var encoder = (com.idorsia.research.chem.hyperspace3d.model.ProductEmbeddingBatchEncoder) batch -> {
            float[][] output = new float[batch.batchSize()][16];
            for (int i = 0; i < output.length; i++) output[i][i % 16] = 1;
            return output;
        };
        Path output = temporary.resolve("index");
        var built = new CompactProductVectorIndexBuilder(space, space, encoder, config, provenance).build(output);
        assertEquals(4, built.recordCount());
        assertEquals(2, built.shardCount());
        var manifest = ProductVectorIndexReader.loadManifest(output, true);
        assertEquals(List.of(2L, 2L), manifest.shards.stream().map(shard -> shard.recordCount()).toList());
        assertTrue(manifest.buildStatistics.containsKey("accepted"));
    }

    private static RawSynthonSpace space() throws Exception {
        return RawSynthonSpace.builder("space")
                .addRawFragments("r", 0, List.of(synthon(0, "a", "[U]CCC"), synthon(0, "b", "[U]CCCC")))
                .addRawFragments("r", 1, List.of(synthon(1, "x", "[U]CCC"), synthon(1, "y", "[U]CCO")))
                .build();
    }
    private static RawSynthon synthon(int position, String id, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule(); new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule("r", position, id, molecule);
    }
}
