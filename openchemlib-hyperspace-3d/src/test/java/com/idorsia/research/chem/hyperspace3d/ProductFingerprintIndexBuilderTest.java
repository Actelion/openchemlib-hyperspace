package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.screening.ReactionScheduler;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintRecord;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductTupleSampler;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class ProductFingerprintIndexBuilderTest {
    @TempDir Path temporary;

    @Test void samplerIsDeterministicUniqueAndUsesFullOrdinals() {
        RawSynthonSpace full = arbitrarySpace("space", List.of("a", "b", "c"), List.of("x", "y"));
        RawSynthonSpace reduced = arbitrarySpace("space", List.of("c", "a"), List.of("y", "x"));
        var weighting = ReactionScheduler.Weighting.exponent(0.01, 1.0);
        ProductTupleSampler first = new ProductTupleSampler(full, reduced, 17, weighting);
        ProductTupleSampler second = new ProductTupleSampler(full, reduced, 17, weighting);
        Set<ProductTuple> observed = new HashSet<>();
        for (int i = 0; i < 4; i++) {
            ProductTuple one = first.next("r").tuple();
            ProductTuple two = second.next("r").tuple();
            assertEquals(one, two);
            assertTrue(observed.add(one));
            assertEquals(one.synthonIds().get(0).equals("c") ? 2 : 0,
                    one.synthonOrdinals().get(0));
            assertEquals(one.synthonIds().get(1).equals("y") ? 1 : 0,
                    one.synthonOrdinals().get(1));
        }
        assertTrue(first.exhausted("r"));
        assertNull(first.next("r"));
    }

    @Test void rejectsDownsampledStructureMismatch() {
        RawSynthonSpace full = arbitrarySpace("space", List.of("a"), List.of("x"));
        RawSynthonSpace reduced = RawSynthonSpace.builder("space")
                .addRawFragments("r", 0, List.of(raw("r", 0, "a", "different")))
                .addRawFragments("r", 1, List.of(raw("r", 1, "x", "x"))).build();
        assertThrows(IllegalArgumentException.class, () -> new ProductTupleSampler(full, reduced,
                1, ReactionScheduler.Weighting.exponent(0.01, 1.0)));
    }

    @Test void buildsShardsAndCompletedRunIsIdempotent() throws Exception {
        RawSynthonSpace full = chemicalSpace();
        ProductFingerprintIndexBuildConfig config = config(4, 2);
        String configHash = ProductFingerprintIndexBuilder.configurationHash(config);
        var provenance = new ProductFingerprintIndexBuilder.Provenance(
                "chemical:1.0", "raw-hash", "chemical:1.0", "down-hash", "bundle-hash", configHash);
        var fakeEncoder = (com.idorsia.research.chem.hyperspace3d.model.ProductEmbeddingBatchEncoder) batch -> {
            float[][] output = new float[batch.batchSize()][128];
            for (int i = 0; i < output.length; i++) {
                int atoms = 0;
                for (int a = 0; a < 32; a++) if (batch.atomMask()[i * 32 + a]) atoms++;
                output[i][0] = atoms;
            }
            return output;
        };
        Path output = temporary.resolve("index");
        var result = new ProductFingerprintIndexBuilder(full, full, fakeEncoder,
                config, provenance).build(output);
        assertEquals(4, result.recordCount());
        assertEquals(2, result.shardCount());
        ProductFingerprintIndexManifest manifest = new ObjectMapper().readValue(
                result.manifestPath().toFile(), ProductFingerprintIndexManifest.class);
        manifest.validate();
        assertEquals(4, manifest.recordCount);
        assertEquals(ProductTupleSampler.ALGORITHM, manifest.samplingAlgorithm);
        Set<ProductTuple> tuples = new HashSet<>();
        for (var shard : manifest.shards) {
            try (ProductFingerprintIndexReader reader = new ProductFingerprintIndexReader(
                    output.resolve(shard.path()))) {
                List<ProductFingerprintRecord> records = reader.readBatch(10);
                assertEquals(2, records.size());
                records.forEach(record -> {
                    assertTrue(tuples.add(record.tuple()));
                    assertTrue(record.embedding()[0] >= 6);
                });
            }
        }
        var repeated = new ProductFingerprintIndexBuilder(full, full, fakeEncoder,
                config, provenance).build(output);
        assertEquals(4, repeated.recordCount());
        assertTrue(repeated.resumed());
    }

    @Test void configurationAppliesDeviceBatchDefaultsAndValidatesFloor() {
        ProductFingerprintIndexBuildConfig config = config(1, 1);
        config.runtime.encoderBatchSize = 0;
        config.runtime.device = "CUDA";
        config.validate();
        assertEquals(512, config.runtime.encoderBatchSize);
        config.filters.maxHeavyAtoms = 33;
        assertThrows(IllegalArgumentException.class, config::validate);
    }

    private static ProductFingerprintIndexBuildConfig config(long target, long shard) {
        ProductFingerprintIndexBuildConfig config = new ProductFingerprintIndexBuildConfig();
        config.inputs.rawFull = "full";
        config.inputs.rawDownsampled = "reduced";
        config.inputs.modelBundle = "bundle";
        config.output.directory = "output";
        config.output.targetRecordCount = target;
        config.output.recordsPerShard = shard;
        config.runtime.cpuWorkers = 2;
        config.runtime.encoderBatchSize = 2;
        config.runtime.progressIntervalSeconds = 0;
        config.filters.maxRotatableBonds = 0;
        config.validate();
        return config;
    }

    private static RawSynthonSpace arbitrarySpace(String name, List<String> zero, List<String> one) {
        List<RawSynthon> a = new ArrayList<>();
        List<RawSynthon> b = new ArrayList<>();
        zero.forEach(id -> a.add(raw("r", 0, id, id)));
        one.forEach(id -> b.add(raw("r", 1, id, id)));
        return RawSynthonSpace.builder(name).addRawFragments("r", 0, a)
                .addRawFragments("r", 1, b).build();
    }

    private static RawSynthonSpace chemicalSpace() throws Exception {
        List<RawSynthon> zero = List.of(chemical("r", 0, "a", "[U]CCC"),
                chemical("r", 0, "b", "[U]CCCC"));
        List<RawSynthon> one = List.of(chemical("r", 1, "x", "[U]CCC"),
                chemical("r", 1, "y", "[U]CCO"));
        return RawSynthonSpace.builder("chemical").addRawFragments("r", 0, zero)
                .addRawFragments("r", 1, one).build();
    }

    private static RawSynthon chemical(String reaction, int position, String id, String smiles)
            throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule(reaction, position, id, molecule);
    }

    private static RawSynthon raw(String reaction, int position, String id, String idcode) {
        BitSet connectors = new BitSet();
        connectors.set(0);
        return new RawSynthon(reaction, position, id, idcode, connectors);
    }
}
