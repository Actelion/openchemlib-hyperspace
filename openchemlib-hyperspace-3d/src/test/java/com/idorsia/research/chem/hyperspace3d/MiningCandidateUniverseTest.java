package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintProvenance;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintShardReader;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorBatch;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.mining.MiningCandidateUniverse;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MiningCandidateUniverseTest {
    private static final String HASH_A = "a".repeat(64);
    private static final String HASH_B = "b".repeat(64);

    @Test
    void materializesVectorsAndChemistryMetadata(@TempDir Path temporary) throws Exception {
        try (var source = new Source()) {
            MiningCandidateUniverse.build(temporary.resolve("universe"), source, 3, 1, 17,
                    Set.of("CCC"), HASH_A, HASH_B);
        }
        try (var universe = MiningCandidateUniverse.open(temporary.resolve("universe"))) {
            assertEquals(3, universe.candidates().size());
            assertEquals(HASH_B, universe.manifest().fingerprintModelHash);
            for (var candidate : universe.candidates()) {
                assertTrue(candidate.heavyAtoms() > 0);
                assertTrue(candidate.rotatableBonds() >= 0);
                assertTrue(candidate.sp3Fraction() >= 0.0 && candidate.sp3Fraction() <= 1.0);
                assertArrayEquals(vector(candidate.reference().localRow()),
                        universe.fingerprint(candidate.id()), 0.002f);
            }
        }
    }

    private static float[] vector(float value) {
        float[] result = new float[128];
        java.util.Arrays.fill(result, value);
        return result;
    }

    private static final class Source implements MoleculeFingerprintDataSource {
        private static final List<String> SMILES =
                List.of("CC", "CCC", "C[C@H](O)Cl", "c1ccccc1", "CCCCCCCC");
        @Override public String artifactType() { return "test-source"; }
        @Override public long recordCount() { return SMILES.size(); }
        @Override public int shardCount() { return 1; }
        @Override public MoleculeFingerprintProvenance provenance() {
            return new MoleculeFingerprintProvenance(HASH_A, HASH_B, HASH_A);
        }
        @Override public MoleculeFingerprintShardReader openShard(
                int shard, MoleculeFingerprintColumn column) {
            if (shard != 0 || column != MoleculeFingerprintColumn.BASE_128) {
                throw new IllegalArgumentException();
            }
            return new MoleculeFingerprintShardReader() {
                private int next;
                @Override public long recordCount() { return SMILES.size(); }
                @Override public float[] readVector(long row) { return vector(row); }
                @Override public MoleculeVectorBatch readBatch(int maximum) {
                    int count = Math.min(maximum, SMILES.size() - next);
                    float[] values = new float[count * 128];
                    for (int row = 0; row < count; row++) {
                        System.arraycopy(vector(next + row), 0, values, row * 128, 128);
                    }
                    List<MoleculeVectorReference> references = new java.util.ArrayList<>(count);
                    for (int row = 0; row < count; row++) {
                        references.add(new MoleculeVectorReference(0, next + row));
                    }
                    var result = new MoleculeVectorBatch(values, 128, references);
                    next += count;
                    return result;
                }
                @Override public void close() {}
            };
        }
        @Override public Map<MoleculeVectorReference, MoleculeFingerprintMetadata> resolveMetadata(
                java.util.Collection<MoleculeVectorReference> references) {
            Map<MoleculeVectorReference, MoleculeFingerprintMetadata> result = new HashMap<>();
            for (var ref : references) {
                int row = (int) ref.localRow(); String smiles = SMILES.get(row);
                result.put(ref, new MoleculeFingerprintMetadata(row, -1, "m" + row,
                        smiles, smiles));
            }
            return result;
        }
        @Override public void validateCompatibility(DeepSpaceModelBundle model,
                CompactSkelSpheresModelBundle compact, boolean compactRequired) {}
        @Override public void close() {}
    }
}
