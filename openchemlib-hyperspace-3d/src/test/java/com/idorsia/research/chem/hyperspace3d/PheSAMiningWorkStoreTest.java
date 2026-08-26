package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace3d.mining.ExactPheSALabeler;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAMiningWorkStore;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAQueryPairRecord;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PheSAMiningWorkStoreTest {
    @Test
    void roundTripsAndValidatesRestartState(@TempDir Path temporary) throws Exception {
        String identity = "a".repeat(64); Path work = temporary.resolve("work");
        var store = new PheSAMiningWorkStore(work, identity, 2);
        var record = new PheSAQueryPairRecord(2, .5f, .6f, .4f, 1, 2, 3,
                .5f, .6f, .4f, 7, 0, PheSAQueryPairRecord.EXACT_OK);
        var descriptors = List.of(
                new PheSAMiningWorkStore.DescriptorEntry(10, 2, "encoded-value"),
                new PheSAMiningWorkStore.DescriptorEntry(11, 3, null));
        store.save(0, 0, 1_000_000, 0, List.of(record), descriptors);

        var resumed = new PheSAMiningWorkStore(work, identity, 2).load(0);
        assertEquals(List.of(record), resumed.shard().records());
        assertEquals("encoded-value", resumed.descriptorEntries().get(0).descriptor());
        assertNull(resumed.descriptorEntries().get(1).descriptor());
        assertNull(store.load(1));
        assertThrows(java.io.IOException.class,
                () -> new PheSAMiningWorkStore(work, "b".repeat(64), 2));

        Files.write(work.resolve("queries/q000000.bin"), new byte[]{1},
                StandardOpenOption.APPEND);
        assertThrows(java.io.IOException.class, () -> store.load(0));
    }

    @Test
    void threadConfinedDescriptorGenerationIsDeterministic() throws Exception {
        List<String> smiles = List.of("CCO", "c1ccccc1", "C[C@H](O)Cl", "CCN(CC)CC",
                "O=C(O)c1ccccc1", "C1CCCCC1", "FC(F)(F)c1ccccc1", "CC(=O)NC");
        List<String> sequential = new ArrayList<>();
        for (String value : smiles) sequential.add(descriptor(value));
        var executor = Executors.newFixedThreadPool(4);
        try {
            var futures = smiles.stream().map(value -> executor.submit(() -> descriptor(value))).toList();
            for (int index = 0; index < smiles.size(); index++) {
                assertEquals(sequential.get(index), futures.get(index).get());
            }
        } finally { executor.shutdownNow(); }
    }

    private static String descriptor(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule(); new SmilesParser().parse(molecule, smiles);
        var labeler = new ExactPheSALabeler(8, .5);
        return labeler.encode(labeler.describe(molecule));
    }
}
