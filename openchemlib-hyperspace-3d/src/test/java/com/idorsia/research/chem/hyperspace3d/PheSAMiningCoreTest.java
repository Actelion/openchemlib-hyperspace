package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace3d.mining.ExactPheSALabeler;
import com.idorsia.research.chem.hyperspace3d.mining.MiningSelection;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAQueryPairRecord;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAQueryShardIO;
import com.idorsia.research.chem.hyperspace3d.mining.PredictedRankMiner;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import static org.junit.jupiter.api.Assertions.*;

class PheSAMiningCoreTest {
    @TempDir Path temporary;

    @Test void queryShardRoundTripPreservesTheCrossLanguageRecord() throws Exception {
        var record = new PheSAQueryPairRecord(42, 0.7f, 0.8f, 0.6f,
                1, 2, 3, 0.75f, 0.8f, 0.7f,
                MiningSelection.RANDOM | MiningSelection.rankBucket(
                        MiningSelection.TOTAL_OFFSET, 0), 2, PheSAQueryPairRecord.EXACT_OK);
        Path path = temporary.resolve("query.bin");
        PheSAQueryShardIO.write(path, 7, 1_000_000, 2, List.of(record));
        var loaded = PheSAQueryShardIO.read(path);
        assertEquals(7, loaded.queryId());
        assertEquals(1_000_000, loaded.screenedCount());
        assertEquals(record, loaded.records().get(0));
        assertEquals(PheSAQueryShardIO.HEADER_BYTES + PheSAQueryShardIO.RECORD_BYTES,
                java.nio.file.Files.size(path));
    }

    @Test void independentRanksAreStableAndSelectionsOverlapByMask() {
        float[][] scores = new float[10_000][3];
        for (int row = 0; row < scores.length; row++) {
            scores[row][0] = row / 10_000f;
            scores[row][1] = (scores.length - row) / 10_000f;
            scores[row][2] = row % 11;
        }
        var miner = new PredictedRankMiner(
                PredictedRankMiner.DEFAULT_BOUNDARIES,
                new int[]{2, 3, 4, 5, 6, 7}, 10, 17);
        var first = miner.select(scores);
        var second = miner.select(scores);
        assertEquals(first.stream().map(PredictedRankMiner.Selection::candidateOrdinal).toList(),
                second.stream().map(PredictedRankMiner.Selection::candidateOrdinal).toList());
        assertTrue(first.stream().anyMatch(value ->
                (value.selectionMask() & MiningSelection.RANDOM) != 0));
        assertTrue(first.stream().allMatch(value -> value.rankTotal() >= 1
                && value.rankTotal() <= scores.length));
        assertEquals(1, first.stream().filter(value -> value.candidateOrdinal() == 9999)
                .findFirst().orElseThrow().rankTotal());
        assertEquals(1, first.stream().filter(value -> value.candidateOrdinal() == 0)
                .findFirst().orElseThrow().rankShape());
    }

    @Test void exactPheSAExposesAnalyticallyConsistentComponents() throws Exception {
        var parser = new SmilesParser();
        StereoMolecule molecule = new StereoMolecule();
        parser.parse(molecule, "c1ccccc1");
        var labeler = new ExactPheSALabeler(1, 0.5);
        var descriptor = labeler.describe(molecule);
        assertNotNull(descriptor);
        var label = labeler.score(descriptor, descriptor);
        assertEquals(PheSAQueryPairRecord.EXACT_OK, label.status());
        assertEquals(0.5 * (label.shape() + label.pharmacophore()), label.total(), 2e-6);
        assertEquals(1.0, label.total(), 1e-4);
    }
}
