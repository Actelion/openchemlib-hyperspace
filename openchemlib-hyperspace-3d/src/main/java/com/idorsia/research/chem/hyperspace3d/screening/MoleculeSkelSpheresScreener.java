package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/** Exact flat scan of molecule compact vectors with delayed metadata resolution. */
public final class MoleculeSkelSpheresScreener {
    private static final Comparator<MoleculeSkelSpheresHit> BEST_FIRST =
            Comparator.comparingDouble(MoleculeSkelSpheresHit::predictedSimilarity).reversed()
                    .thenComparingInt(hit -> hit.reference().shardIndex())
                    .thenComparingLong(hit -> hit.reference().localRow());
    private final CompactSkelSpheresScorer scorer;

    public MoleculeSkelSpheresScreener(CompactSkelSpheresScorer scorer) {
        this.scorer = scorer;
    }

    public MoleculeSkelSpheresScreeningResult screen(Path directory,
            MoleculeFingerprintIndexManifest manifest, float[] query,
            int topK, int batchSize) throws IOException {
        if (query == null || query.length != 16 || topK < 1 || batchSize < 1) {
            throw new IllegalArgumentException("invalid molecule screening arguments");
        }
        PriorityQueue<MoleculeSkelSpheresHit> retained =
                new PriorityQueue<>(BEST_FIRST.reversed());
        long records = 0;
        long scanStarted = System.nanoTime();
        for (int shard = 0; shard < manifest.shards.size(); shard++) {
            try (var reader = new MoleculeFingerprintIndexReader(directory, manifest, shard,
                    MoleculeFingerprintIndexReader.VectorColumns.COMPACT_16)) {
                while (true) {
                    var batch = reader.readCompactBatch(batchSize);
                    if (batch.vectors().length == 0) break;
                    for (int row = 0; row < batch.vectors().length; row++) {
                        var score = scorer.score(query, batch.vectors()[row]);
                        retain(retained, new MoleculeSkelSpheresHit(
                                batch.references().get(row), null, score.dotProduct(),
                                score.calibratedSimilarity(), null), topK);
                        records++;
                    }
                }
            }
        }
        long scanNanos = System.nanoTime() - scanStarted;
        List<MoleculeSkelSpheresHit> hits = new ArrayList<>(retained);
        hits.sort(BEST_FIRST);
        long resolveStarted = System.nanoTime();
        Map<MoleculeVectorReference, MoleculeSkelSpheresHit> resolved =
                resolve(directory, manifest, hits);
        hits = hits.stream().map(hit -> resolved.get(hit.reference())).toList();
        return new MoleculeSkelSpheresScreeningResult(hits, records, scanNanos,
                System.nanoTime() - resolveStarted);
    }

    private static Map<MoleculeVectorReference, MoleculeSkelSpheresHit> resolve(
            Path directory, MoleculeFingerprintIndexManifest manifest,
            List<MoleculeSkelSpheresHit> hits) throws IOException {
        Map<Integer, List<MoleculeSkelSpheresHit>> byShard = new HashMap<>();
        for (var hit : hits) byShard.computeIfAbsent(hit.reference().shardIndex(),
                ignored -> new ArrayList<>()).add(hit);
        Map<MoleculeVectorReference, MoleculeSkelSpheresHit> result = new HashMap<>();
        for (var entry : byShard.entrySet()) {
            try (var reader = new MoleculeFingerprintIndexReader(directory, manifest,
                    entry.getKey(), MoleculeFingerprintIndexReader.VectorColumns.COMPACT_16)) {
                for (var hit : entry.getValue()) {
                    result.put(hit.reference(), hit.withMolecule(
                            reader.readMetadata(hit.reference().localRow())));
                }
            }
        }
        return result;
    }

    private static void retain(PriorityQueue<MoleculeSkelSpheresHit> queue,
            MoleculeSkelSpheresHit hit, int limit) {
        if (queue.size() < limit) queue.add(hit);
        else if (BEST_FIRST.compare(hit, queue.peek()) < 0) {
            queue.poll();
            queue.add(hit);
        }
    }
}
