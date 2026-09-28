package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.ProductVectorBatch;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorReference;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/** Exact flat scan of compact vectors; tuple payloads are touched only for retained elites. */
public final class CompactSkelSpheresScreener {
    private static final Comparator<CompactSkelSpheresHit> BEST_FIRST =
            Comparator.comparingDouble(CompactSkelSpheresHit::predictedSimilarity).reversed()
                    .thenComparingInt(hit -> hit.reference().shardIndex())
                    .thenComparingLong(hit -> hit.reference().rowIndex());
    private final CompactSkelSpheresScorer scorer;

    public CompactSkelSpheresScreener(CompactSkelSpheresScorer scorer) { this.scorer = scorer; }

    public CompactSkelSpheresScreeningResult screen(Path directory, ProductVectorIndexManifest manifest,
            float[] query, int globalK, int perReactionK, int batchSize) throws IOException {
        if (query.length != 16 || globalK < 1 || perReactionK < 0 || batchSize < 1)
            throw new IllegalArgumentException("invalid screening arguments");
        PriorityQueue<CompactSkelSpheresHit> global = new PriorityQueue<>(BEST_FIRST.reversed());
        Map<Integer, PriorityQueue<CompactSkelSpheresHit>> reactions = new HashMap<>();
        long records = 0; long scanStarted = System.nanoTime();
        for (int shard = 0; shard < manifest.shards.size(); shard++) {
            try (ProductVectorIndexReader reader = new ProductVectorIndexReader(directory, manifest, shard)) {
                while (true) {
                    ProductVectorBatch batch = reader.readBatch(batchSize);
                    if (batch.vectors().length == 0) break;
                    for (int i = 0; i < batch.vectors().length; i++) {
                        ProductVectorReference ref = batch.references().get(i);
                        var score = scorer.score(query, batch.vectors()[i]);
                        var hit = new CompactSkelSpheresHit(ref, null, score.dotProduct(),
                                score.calibratedSimilarity(), null);
                        retain(global, hit, globalK);
                        if (perReactionK > 0) retain(reactions.computeIfAbsent(
                                ref.reactionOrdinal(), ignored -> new PriorityQueue<>(BEST_FIRST.reversed())),
                                hit, perReactionK);
                        records++;
                    }
                }
            }
        }
        long scanNanos = System.nanoTime() - scanStarted;
        List<CompactSkelSpheresHit> globalHits = sorted(global);
        Map<Integer, List<CompactSkelSpheresHit>> ordinalHits = new HashMap<>();
        reactions.forEach((ordinal, queue) -> ordinalHits.put(ordinal, sorted(queue)));
        long resolveStarted = System.nanoTime();
        List<CompactSkelSpheresHit> finalists = new ArrayList<>(globalHits);
        for (List<CompactSkelSpheresHit> hits : ordinalHits.values()) finalists.addAll(hits);
        Map<ProductVectorReference, CompactSkelSpheresHit> resolved =
                resolveAll(directory, manifest, finalists);
        globalHits = globalHits.stream().map(hit -> resolved.get(hit.reference())).toList();
        Map<String, List<CompactSkelSpheresHit>> byReaction = new LinkedHashMap<>();
        ordinalHits.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(entry ->
                byReaction.put(manifest.reactionDictionary.get(entry.getKey()),
                        entry.getValue().stream().map(hit -> resolved.get(hit.reference())).toList()));
        return new CompactSkelSpheresScreeningResult(globalHits, byReaction, records, scanNanos,
                System.nanoTime() - resolveStarted);
    }

    private static Map<ProductVectorReference, CompactSkelSpheresHit> resolveAll(
            Path directory, ProductVectorIndexManifest manifest,
            List<CompactSkelSpheresHit> finalists) throws IOException {
        Map<Integer, Map<ProductVectorReference, CompactSkelSpheresHit>> byShard = new HashMap<>();
        for (CompactSkelSpheresHit hit : finalists) {
            byShard.computeIfAbsent(hit.reference().shardIndex(), ignored -> new HashMap<>())
                    .putIfAbsent(hit.reference(), hit);
        }
        Map<ProductVectorReference, CompactSkelSpheresHit> resolved = new HashMap<>();
        for (var entry : byShard.entrySet()) {
            try (ProductVectorIndexReader reader = new ProductVectorIndexReader(
                    directory, manifest, entry.getKey())) {
                for (CompactSkelSpheresHit hit : entry.getValue().values())
                    resolved.put(hit.reference(), hit.withTuple(reader.resolveTuple(hit.reference())));
            }
        }
        return resolved;
    }
    private static void retain(PriorityQueue<CompactSkelSpheresHit> queue,
            CompactSkelSpheresHit hit, int limit) {
        if (queue.size() < limit) queue.add(hit);
        else if (BEST_FIRST.compare(hit, queue.peek()) < 0) { queue.poll(); queue.add(hit); }
    }
    private static List<CompactSkelSpheresHit> sorted(PriorityQueue<CompactSkelSpheresHit> queue) {
        List<CompactSkelSpheresHit> result = new ArrayList<>(queue); result.sort(BEST_FIRST); return result;
    }
}
