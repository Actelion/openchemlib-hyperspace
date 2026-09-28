package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorBatchPrefetcher;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;
import com.idorsia.research.chem.hyperspace3d.model.EmbeddingComparator;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/** Streaming learned PheSA screening and compact-shortlist reranking. */
public final class MoleculePheSAScreener {
    private static final Comparator<MoleculePheSAHit> BEST_FIRST =
            Comparator.comparingDouble(MoleculePheSAHit::objectiveScore).reversed()
                    .thenComparingInt(hit -> hit.reference().shardIndex())
                    .thenComparingLong(hit -> hit.reference().localRow());
    private final EmbeddingComparator comparator;
    private final DeepSpaceModelManifest manifest;
    private final ScreeningObjective objective;

    public MoleculePheSAScreener(EmbeddingComparator comparator,
            DeepSpaceModelManifest manifest, ScreeningObjective objective) {
        this.comparator = comparator;
        this.manifest = manifest;
        this.objective = objective;
    }

    public MoleculePheSAScreeningResult screen(MoleculeFingerprintDataSource source,
            float[] queryEmbedding, int topK, int batchSize) throws IOException {
        return screen(source, queryEmbedding, topK, batchSize, 2);
    }

    public MoleculePheSAScreeningResult screen(MoleculeFingerprintDataSource source,
            float[] queryEmbedding, int topK, int batchSize, int prefetchDepth) throws IOException {
        validate(queryEmbedding, topK, batchSize);
        if (prefetchDepth < 1) throw new IllegalArgumentException("prefetch depth must be positive");
        PriorityQueue<MoleculePheSAHit> retained = new PriorityQueue<>(BEST_FIRST.reversed());
        long records = 0;
        long comparatorNanos = 0;
        long started = System.nanoTime();
        try (var batches = new MoleculeVectorBatchPrefetcher(source,
                MoleculeFingerprintColumn.BASE_128, batchSize, prefetchDepth)) {
            for (var batch = batches.take(); batch != null; batch = batches.take()) {
                long inferenceStarted = System.nanoTime();
                float[][] scores = comparator.compareFlat(
                        queryEmbedding, batch.values(), batch.size());
                comparatorNanos += System.nanoTime() - inferenceStarted;
                for (int row = 0; row < batch.size(); row++) {
                    retain(retained, batch.references().get(row), scores[row], null, topK);
                }
                records += batch.size();
            }
        }
        long vectorAndInferenceNanos = System.nanoTime() - started;
        return resolve(source, retained, records, vectorAndInferenceNanos, comparatorNanos);
    }

    public MoleculePheSAScreeningResult rerank(MoleculeFingerprintDataSource source,
            float[] queryEmbedding, List<MoleculeSkelSpheresHit> compactHits,
            int topK, int batchSize) throws IOException {
        validate(queryEmbedding, topK, batchSize);
        if (compactHits == null) throw new IllegalArgumentException("compact hits are required");
        Map<MoleculeVectorReference, CompactProvenance> compact = new HashMap<>();
        Map<Integer, List<MoleculeVectorReference>> byShard = new HashMap<>();
        for (int rank = 0; rank < compactHits.size(); rank++) {
            var hit = compactHits.get(rank);
            compact.put(hit.reference(), new CompactProvenance(
                    hit.predictedSimilarity(), hit.dotProduct(), rank + 1));
            byShard.computeIfAbsent(hit.reference().shardIndex(), ignored -> new ArrayList<>())
                    .add(hit.reference());
        }
        for (var references : byShard.values()) {
            references.sort(Comparator.comparingLong(MoleculeVectorReference::localRow));
        }

        PriorityQueue<MoleculePheSAHit> retained = new PriorityQueue<>(BEST_FIRST.reversed());
        long comparatorNanos = 0;
        long started = System.nanoTime();
        for (int shard = 0; shard < source.shardCount(); shard++) {
            List<MoleculeVectorReference> references = byShard.get(shard);
            if (references == null) continue;
            try (var reader = source.openShard(shard, MoleculeFingerprintColumn.BASE_128)) {
                for (int start = 0; start < references.size(); start += batchSize) {
                    int count = Math.min(batchSize, references.size() - start);
                    float[] vectors = new float[Math.multiplyExact(count, 128)];
                    for (int row = 0; row < count; row++) {
                        float[] vector = reader.readVector(references.get(start + row).localRow());
                        System.arraycopy(vector, 0, vectors, row * 128, 128);
                    }
                    long inferenceStarted = System.nanoTime();
                    float[][] scores = comparator.compareFlat(queryEmbedding, vectors, count);
                    comparatorNanos += System.nanoTime() - inferenceStarted;
                    for (int row = 0; row < count; row++) {
                        var reference = references.get(start + row);
                        var provenance = compact.get(reference);
                        retain(retained, reference, scores[row], provenance, topK);
                    }
                }
            }
        }
        long vectorAndInferenceNanos = System.nanoTime() - started;
        return resolve(source, retained, compactHits.size(),
                vectorAndInferenceNanos, comparatorNanos);
    }


    private static MoleculePheSAScreeningResult resolve(MoleculeFingerprintDataSource source,
            PriorityQueue<MoleculePheSAHit> retained, long records,
            long vectorAndInferenceNanos, long comparatorNanos) throws IOException {
        List<MoleculePheSAHit> hits = new ArrayList<>(retained);
        hits.sort(BEST_FIRST);
        long metadataStarted = System.nanoTime();
        var metadata = source.resolveMetadata(hits.stream()
                .map(MoleculePheSAHit::reference).toList());
        hits = hits.stream().map(hit -> hit.withMolecule(metadata.get(hit.reference()))).toList();
        return new MoleculePheSAScreeningResult(hits, records, vectorAndInferenceNanos,
                comparatorNanos, System.nanoTime() - metadataStarted);
    }

    private void retain(PriorityQueue<MoleculePheSAHit> queue,
            MoleculeVectorReference reference, float[] scores,
            CompactProvenance compact, int limit) {
        double objectiveScore = objective.score(scores, manifest);
        MoleculePheSAHit worst = queue.peek();
        boolean accepted = queue.size() < limit || objectiveScore > worst.objectiveScore()
                || (Double.compare(objectiveScore, worst.objectiveScore()) == 0
                && (reference.shardIndex() < worst.reference().shardIndex()
                || (reference.shardIndex() == worst.reference().shardIndex()
                && reference.localRow() < worst.reference().localRow())));
        if (!accepted) return;
        if (queue.size() == limit) queue.poll();
        queue.add(new MoleculePheSAHit(reference, null, objectiveScore, scores,
                compact == null ? null : compact.predicted,
                compact == null ? null : compact.dot,
                compact == null ? null : compact.rank, null));
    }

    private static void validate(float[] query, int topK, int batchSize) {
        if (query == null || query.length != 128 || topK < 1 || batchSize < 1) {
            throw new IllegalArgumentException("invalid learned PheSA screening arguments");
        }
    }

    private record CompactProvenance(double predicted, double dot, int rank) {}
}
