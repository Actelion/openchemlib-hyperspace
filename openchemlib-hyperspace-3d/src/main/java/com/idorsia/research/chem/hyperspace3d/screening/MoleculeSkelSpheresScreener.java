package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.JavaMoleculeFingerprintDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexManifest;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorBatchPrefetcher;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
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

    /** Backward-compatible entry point for native Java indexes. */
    public MoleculeSkelSpheresScreeningResult screen(Path directory,
            MoleculeFingerprintIndexManifest ignoredManifest, float[] query,
            int topK, int batchSize) throws IOException {
        try (var source = new JavaMoleculeFingerprintDataSource(directory)) {
            return screen(source, query, topK, batchSize);
        }
    }

    public MoleculeSkelSpheresScreeningResult screen(MoleculeFingerprintDataSource source,
            float[] query, int topK, int batchSize) throws IOException {
        return screen(source, query, topK, batchSize, 2);
    }

    public MoleculeSkelSpheresScreeningResult screen(MoleculeFingerprintDataSource source,
            float[] query, int topK, int batchSize, int prefetchDepth) throws IOException {
        var unresolved = scanUnresolved(source, query, topK, batchSize, prefetchDepth);
        long resolveStarted = System.nanoTime();
        var metadata = source.resolveMetadata(unresolved.hits().stream()
                .map(MoleculeSkelSpheresHit::reference).toList());
        var hits = unresolved.hits().stream()
                .map(hit -> hit.withMolecule(metadata.get(hit.reference()))).toList();
        return new MoleculeSkelSpheresScreeningResult(hits, unresolved.recordsScanned(),
                unresolved.scanNanos(), System.nanoTime() - resolveStarted);
    }

    public MoleculeSkelSpheresScreeningResult scanUnresolved(
            MoleculeFingerprintDataSource source, float[] query,
            int topK, int batchSize) throws IOException {
        return scanUnresolved(source, query, topK, batchSize, 2);
    }

    /** Compact scan that deliberately leaves metadata unresolved for a later cascade stage. */
    public MoleculeSkelSpheresScreeningResult scanUnresolved(
            MoleculeFingerprintDataSource source, float[] query,
            int topK, int batchSize, int prefetchDepth) throws IOException {
        if (query == null || query.length != 16 || topK < 1 || batchSize < 1
                || prefetchDepth < 1) {
            throw new IllegalArgumentException("invalid molecule screening arguments");
        }
        PriorityQueue<MoleculeSkelSpheresHit> retained =
                new PriorityQueue<>(BEST_FIRST.reversed());
        long records = 0;
        long scanStarted = System.nanoTime();
        try (var batches = new MoleculeVectorBatchPrefetcher(source,
                MoleculeFingerprintColumn.COMPACT_16, batchSize, prefetchDepth)) {
            for (var batch = batches.take(); batch != null; batch = batches.take()) {
                for (int row = 0; row < batch.size(); row++) {
                    int offset = row * 16;
                    double dot = 0;
                    for (int column = 0; column < 16; column++) {
                        dot += (double) query[column] * batch.values()[offset + column];
                    }
                    retain(retained, batch.references().get(row), dot,
                            scorer.calibrated(dot), topK);
                    records++;
                }
            }
        }
        List<MoleculeSkelSpheresHit> hits = new ArrayList<>(retained);
        hits.sort(BEST_FIRST);
        return new MoleculeSkelSpheresScreeningResult(hits, records,
                System.nanoTime() - scanStarted, 0);
    }

    private static void retain(PriorityQueue<MoleculeSkelSpheresHit> queue,
            MoleculeVectorReference reference,
            double dot, double predicted, int limit) {
        MoleculeSkelSpheresHit worst = queue.peek();
        boolean accepted = queue.size() < limit || predicted > worst.predictedSimilarity()
                || (Double.compare(predicted, worst.predictedSimilarity()) == 0
                && (reference.shardIndex() < worst.reference().shardIndex()
                || (reference.shardIndex() == worst.reference().shardIndex()
                && reference.localRow() < worst.reference().localRow())));
        if (!accepted) return;
        if (queue.size() == limit) queue.poll();
        queue.add(new MoleculeSkelSpheresHit(reference, null, dot, predicted, null));
    }
}
