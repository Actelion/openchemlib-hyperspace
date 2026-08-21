package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintRecord;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;
import com.idorsia.research.chem.hyperspace3d.model.EmbeddingComparator;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class FingerprintIndexScreener {
    public record Result(List<ScoredProductFingerprint> global,
                         Map<String, List<ScoredProductFingerprint>> byReaction,
                         long uniqueTupleCount) {}

    private final EmbeddingComparator comparator;
    private final DeepSpaceModelManifest manifest;

    public FingerprintIndexScreener(EmbeddingComparator comparator, DeepSpaceModelManifest manifest) {
        this.comparator = comparator;
        this.manifest = manifest;
    }

    public Result screen(ProductFingerprintIndexReader reader, float[] queryEmbedding,
                         ScreeningObjective objective, int batchSize, int globalK,
                         int reactionK) throws IOException {
        GlobalEliteCollector global = new GlobalEliteCollector(globalK);
        ReactionEliteCollector reactions = new ReactionEliteCollector(reactionK);
        Set<ProductTuple> seen = new HashSet<>();
        long count = 0;
        for (;;) {
            List<ProductFingerprintRecord> batch = reader.readBatch(batchSize);
            if (batch.isEmpty()) break;
            List<ProductFingerprintRecord> unique = batch.stream()
                    .filter(record -> seen.add(record.tuple())).toList();
            if (unique.isEmpty()) continue;
            float[][] embeddings = unique.stream().map(ProductFingerprintRecord::embedding)
                    .toArray(float[][]::new);
            float[][] scores = comparator.compare(queryEmbedding, embeddings);
            for (int i = 0; i < unique.size(); i++) {
                ScoredProductFingerprint hit = new ScoredProductFingerprint(unique.get(i).tuple(),
                        objective.score(scores[i], manifest), scores[i]);
                global.offer(hit);
                reactions.offer(hit);
                count++;
            }
        }
        return new Result(global.results(), reactions.results(), count);
    }
}
