package com.idorsia.research.chem.hyperspace3d.batch;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.FeaturizationResult;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Comparator;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import com.idorsia.research.chem.hyperspace3d.screening.ScreeningObjective;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Correctness-first vertical slice. A production coordinator may feed this
 * scorer through bounded CPU/GPU queues without changing its batch contract.
 */
public final class DeepSpaceBatchAssemblyScorer implements BatchAssemblyScorer {
    private final SynthonTupleResolver resolver;
    private final OCLDeepSpaceFeaturizer featurizer;
    private final DeepSpaceTensorBatchBuilder batchBuilder;
    private final DeepSpaceV1Encoder encoder;
    private final DeepSpaceV1Comparator comparator;
    private final DeepSpaceModelManifest manifest;
    private final ScreeningObjective objective;
    private final float[] queryEmbedding;
    private final Map<ProductTuple, float[]> runEmbeddingCache = new ConcurrentHashMap<>();

    public DeepSpaceBatchAssemblyScorer(SynthonTupleResolver resolver,
            OCLDeepSpaceFeaturizer featurizer, DeepSpaceV1Encoder encoder,
            DeepSpaceV1Comparator comparator, DeepSpaceModelManifest manifest,
            ScreeningObjective objective, float[] queryEmbedding) {
        this.resolver = resolver;
        this.featurizer = featurizer;
        this.batchBuilder = new DeepSpaceTensorBatchBuilder(featurizer);
        this.encoder = encoder;
        this.comparator = comparator;
        this.manifest = manifest;
        this.objective = objective;
        this.queryEmbedding = queryEmbedding.clone();
    }

    @Override
    public List<ScoredAssemblyCandidate> scoreBatch(List<AssemblyCandidate> candidates) {
        Map<ProductTuple, List<AssemblyCandidate>> deduplicated = new LinkedHashMap<>();
        candidates.forEach(candidate -> deduplicated.computeIfAbsent(candidate.tuple(),
                ignored -> new ArrayList<>()).add(candidate));
        Map<ProductTuple, float[]> embeddings = new LinkedHashMap<>();
        Map<ProductTuple, String> rejected = new LinkedHashMap<>();
        List<ProductTuple> pendingTuples = new ArrayList<>();
        List<StereoMolecule> pendingMolecules = new ArrayList<>();
        deduplicated.forEach((tuple, copies) -> {
            float[] cached = runEmbeddingCache.get(tuple);
            if (cached != null) {
                embeddings.put(tuple, cached);
                return;
            }
            try {
                StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(resolver.resolve(tuple));
                FeaturizationResult features = featurizer.featurize(assembled);
                if (!features.accepted()) rejected.put(tuple, features.rejectionReason());
                else {
                    pendingTuples.add(tuple);
                    pendingMolecules.add(assembled);
                }
            } catch (RuntimeException failure) {
                rejected.put(tuple, "ASSEMBLY_FAILED:" + failure.getClass().getSimpleName());
            }
        });
        if (!pendingMolecules.isEmpty()) {
            float[][] encoded = encoder.encode(batchBuilder.build(pendingMolecules));
            for (int i = 0; i < encoded.length; i++) {
                runEmbeddingCache.put(pendingTuples.get(i), encoded[i]);
                embeddings.put(pendingTuples.get(i), encoded[i]);
            }
        }
        List<ProductTuple> scoredTuples = embeddings.keySet().stream().toList();
        float[][] scores = scoredTuples.isEmpty() ? new float[0][]
                : comparator.compare(queryEmbedding,
                        scoredTuples.stream().map(embeddings::get).toArray(float[][]::new));
        Map<ProductTuple, float[]> scoreByTuple = new LinkedHashMap<>();
        for (int i = 0; i < scoredTuples.size(); i++) scoreByTuple.put(scoredTuples.get(i), scores[i]);

        List<ScoredAssemblyCandidate> result = new ArrayList<>(candidates.size());
        for (AssemblyCandidate candidate : candidates) {
            float[] components = scoreByTuple.get(candidate.tuple());
            if (components == null) {
                result.add(ScoredAssemblyCandidate.rejected(candidate,
                        rejected.getOrDefault(candidate.tuple(), "SCORING_FAILED")));
            } else {
                result.add(new ScoredAssemblyCandidate(candidate,
                        ScoredAssemblyCandidate.Status.SCORED, null,
                        objective.score(components, manifest), components));
            }
        }
        return result;
    }

    public int cachedTupleCount() { return runEmbeddingCache.size(); }
}
