package com.idorsia.research.chem.hyperspace3d.feature;

import com.actelion.research.chem.StereoMolecule;
import java.util.ArrayList;
import java.util.List;

public final class DeepSpaceTensorBatchBuilder {
    private final OCLDeepSpaceFeaturizer featurizer;

    public DeepSpaceTensorBatchBuilder(OCLDeepSpaceFeaturizer featurizer) {
        this.featurizer = featurizer;
    }

    public DeepSpaceTensorBatch build(List<StereoMolecule> molecules) {
        if (molecules.isEmpty()) throw new IllegalArgumentException("batch must not be empty");
        List<FeaturizationResult> features = new ArrayList<>(molecules.size());
        for (StereoMolecule molecule : molecules) features.add(featurizer.featurize(molecule));
        return buildFromFeatures(features);
    }

    /** Packs already calculated features without running the featurizer again. */
    public DeepSpaceTensorBatch buildFromFeatures(List<FeaturizationResult> features) {
        if (features.isEmpty()) throw new IllegalArgumentException("batch must not be empty");
        int atoms = DeepSpaceFeatureSchema.MAX_ATOMS;
        float[] atom = new float[features.size() * atoms * DeepSpaceFeatureSchema.ATOM_FEATURE_DIM];
        float[] pair = new float[features.size() * atoms * atoms * DeepSpaceFeatureSchema.PAIR_FEATURE_DIM];
        boolean[] mask = new boolean[features.size() * atoms];
        for (int i = 0; i < features.size(); i++) {
            FeaturizationResult result = features.get(i);
            if (!result.accepted()) {
                throw new IllegalArgumentException("molecule " + i + " rejected: " + result.rejectionReason());
            }
            System.arraycopy(result.atomFeatures(), 0, atom,
                    i * atoms * DeepSpaceFeatureSchema.ATOM_FEATURE_DIM, result.atomFeatures().length);
            System.arraycopy(result.pairFeatures(), 0, pair,
                    i * atoms * atoms * DeepSpaceFeatureSchema.PAIR_FEATURE_DIM, result.pairFeatures().length);
            System.arraycopy(result.atomMask(), 0, mask, i * atoms, result.atomMask().length);
        }
        return new DeepSpaceTensorBatch(features.size(), atom, pair, mask);
    }
}
