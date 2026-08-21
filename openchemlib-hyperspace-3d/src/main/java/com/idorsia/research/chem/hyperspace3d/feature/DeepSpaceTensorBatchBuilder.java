package com.idorsia.research.chem.hyperspace3d.feature;

import com.actelion.research.chem.StereoMolecule;
import java.util.List;

public final class DeepSpaceTensorBatchBuilder {
    private final OCLDeepSpaceFeaturizer featurizer;

    public DeepSpaceTensorBatchBuilder(OCLDeepSpaceFeaturizer featurizer) {
        this.featurizer = featurizer;
    }

    public DeepSpaceTensorBatch build(List<StereoMolecule> molecules) {
        if (molecules.isEmpty()) throw new IllegalArgumentException("batch must not be empty");
        int atoms = DeepSpaceFeatureSchema.MAX_ATOMS;
        float[] atom = new float[molecules.size() * atoms * DeepSpaceFeatureSchema.ATOM_FEATURE_DIM];
        float[] pair = new float[molecules.size() * atoms * atoms * DeepSpaceFeatureSchema.PAIR_FEATURE_DIM];
        boolean[] mask = new boolean[molecules.size() * atoms];
        for (int i = 0; i < molecules.size(); i++) {
            FeaturizationResult result = featurizer.featurize(molecules.get(i));
            if (!result.accepted()) {
                throw new IllegalArgumentException("molecule " + i + " rejected: " + result.rejectionReason());
            }
            System.arraycopy(result.atomFeatures(), 0, atom,
                    i * atoms * DeepSpaceFeatureSchema.ATOM_FEATURE_DIM, result.atomFeatures().length);
            System.arraycopy(result.pairFeatures(), 0, pair,
                    i * atoms * atoms * DeepSpaceFeatureSchema.PAIR_FEATURE_DIM, result.pairFeatures().length);
            System.arraycopy(result.atomMask(), 0, mask, i * atoms, result.atomMask().length);
        }
        return new DeepSpaceTensorBatch(molecules.size(), atom, pair, mask);
    }
}
