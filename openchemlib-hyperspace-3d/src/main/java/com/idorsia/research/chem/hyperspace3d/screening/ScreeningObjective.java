package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelManifest;

@FunctionalInterface
public interface ScreeningObjective {
    double score(float[] scores, DeepSpaceModelManifest manifest);

    static ScreeningObjective direct(String target) {
        return (scores, manifest) -> scores[manifest.targetIndex(target)];
    }

    static ScreeningObjective phesaComposite(double shapeWeight, double pharmacophoreWeight) {
        if (!Double.isFinite(shapeWeight) || !Double.isFinite(pharmacophoreWeight)) {
            throw new IllegalArgumentException("composite weights must be finite");
        }
        return (scores, manifest) ->
                shapeWeight * scores[manifest.targetIndex("phesa_shape")]
                + pharmacophoreWeight * scores[manifest.targetIndex("phesa_pharmacophore")];
    }
}
