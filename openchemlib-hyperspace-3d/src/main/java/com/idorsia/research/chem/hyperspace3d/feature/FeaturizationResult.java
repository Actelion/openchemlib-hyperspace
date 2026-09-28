package com.idorsia.research.chem.hyperspace3d.feature;

import java.util.Objects;
import java.util.Optional;

public record FeaturizationResult(float[] atomFeatures, float[] pairFeatures,
                                  boolean[] atomMask, String rejectionReason) {
    public FeaturizationResult {
        if (rejectionReason == null) {
            Objects.requireNonNull(atomFeatures);
            Objects.requireNonNull(pairFeatures);
            Objects.requireNonNull(atomMask);
        }
    }

    public static FeaturizationResult rejected(String reason) {
        return new FeaturizationResult(null, null, null, Objects.requireNonNull(reason));
    }

    public boolean accepted() { return rejectionReason == null; }
    public Optional<String> rejection() { return Optional.ofNullable(rejectionReason); }
}
