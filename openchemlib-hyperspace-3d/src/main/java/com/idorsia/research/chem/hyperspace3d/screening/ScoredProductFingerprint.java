package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import java.util.Arrays;

public record ScoredProductFingerprint(ProductTuple tuple, double objective, float[] components) {
    public ScoredProductFingerprint {
        components = Arrays.copyOf(components, components.length);
    }
}
