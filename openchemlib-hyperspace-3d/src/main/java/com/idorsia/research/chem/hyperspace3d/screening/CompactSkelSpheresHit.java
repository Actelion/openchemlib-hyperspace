package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorReference;

public record CompactSkelSpheresHit(ProductVectorReference reference, ProductTuple tuple,
        double dotProduct, double predictedSimilarity, Double exactSimilarity) {
    public CompactSkelSpheresHit withTuple(ProductTuple value) {
        return new CompactSkelSpheresHit(reference, value, dotProduct, predictedSimilarity, exactSimilarity);
    }
    public CompactSkelSpheresHit withExact(double value) {
        return new CompactSkelSpheresHit(reference, tuple, dotProduct, predictedSimilarity, value);
    }
}
