package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;

public record ProductVectorBatch(List<ProductVectorReference> references, float[][] vectors) {
    public ProductVectorBatch {
        references = List.copyOf(references);
        if (references.size() != vectors.length) throw new IllegalArgumentException("batch size mismatch");
    }
}
