package com.idorsia.research.chem.hyperspace3d.model;

import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;

/** Query-independent molecular embedding boundary used by index construction. */
@FunctionalInterface
public interface ProductEmbeddingBatchEncoder {
    float[][] encode(DeepSpaceTensorBatch batch);
}
