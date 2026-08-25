package com.idorsia.research.chem.hyperspace3d.model;

import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;

/** Encodes one tensor batch into the two persistent molecule representations. */
@FunctionalInterface
public interface MoleculeFingerprintBatchEncoder {
    MoleculeFingerprintBatch encode(DeepSpaceTensorBatch batch);
}
