package com.idorsia.research.chem.hyperspace3d.screening;

import java.util.List;

public record MoleculePheSAScreeningResult(List<MoleculePheSAHit> hits,
        long recordsScanned, long vectorAndInferenceNanos, long comparatorNanos,
        long metadataResolutionNanos) {
    public long vectorPreparationNanos() { return vectorAndInferenceNanos - comparatorNanos; }
}
