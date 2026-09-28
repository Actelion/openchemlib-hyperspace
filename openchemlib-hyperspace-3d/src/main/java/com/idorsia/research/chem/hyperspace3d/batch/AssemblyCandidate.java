package com.idorsia.research.chem.hyperspace3d.batch;

import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import java.util.Map;
import java.util.Objects;

public record AssemblyCandidate(long candidateId, ProductTuple tuple, String basinId,
                                int beamRound, Map<String, String> provenance) {
    public AssemblyCandidate {
        Objects.requireNonNull(tuple);
        Objects.requireNonNull(basinId);
        provenance = provenance == null ? Map.of() : Map.copyOf(provenance);
        if (beamRound < 0) throw new IllegalArgumentException("beamRound must be non-negative");
    }
}
