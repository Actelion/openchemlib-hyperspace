package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;
import java.util.Objects;

public record ProductTuple(String reactionId, List<String> synthonIds,
                           List<Integer> synthonOrdinals) {
    public ProductTuple {
        Objects.requireNonNull(reactionId);
        synthonIds = List.copyOf(synthonIds);
        synthonOrdinals = List.copyOf(synthonOrdinals);
        if (synthonIds.isEmpty() || synthonIds.size() != synthonOrdinals.size()) {
            throw new IllegalArgumentException("synthon IDs and ordinals must be non-empty and aligned");
        }
    }
}
