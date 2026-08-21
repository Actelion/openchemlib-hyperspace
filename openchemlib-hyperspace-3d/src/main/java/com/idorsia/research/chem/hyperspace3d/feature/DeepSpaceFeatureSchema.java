package com.idorsia.research.chem.hyperspace3d.feature;

import java.util.List;
import java.util.Set;

public final class DeepSpaceFeatureSchema {
    public static final int MAX_ATOMS = 32;
    public static final int BASE_ATOM_FEATURE_DIM = 50;
    public static final int ATOM_FEATURE_DIM = 56;
    public static final int PAIR_FEATURE_DIM = 36;
    public static final int EMBEDDING_DIM = 128;
    public static final String GRAPH_FEATURE_CONTRACT = "deepspace7-v3";
    public static final List<Integer> ELEMENTS =
            List.of(5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53);
    public static final Set<Integer> ELEMENT_SET = Set.copyOf(ELEMENTS);

    private DeepSpaceFeatureSchema() {}
}
