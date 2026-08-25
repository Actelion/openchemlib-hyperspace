package com.idorsia.research.chem.hyperspace3d.index;

import java.util.List;

/** Compact-only scan batch; metadata strings are resolved only for retained rows. */
public record MoleculeCompactVectorBatch(float[][] vectors,
        List<MoleculeVectorReference> references) {}
