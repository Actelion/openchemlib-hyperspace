package com.idorsia.research.chem.hyperspace3d.index;

/** Resolved non-vector fields for one flat molecule-index row. */
public record MoleculeFingerprintMetadata(long sourceRow, int heavyAtomCount,
        String moleculeId, String smiles) {}
