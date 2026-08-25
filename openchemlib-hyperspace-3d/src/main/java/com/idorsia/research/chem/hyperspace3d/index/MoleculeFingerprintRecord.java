package com.idorsia.research.chem.hyperspace3d.index;

/** One resolved row from a flat molecule fingerprint index. */
public record MoleculeFingerprintRecord(
        long sourceRow, int heavyAtomCount, String moleculeId, String smiles,
        float[] base128, float[] compact16) {
}
