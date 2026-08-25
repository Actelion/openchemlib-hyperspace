package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;

/** One learned flat-library hit, optionally resolved and exact-reranked. */
public record MoleculeSkelSpheresHit(MoleculeVectorReference reference,
        MoleculeFingerprintMetadata molecule, double dotProduct,
        double predictedSimilarity, Double exactSimilarity) {
    public MoleculeSkelSpheresHit withMolecule(MoleculeFingerprintMetadata value) {
        return new MoleculeSkelSpheresHit(reference, value, dotProduct,
                predictedSimilarity, exactSimilarity);
    }
    public MoleculeSkelSpheresHit withExact(double value) {
        return new MoleculeSkelSpheresHit(reference, molecule, dotProduct,
                predictedSimilarity, value);
    }
}
