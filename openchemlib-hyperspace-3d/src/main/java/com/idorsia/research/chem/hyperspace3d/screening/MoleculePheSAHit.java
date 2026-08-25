package com.idorsia.research.chem.hyperspace3d.screening;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;

/** One learned 128D similarity hit with optional compact-screen provenance. */
public record MoleculePheSAHit(MoleculeVectorReference reference,
        MoleculeFingerprintMetadata molecule, double objectiveScore, float[] modelScores,
        Double predictedSkelSpheres, Double skelSpheresDotProduct, Integer compactRank,
        Double exactSkelSpheres) {
    public MoleculePheSAHit withMolecule(MoleculeFingerprintMetadata value) {
        return new MoleculePheSAHit(reference, value, objectiveScore, modelScores,
                predictedSkelSpheres, skelSpheresDotProduct, compactRank, exactSkelSpheres);
    }
    public MoleculePheSAHit withExactSkelSpheres(Double value) {
        return new MoleculePheSAHit(reference, molecule, objectiveScore, modelScores,
                predictedSkelSpheres, skelSpheresDotProduct, compactRank, value);
    }
}
