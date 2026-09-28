package com.idorsia.research.chem.hyperspace3d.archive;

import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;

public record CandidateArchiveRecord(int formatVersion, ProductTuple tuple,
        double predictedTotalScore, double predictedShapeScore,
        double predictedPharmacophoreScore, String source, String basinId,
        int beamRound, String modelBundleHash, String queryIdentifier,
        String querySha256) {
    public CandidateArchiveRecord {
        if (formatVersion != 1) throw new IllegalArgumentException("unsupported archive version");
    }
}
