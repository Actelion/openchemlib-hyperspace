package com.idorsia.research.chem.hyperspace3d.mining;

/** One immutable learned/exact observation in a query shard. */
public record PheSAQueryPairRecord(
        long candidateId,
        float predictedTotal,
        float predictedShape,
        float predictedPharmacophore,
        long rankTotal,
        long rankShape,
        long rankPharmacophore,
        float exactTotal,
        float exactShape,
        float exactPharmacophore,
        long selectionMask,
        int miningRound,
        int exactStatus) {

    public static final int EXACT_OK = 0;
    public static final int EXACT_DESCRIPTOR_FAILED = 1;
    public static final int EXACT_ALIGNMENT_FAILED = 2;
    public static final int EXACT_NONFINITE = 3;
    public static final int EXACT_COMPONENT_MISMATCH = 4;

    public boolean exactValid() {
        return exactStatus == EXACT_OK && Float.isFinite(exactTotal)
                && Float.isFinite(exactShape) && Float.isFinite(exactPharmacophore);
    }

    public PheSAQueryPairRecord withExact(ExactPheSALabel value) {
        return new PheSAQueryPairRecord(candidateId, predictedTotal, predictedShape,
                predictedPharmacophore, rankTotal, rankShape, rankPharmacophore,
                value.total(), value.shape(), value.pharmacophore(), selectionMask,
                miningRound, value.status());
    }
}
