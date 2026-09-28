package com.idorsia.research.chem.hyperspace3d.mining;

/** Exact OCL PheSA components and a stable status code. */
public record ExactPheSALabel(float total, float shape, float pharmacophore,
        int status, long elapsedNanos) {
    public static ExactPheSALabel failure(int status, long elapsedNanos) {
        return new ExactPheSALabel(Float.NaN, Float.NaN, Float.NaN, status, elapsedNanos);
    }
}
