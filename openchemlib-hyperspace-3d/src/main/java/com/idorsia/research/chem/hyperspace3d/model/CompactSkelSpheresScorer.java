package com.idorsia.research.chem.hyperspace3d.model;

public final class CompactSkelSpheresScorer {
    private final CompactSkelSpheresManifest manifest;

    public CompactSkelSpheresScorer(CompactSkelSpheresManifest manifest) {
        manifest.validate();
        this.manifest = manifest;
    }

    public Score score(float[] query, float[] candidate) {
        if (query.length != 16 || candidate.length != 16) {
            throw new IllegalArgumentException("compact embeddings must be 16D");
        }
        double dot = 0.0;
        for (int i = 0; i < 16; i++) dot += (double) query[i] * candidate[i];
        return new Score(dot, manifest.calibrate(dot));
    }

    public double calibrated(double dotProduct) {
        return manifest.calibrate(dotProduct);
    }

    public record Score(double dotProduct, double calibratedSimilarity) {}
}
