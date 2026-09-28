package com.idorsia.research.chem.hyperspace3d.model;

public interface EmbeddingComparator {
    float[][] compare(float[] queryEmbedding, float[][] candidateEmbeddings);

    default float[][] compareFlat(float[] queryEmbedding, float[] candidateEmbeddings,
            int candidateCount) {
        if (candidateEmbeddings.length != Math.multiplyExact(candidateCount, 128)) {
            throw new IllegalArgumentException("flat candidate embeddings have invalid length");
        }
        float[][] rows = new float[candidateCount][128];
        for (int row = 0; row < candidateCount; row++) {
            System.arraycopy(candidateEmbeddings, row * 128, rows[row], 0, 128);
        }
        return compare(queryEmbedding, rows);
    }
}
