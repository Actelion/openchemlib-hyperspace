package com.idorsia.research.chem.hyperspace3d.model;

@FunctionalInterface
public interface EmbeddingComparator {
    float[][] compare(float[] queryEmbedding, float[][] candidateEmbeddings);
}
