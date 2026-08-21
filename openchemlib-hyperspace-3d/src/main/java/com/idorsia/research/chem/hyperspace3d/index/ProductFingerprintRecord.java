package com.idorsia.research.chem.hyperspace3d.index;

import java.util.Arrays;

public record ProductFingerprintRecord(ProductTuple tuple, float[] embedding) {
    public ProductFingerprintRecord {
        if (embedding.length != 128) throw new IllegalArgumentException("embedding must be 128D");
        embedding = Arrays.copyOf(embedding, embedding.length);
    }
}
