package com.idorsia.research.chem.hyperspace3d.model;

/** Universal and compact fingerprints produced from the same molecule batch. */
public record MoleculeFingerprintBatch(float[][] base128, float[][] compact16) {
    public MoleculeFingerprintBatch {
        if (base128 == null || compact16 == null || base128.length != compact16.length) {
            throw new IllegalArgumentException("fingerprint batch row counts must match");
        }
        for (int row = 0; row < base128.length; row++) {
            if (base128[row] == null || base128[row].length != 128) {
                throw new IllegalArgumentException("base fingerprint must be 128D");
            }
            if (compact16[row] == null || compact16[row].length != 16) {
                throw new IllegalArgumentException("compact fingerprint must be 16D");
            }
        }
    }

    public int size() { return base128.length; }
}
