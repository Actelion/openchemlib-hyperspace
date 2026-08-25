package com.idorsia.research.chem.hyperspace3d.index;

import java.io.IOException;

/** Sequential and random vector access for one molecule-fingerprint shard. */
public interface MoleculeFingerprintShardReader extends AutoCloseable {
    MoleculeVectorBatch readBatch(int maximum) throws IOException;
    float[] readVector(long localRow) throws IOException;
    long recordCount();
    @Override void close() throws IOException;
}
