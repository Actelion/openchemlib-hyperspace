package com.idorsia.research.chem.hyperspace3d.index;

import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/** Format-neutral access to a sharded flat-molecule fingerprint collection. */
public interface MoleculeFingerprintDataSource extends AutoCloseable {
    String artifactType();
    long recordCount();
    int shardCount();
    MoleculeFingerprintProvenance provenance();
    MoleculeFingerprintShardReader openShard(int shardIndex,
            MoleculeFingerprintColumn column) throws IOException;
    Map<MoleculeVectorReference, MoleculeFingerprintMetadata> resolveMetadata(
            Collection<MoleculeVectorReference> references) throws IOException;
    void validateCompatibility(DeepSpaceModelBundle model,
            CompactSkelSpheresModelBundle compact, boolean compactRequired);
    @Override default void close() throws IOException {}
}
