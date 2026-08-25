package com.idorsia.research.chem.hyperspace3d.index;

import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Adapter exposing the native Java dual-vector index through the shared API. */
public final class JavaMoleculeFingerprintDataSource implements MoleculeFingerprintDataSource {
    private final Path directory;
    private final MoleculeFingerprintIndexManifest manifest;

    public JavaMoleculeFingerprintDataSource(Path directory) throws IOException {
        this.directory = directory.toAbsolutePath().normalize();
        this.manifest = MoleculeFingerprintIndexReader.loadManifest(this.directory);
    }

    @Override public String artifactType() { return manifest.artifactType; }
    @Override public long recordCount() { return manifest.recordCount; }
    @Override public int shardCount() { return manifest.shards.size(); }
    @Override public MoleculeFingerprintProvenance provenance() {
        return new MoleculeFingerprintProvenance(null, null, null);
    }

    @Override public MoleculeFingerprintShardReader openShard(int shardIndex,
            MoleculeFingerprintColumn column) throws IOException {
        var selected = column == MoleculeFingerprintColumn.BASE_128
                ? MoleculeFingerprintIndexReader.VectorColumns.BASE_128
                : MoleculeFingerprintIndexReader.VectorColumns.COMPACT_16;
        var reader = new MoleculeFingerprintIndexReader(
                directory, manifest, shardIndex, selected);
        return new MoleculeFingerprintShardReader() {
            @Override public MoleculeVectorBatch readBatch(int maximum) {
                return reader.readVectorBatch(maximum, column.dimension());
            }
            @Override public float[] readVector(long localRow) {
                return column == MoleculeFingerprintColumn.BASE_128
                        ? reader.base128(localRow) : reader.compact16(localRow);
            }
            @Override public long recordCount() { return reader.recordCount(); }
            @Override public void close() throws IOException { reader.close(); }
        };
    }

    @Override public Map<MoleculeVectorReference, MoleculeFingerprintMetadata> resolveMetadata(
            Collection<MoleculeVectorReference> references) throws IOException {
        Map<Integer, List<MoleculeVectorReference>> byShard = new HashMap<>();
        for (var reference : references) {
            requireReference(reference);
            byShard.computeIfAbsent(reference.shardIndex(), ignored -> new ArrayList<>())
                    .add(reference);
        }
        Map<MoleculeVectorReference, MoleculeFingerprintMetadata> result = new HashMap<>();
        for (var entry : byShard.entrySet()) {
            try (var reader = new MoleculeFingerprintIndexReader(directory, manifest,
                    entry.getKey(), MoleculeFingerprintIndexReader.VectorColumns.COMPACT_16)) {
                for (var reference : entry.getValue()) {
                    result.put(reference, reader.readMetadata(reference.localRow()));
                }
            }
        }
        return result;
    }

    @Override public void validateCompatibility(DeepSpaceModelBundle model,
            CompactSkelSpheresModelBundle compact, boolean compactRequired) {
        requirePath(manifest.modelBundle, model.directory(), "primary model bundle");
        if (compactRequired) {
            if (compact == null) throw new IllegalArgumentException("compact model bundle is required");
            requirePath(manifest.compactBundle, compact.directory(), "compact model bundle");
        }
    }

    private void requireReference(MoleculeVectorReference reference) {
        if (reference.shardIndex() < 0 || reference.shardIndex() >= manifest.shards.size()
                || reference.localRow() < 0
                || reference.localRow() >= manifest.shards.get(reference.shardIndex()).recordCount) {
            throw new IllegalArgumentException("molecule vector reference is outside the index");
        }
    }

    private static void requirePath(String indexed, Path selected, String label) {
        if (!Path.of(indexed).toAbsolutePath().normalize().equals(selected.toAbsolutePath().normalize())) {
            throw new IllegalArgumentException(label + " does not match index manifest");
        }
    }
}
