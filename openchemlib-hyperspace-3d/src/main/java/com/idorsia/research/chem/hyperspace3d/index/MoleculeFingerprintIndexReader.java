package com.idorsia.research.chem.hyperspace3d.index;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.Closeable;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Memory-mapped reader for one flat molecule fingerprint shard. */
public final class MoleculeFingerprintIndexReader implements Closeable {
    private final MoleculeFingerprintIndexShard shard;
    private final int shardIndex;
    public enum VectorColumns { BOTH, BASE_128, COMPACT_16 }
    private final FileChannel baseChannel;
    private final FileChannel compactChannel;
    private final FileChannel rowsChannel;
    private final ByteBuffer base;
    private final ByteBuffer compact;
    private final ByteBuffer rows;
    private final RandomAccessFile strings;
    private long rowIndex;

    public static MoleculeFingerprintIndexManifest loadManifest(Path directory)
            throws IOException {
        MoleculeFingerprintIndexManifest manifest = new ObjectMapper().readValue(
                directory.resolve("manifest.json").toFile(),
                MoleculeFingerprintIndexManifest.class);
        manifest.validate();
        for (MoleculeFingerprintIndexShard shard : manifest.shards) {
            Path path = directory.resolve(shard.directory);
            if (!Files.isRegularFile(path.resolve(".complete"))) {
                throw new IOException("incomplete molecule-index shard: " + shard.directory);
            }
        }
        return manifest;
    }

    public MoleculeFingerprintIndexReader(Path directory,
            MoleculeFingerprintIndexManifest manifest, int shardIndex) throws IOException {
        this(directory, manifest, shardIndex, VectorColumns.BOTH);
    }

    public MoleculeFingerprintIndexReader(Path directory,
            MoleculeFingerprintIndexManifest manifest, int shardIndex,
            VectorColumns columns) throws IOException {
        if (shardIndex < 0 || shardIndex >= manifest.shards.size()) {
            throw new IllegalArgumentException("invalid shard index");
        }
        if (columns == null) throw new IllegalArgumentException("vector columns are required");
        this.shardIndex = shardIndex;
        shard = manifest.shards.get(shardIndex);
        Path shardDirectory = directory.resolve(shard.directory);
        baseChannel = columns == VectorColumns.COMPACT_16 ? null
                : FileChannel.open(shardDirectory.resolve("vectors-128.fp16"));
        compactChannel = columns == VectorColumns.BASE_128 ? null
                : FileChannel.open(shardDirectory.resolve("vectors-16.fp16"));
        rowsChannel = FileChannel.open(shardDirectory.resolve("rows.bin"));
        base = baseChannel == null ? null : map(baseChannel);
        compact = compactChannel == null ? null : map(compactChannel);
        rows = map(rowsChannel);
        strings = new RandomAccessFile(shardDirectory.resolve("strings.bin").toFile(), "r");
        if (base != null) {
            checkHeader(base, MoleculeFingerprintIndexWriter.BASE_MAGIC,
                    128, 1, shard.recordCount);
            requireSize(baseChannel, MoleculeFingerprintIndexWriter.HEADER_BYTES
                    + Math.multiplyExact(shard.recordCount, 256L), "base vector");
        }
        if (compact != null) {
            checkHeader(compact, MoleculeFingerprintIndexWriter.COMPACT_MAGIC,
                    16, 1, shard.recordCount);
            requireSize(compactChannel, MoleculeFingerprintIndexWriter.HEADER_BYTES
                    + Math.multiplyExact(shard.recordCount, 32L), "compact vector");
        }
        checkHeader(rows, MoleculeFingerprintIndexWriter.ROW_MAGIC,
                MoleculeFingerprintIndexWriter.ROW_BYTES, 0, shard.recordCount);
        checkHeader(strings, MoleculeFingerprintIndexWriter.STRING_MAGIC, 0, 2,
                shard.recordCount);
        requireSize(rowsChannel, MoleculeFingerprintIndexWriter.HEADER_BYTES
                + Math.multiplyExact(shard.recordCount,
                        MoleculeFingerprintIndexWriter.ROW_BYTES), "row");
    }

    public List<MoleculeFingerprintRecord> readBatch(int maximum) throws IOException {
        if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
        List<MoleculeFingerprintRecord> result = new ArrayList<>(
                (int) Math.min(maximum, shard.recordCount - rowIndex));
        while (result.size() < maximum && rowIndex < shard.recordCount) {
            result.add(read(rowIndex++));
        }
        return result;
    }

    /** Reads only compact vectors and stable row references; molecule strings stay untouched. */
    public MoleculeCompactVectorBatch readCompactBatch(int maximum) {
        if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
        if (compact == null) {
            throw new IllegalStateException("compact vector column was not opened");
        }
        int count = (int) Math.min(maximum, shard.recordCount - rowIndex);
        float[][] vectors = new float[count][];
        List<MoleculeVectorReference> references = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            long localRow = rowIndex++;
            vectors[i] = vector(compact, localRow, 16);
            references.add(new MoleculeVectorReference(shardIndex, localRow));
        }
        return new MoleculeCompactVectorBatch(vectors, references);
    }

    /** Reads one selected vector column into a flat row-major batch. */
    public MoleculeVectorBatch readVectorBatch(int maximum, int dimension) {
        if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
        ByteBuffer data = dimension == 128 ? base : dimension == 16 ? compact : null;
        if (data == null) throw new IllegalStateException("requested vector column was not opened");
        int count = (int) Math.min(maximum, shard.recordCount - rowIndex);
        float[] values = new float[Math.multiplyExact(count, dimension)];
        List<MoleculeVectorReference> references = new ArrayList<>(count);
        for (int row = 0; row < count; row++) {
            long localRow = rowIndex++;
            copyVector(data, localRow, dimension, values, row * dimension);
            references.add(new MoleculeVectorReference(shardIndex, localRow));
        }
        return new MoleculeVectorBatch(values, dimension, references);
    }

    public MoleculeFingerprintRecord read(long localRow) throws IOException {
        if (base == null || compact == null) {
            throw new IllegalStateException("resolved records require both vector columns");
        }
        if (localRow < 0 || localRow >= shard.recordCount) {
            throw new IndexOutOfBoundsException("local row outside shard");
        }
        MoleculeFingerprintMetadata metadata = readMetadata(localRow);
        return new MoleculeFingerprintRecord(metadata.sourceRow(), metadata.heavyAtomCount(),
                metadata.moleculeId(), metadata.smiles(),
                vector(base, localRow, 128), vector(compact, localRow, 16));
    }

    /** Resolves row metadata without loading either fingerprint column. */
    public MoleculeFingerprintMetadata readMetadata(long localRow) throws IOException {
        if (localRow < 0 || localRow >= shard.recordCount) {
            throw new IndexOutOfBoundsException("local row outside shard");
        }
        int rowOffset = Math.toIntExact(MoleculeFingerprintIndexWriter.HEADER_BYTES
                + localRow * MoleculeFingerprintIndexWriter.ROW_BYTES);
        long sourceRow = rows.getLong(rowOffset);
        int heavyAtoms = rows.getInt(rowOffset + 8);
        long stringOffset = rows.getLong(rowOffset + 16);
        String[] metadata = readStrings(stringOffset);
        return new MoleculeFingerprintMetadata(sourceRow, heavyAtoms, metadata[0], metadata[1]);
    }

    public float[] base128(long localRow) { return vectorChecked(base, localRow, 128); }
    public float[] compact16(long localRow) { return vectorChecked(compact, localRow, 16); }
    public long recordCount() { return shard.recordCount; }

    private float[] vectorChecked(ByteBuffer data, long localRow, int width) {
        if (data == null) {
            throw new IllegalStateException("requested vector column was not opened");
        }
        if (localRow < 0 || localRow >= shard.recordCount) {
            throw new IndexOutOfBoundsException("local row outside shard");
        }
        return vector(data, localRow, width);
    }

    private static float[] vector(ByteBuffer data, long localRow, int width) {
        float[] result = new float[width];
        copyVector(data, localRow, width, result, 0);
        return result;
    }

    private static void copyVector(ByteBuffer data, long localRow, int width,
            float[] output, int outputOffset) {
        int offset = Math.toIntExact(MoleculeFingerprintIndexWriter.HEADER_BYTES
                + localRow * width * 2L);
        for (int i = 0; i < width; i++) {
            output[outputOffset + i] = Float.float16ToFloat(data.getShort(offset + i * 2));
        }
    }

    private String[] readStrings(long offset) throws IOException {
        if (offset < MoleculeFingerprintIndexWriter.HEADER_BYTES
                || offset > strings.length() - 8) throw new IOException("invalid metadata offset");
        strings.seek(offset);
        int idLength = readIntLE(strings);
        int smilesLength = readIntLE(strings);
        if (idLength < 0 || smilesLength < 0
                || (long) idLength + smilesLength > strings.length() - offset - 8) {
            throw new IOException("invalid molecule metadata length");
        }
        byte[] id = new byte[idLength];
        byte[] smiles = new byte[smilesLength];
        strings.readFully(id);
        strings.readFully(smiles);
        return new String[]{new String(id, StandardCharsets.UTF_8),
                new String(smiles, StandardCharsets.UTF_8)};
    }

    @Override public void close() throws IOException {
        IOException failure = null;
        try { close(baseChannel); } catch (IOException error) { failure = error; }
        try { close(compactChannel); } catch (IOException error) {
            if (failure == null) failure = error;
        }
        try { rowsChannel.close(); } catch (IOException error) {
            if (failure == null) failure = error;
        }
        try { strings.close(); } catch (IOException error) {
            if (failure == null) failure = error;
        }
        if (failure != null) throw failure;
    }


    private static void close(FileChannel channel) throws IOException {
        if (channel != null) channel.close();
    }
    private static ByteBuffer map(FileChannel channel) throws IOException {
        return channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size())
                .order(ByteOrder.LITTLE_ENDIAN);
    }

    private static void requireSize(FileChannel channel, long expected, String label)
            throws IOException {
        if (channel.size() != expected) {
            throw new IOException(label + " file size does not match its row count");
        }
    }

    private static void checkHeader(ByteBuffer input, int magic, int width, int dtype,
            long count) throws IOException {
        if (input.capacity() < MoleculeFingerprintIndexWriter.HEADER_BYTES
                || input.getInt(0) != magic
                || input.getInt(4) != MoleculeFingerprintIndexWriter.FORMAT_VERSION
                || input.getInt(8) != width || input.getInt(12) != dtype
                || input.getLong(16) != count) {
            throw new IOException("incompatible molecule-index component");
        }
    }

    private static void checkHeader(RandomAccessFile input, int magic, int width,
            int dtype, long count) throws IOException {
        if (input.length() < MoleculeFingerprintIndexWriter.HEADER_BYTES) {
            throw new IOException("truncated molecule-index component");
        }
        input.seek(0);
        if (readIntLE(input) != magic
                || readIntLE(input) != MoleculeFingerprintIndexWriter.FORMAT_VERSION
                || readIntLE(input) != width || readIntLE(input) != dtype
                || readLongLE(input) != count) {
            throw new IOException("incompatible molecule-index component");
        }
    }

    private static int readIntLE(java.io.DataInput input) throws IOException {
        return Integer.reverseBytes(input.readInt());
    }
    private static long readLongLE(java.io.DataInput input) throws IOException {
        return Long.reverseBytes(input.readLong());
    }
}
