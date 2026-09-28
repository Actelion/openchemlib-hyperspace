package com.idorsia.research.chem.hyperspace3d.index;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Sequential vector scan with lazy random tuple resolution. */
public final class ProductVectorIndexReader implements Closeable {
    private final ProductVectorIndexManifest manifest;
    private final int shardIndex;
    private final DataInputStream vectors;
    private final DataInputStream rows;
    private final RandomAccessFile tuples;
    private long rowIndex;

    public static ProductVectorIndexManifest loadManifest(Path directory, boolean verifyHashes) throws IOException {
        ProductVectorIndexManifest manifest = new ObjectMapper().readValue(
                directory.resolve("manifest.json").toFile(), ProductVectorIndexManifest.class);
        manifest.validate();
        if (verifyHashes) for (ProductVectorIndexShard shard : manifest.shards) {
            verify(directory.resolve(shard.vectorsPath()), shard.vectorsSha256());
            verify(directory.resolve(shard.rowsPath()), shard.rowsSha256());
            verify(directory.resolve(shard.tuplesPath()), shard.tuplesSha256());
        }
        return manifest;
    }

    public ProductVectorIndexReader(Path directory, ProductVectorIndexManifest manifest, int shardIndex)
            throws IOException {
        this.manifest = manifest; this.shardIndex = shardIndex;
        ProductVectorIndexShard shard = manifest.shards.get(shardIndex);
        vectors = input(directory.resolve(shard.vectorsPath())); rows = input(directory.resolve(shard.rowsPath()));
        tuples = new RandomAccessFile(directory.resolve(shard.tuplesPath()).toFile(), "r");
        checkHeader(vectors, ProductVectorIndexWriter.VECTOR_MAGIC, 16, 1);
        checkHeader(rows, ProductVectorIndexWriter.ROW_MAGIC, 12, 0);
        checkHeader(tuples, ProductVectorIndexWriter.TUPLE_MAGIC, 0, 0);
    }

    public ProductVectorBatch readBatch(int maximum) throws IOException {
        if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
        List<ProductVectorReference> refs = new ArrayList<>(maximum);
        List<float[]> values = new ArrayList<>(maximum);
        while (refs.size() < maximum) {
            try {
                int reaction = readIntLE(rows); long offset = readLongLE(rows);
                if (reaction < 0 || reaction >= manifest.reactionDictionary.size()) throw new IOException("invalid reaction ordinal");
                float[] vector = new float[16];
                for (int i = 0; i < 16; i++) vector[i] = Float.float16ToFloat(readShortLE(vectors));
                refs.add(new ProductVectorReference(shardIndex, rowIndex++, reaction, offset)); values.add(vector);
            } catch (EOFException end) { break; }
        }
        return new ProductVectorBatch(refs, values.toArray(float[][]::new));
    }

    public ProductTuple resolveTuple(ProductVectorReference ref) throws IOException {
        if (ref.shardIndex() != shardIndex) throw new IllegalArgumentException("reference belongs to another shard");
        tuples.seek(ref.tupleOffset());
        int count = readIntLE(tuples);
        if (count < 1 || count > 16) throw new IOException("invalid synthon count");
        List<String> ids = new ArrayList<>(count); List<Integer> ordinals = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            int length = readIntLE(tuples);
            if (length < 0 || length > 1_000_000) throw new IOException("invalid synthon id length");
            byte[] bytes = new byte[length]; tuples.readFully(bytes);
            ids.add(new String(bytes, StandardCharsets.UTF_8)); ordinals.add(readIntLE(tuples));
        }
        return new ProductTuple(manifest.reactionDictionary.get(ref.reactionOrdinal()), ids, ordinals);
    }

    @Override public void close() throws IOException { vectors.close(); rows.close(); tuples.close(); }
    private static DataInputStream input(Path path) throws IOException { return new DataInputStream(new BufferedInputStream(Files.newInputStream(path))); }
    private static void verify(Path path, String expected) throws IOException {
        if (!Files.isRegularFile(path) || !ProductFingerprintIndexBuilder.sha256(path).equals(expected))
            throw new IOException("index component checksum mismatch: " + path.getFileName());
    }
    private static void checkHeader(DataInputStream in, int magic, int width, int dtype) throws IOException {
        if (readIntLE(in) != magic || readIntLE(in) != 2 || readIntLE(in) != width || readIntLE(in) != dtype) throw new IOException("incompatible index component");
        readLongLE(in);
    }
    private static void checkHeader(RandomAccessFile in, int magic, int width, int dtype) throws IOException {
        if (readIntLE(in) != magic || readIntLE(in) != 2 || readIntLE(in) != width || readIntLE(in) != dtype) throw new IOException("incompatible tuple component");
        readLongLE(in);
    }
    private static int readIntLE(java.io.DataInput in) throws IOException { return Integer.reverseBytes(in.readInt()); }
    private static long readLongLE(java.io.DataInput in) throws IOException { return Long.reverseBytes(in.readLong()); }
    private static short readShortLE(java.io.DataInput in) throws IOException { return Short.reverseBytes(in.readShort()); }
}
