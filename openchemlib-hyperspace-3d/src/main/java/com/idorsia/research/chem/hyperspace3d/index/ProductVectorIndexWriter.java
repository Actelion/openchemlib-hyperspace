package com.idorsia.research.chem.hyperspace3d.index;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Writes independently streamable FP16 vectors, fixed rows, and variable tuples. */
public final class ProductVectorIndexWriter implements Closeable {
    static final int VECTOR_MAGIC = 0x48565032;
    static final int ROW_MAGIC = 0x48525032;
    static final int TUPLE_MAGIC = 0x48545032;
    static final int HEADER_BYTES = 24;
    private final DataOutputStream vectors;
    private final DataOutputStream rows;
    private final DataOutputStream tuples;
    private final Map<String, Integer> reactionOrdinals;
    private long tupleOffset = HEADER_BYTES;
    private long count;

    public ProductVectorIndexWriter(Path vectorsPath, Path rowsPath, Path tuplesPath,
            List<String> reactionDictionary) throws IOException {
        reactionOrdinals = new HashMap<>();
        for (int i = 0; i < reactionDictionary.size(); i++) reactionOrdinals.put(reactionDictionary.get(i), i);
        vectors = output(vectorsPath); rows = output(rowsPath); tuples = output(tuplesPath);
        header(vectors, VECTOR_MAGIC, 16, 1);
        header(rows, ROW_MAGIC, 12, 0);
        header(tuples, TUPLE_MAGIC, 0, 0);
    }

    public void write(ProductTuple tuple, float[] vector) throws IOException {
        if (vector.length != 16) throw new IllegalArgumentException("compact vector must have 16 values");
        Integer reaction = reactionOrdinals.get(tuple.reactionId());
        if (reaction == null) throw new IllegalArgumentException("reaction is absent from dictionary: " + tuple.reactionId());
        writeIntLE(rows, reaction); writeLongLE(rows, tupleOffset);
        for (float value : vector) writeShortLE(vectors, Float.floatToFloat16(value));
        byte[] encoded = encodeTuple(tuple);
        tuples.write(encoded); tupleOffset += encoded.length; count++;
    }

    public long recordCount() { return count; }
    public void flush() throws IOException { vectors.flush(); rows.flush(); tuples.flush(); }
    @Override public void close() throws IOException {
        IOException failure = null;
        try { vectors.close(); } catch (IOException e) { failure = e; }
        try { rows.close(); } catch (IOException e) { if (failure == null) failure = e; }
        try { tuples.close(); } catch (IOException e) { if (failure == null) failure = e; }
        if (failure != null) throw failure;
    }

    private static byte[] encodeTuple(ProductTuple tuple) throws IOException {
        var bytes = new java.io.ByteArrayOutputStream();
        var out = new DataOutputStream(bytes);
        writeIntLE(out, tuple.synthonIds().size());
        for (int i = 0; i < tuple.synthonIds().size(); i++) {
            byte[] id = tuple.synthonIds().get(i).getBytes(StandardCharsets.UTF_8);
            writeIntLE(out, id.length); out.write(id); writeIntLE(out, tuple.synthonOrdinals().get(i));
        }
        out.flush(); return bytes.toByteArray();
    }

    private static DataOutputStream output(Path path) throws IOException {
        return new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(path)));
    }
    private static void header(DataOutputStream out, int magic, int width, int dtype) throws IOException {
        writeIntLE(out, magic); writeIntLE(out, 2); writeIntLE(out, width); writeIntLE(out, dtype);
        writeLongLE(out, -1);
    }
    static void writeIntLE(DataOutputStream out, int value) throws IOException { out.writeInt(Integer.reverseBytes(value)); }
    static void writeLongLE(DataOutputStream out, long value) throws IOException { out.writeLong(Long.reverseBytes(value)); }
    static void writeShortLE(DataOutputStream out, short value) throws IOException { out.writeShort(Short.reverseBytes(value)); }
}
