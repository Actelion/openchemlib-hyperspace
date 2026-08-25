package com.idorsia.research.chem.hyperspace3d.index;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/** Writes independently memory-mappable 128D/16D FP16 columns and molecule metadata. */
public final class MoleculeFingerprintIndexWriter implements Closeable {
    static final int BASE_MAGIC = 0x484d4231;
    static final int COMPACT_MAGIC = 0x484d4331;
    static final int ROW_MAGIC = 0x484d5231;
    static final int STRING_MAGIC = 0x484d5331;
    static final int FORMAT_VERSION = 1;
    static final int HEADER_BYTES = 24;
    static final int ROW_BYTES = 24;
    private final Path basePath;
    private final Path compactPath;
    private final Path rowsPath;
    private final Path stringsPath;
    private final DataOutputStream base;
    private final DataOutputStream compact;
    private final DataOutputStream rows;
    private final DataOutputStream strings;
    private long stringOffset = HEADER_BYTES;
    private long count;
    private boolean closed;

    public MoleculeFingerprintIndexWriter(Path shardDirectory) throws IOException {
        Files.createDirectories(shardDirectory);
        basePath = shardDirectory.resolve("vectors-128.fp16");
        compactPath = shardDirectory.resolve("vectors-16.fp16");
        rowsPath = shardDirectory.resolve("rows.bin");
        stringsPath = shardDirectory.resolve("strings.bin");
        base = output(basePath);
        compact = output(compactPath);
        rows = output(rowsPath);
        strings = output(stringsPath);
        header(base, BASE_MAGIC, 128, 1);
        header(compact, COMPACT_MAGIC, 16, 1);
        header(rows, ROW_MAGIC, ROW_BYTES, 0);
        header(strings, STRING_MAGIC, 0, 2);
    }

    public void write(long sourceRow, int heavyAtomCount, String moleculeId, String smiles,
            float[] base128, float[] compact16) throws IOException {
        requireVector(base128, 128, "base");
        requireVector(compact16, 16, "compact");
        if (sourceRow < 0 || heavyAtomCount < 1 || moleculeId == null || smiles == null) {
            throw new IllegalArgumentException("invalid molecule metadata");
        }
        byte[] id = moleculeId.getBytes(StandardCharsets.UTF_8);
        byte[] structure = smiles.getBytes(StandardCharsets.UTF_8);
        writeLongLE(rows, sourceRow);
        writeIntLE(rows, heavyAtomCount);
        writeIntLE(rows, 0);
        writeLongLE(rows, stringOffset);
        for (float value : base128) writeShortLE(base, Float.floatToFloat16(value));
        for (float value : compact16) writeShortLE(compact, Float.floatToFloat16(value));
        writeIntLE(strings, id.length);
        writeIntLE(strings, structure.length);
        strings.write(id);
        strings.write(structure);
        stringOffset += 8L + id.length + structure.length;
        count++;
    }

    public long recordCount() { return count; }

    @Override public void close() throws IOException {
        if (closed) return;
        closed = true;
        IOException failure = null;
        for (DataOutputStream stream : new DataOutputStream[]{base, compact, rows, strings}) {
            try { stream.close(); } catch (IOException error) {
                if (failure == null) failure = error;
            }
        }
        if (failure == null) {
            patchCount(basePath, count);
            patchCount(compactPath, count);
            patchCount(rowsPath, count);
            patchCount(stringsPath, count);
        }
        if (failure != null) throw failure;
    }

    private static void requireVector(float[] vector, int width, String name) {
        if (vector == null || vector.length != width) {
            throw new IllegalArgumentException(name + " vector must be " + width + "D");
        }
        for (float value : vector) {
            if (!Float.isFinite(value)) throw new IllegalArgumentException(name + " vector is not finite");
        }
    }

    private static void patchCount(Path path, long count) throws IOException {
        try (RandomAccessFile file = new RandomAccessFile(path.toFile(), "rw")) {
            file.seek(16);
            file.writeLong(Long.reverseBytes(count));
        }
    }

    private static DataOutputStream output(Path path) throws IOException {
        return new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(path), 1 << 20));
    }

    private static void header(DataOutputStream out, int magic, int width, int dtype)
            throws IOException {
        writeIntLE(out, magic);
        writeIntLE(out, FORMAT_VERSION);
        writeIntLE(out, width);
        writeIntLE(out, dtype);
        writeLongLE(out, -1);
    }

    static void writeIntLE(DataOutputStream out, int value) throws IOException {
        out.writeInt(Integer.reverseBytes(value));
    }
    static void writeLongLE(DataOutputStream out, long value) throws IOException {
        out.writeLong(Long.reverseBytes(value));
    }
    static void writeShortLE(DataOutputStream out, short value) throws IOException {
        out.writeShort(Short.reverseBytes(value));
    }
}
