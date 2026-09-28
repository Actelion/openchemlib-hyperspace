package com.idorsia.research.chem.hyperspace3d.index;

import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public final class ProductFingerprintIndexReader implements Closeable {
    private final DataInputStream input;

    public ProductFingerprintIndexReader(Path path) throws IOException {
        input = new DataInputStream(new BufferedInputStream(java.nio.file.Files.newInputStream(path)));
        if (input.readInt() != ProductFingerprintIndexWriter.MAGIC
                || input.readInt() != 1 || input.readInt() != 128) {
            throw new IOException("incompatible product-index shard");
        }
    }

    public List<ProductFingerprintRecord> readBatch(int maximum) throws IOException {
        if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
        List<ProductFingerprintRecord> records = new ArrayList<>(maximum);
        while (records.size() < maximum) {
            try {
                String reaction = input.readUTF();
                int count = input.readInt();
                if (count < 1 || count > 16) throw new IOException("invalid synthon count");
                List<String> ids = new ArrayList<>(count);
                List<Integer> ordinals = new ArrayList<>(count);
                for (int i = 0; i < count; i++) {
                    ids.add(input.readUTF());
                    ordinals.add(input.readInt());
                }
                float[] embedding = new float[128];
                for (int i = 0; i < embedding.length; i++) embedding[i] = input.readFloat();
                records.add(new ProductFingerprintRecord(
                        new ProductTuple(reaction, ids, ordinals), embedding));
            } catch (EOFException end) {
                break;
            }
        }
        return records;
    }

    @Override public void close() throws IOException { input.close(); }
}
