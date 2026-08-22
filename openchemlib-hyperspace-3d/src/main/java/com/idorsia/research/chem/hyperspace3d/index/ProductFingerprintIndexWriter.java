package com.idorsia.research.chem.hyperspace3d.index;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.DataOutputStream;
import java.io.IOException;
import java.nio.file.Path;

public final class ProductFingerprintIndexWriter implements Closeable {
    static final int MAGIC = 0x48334449;
    private final DataOutputStream output;
    private long recordCount;

    public ProductFingerprintIndexWriter(Path path) throws IOException {
        output = new DataOutputStream(new BufferedOutputStream(java.nio.file.Files.newOutputStream(path)));
        output.writeInt(MAGIC);
        output.writeInt(1);
        output.writeInt(128);
    }

    public void write(ProductFingerprintRecord record) throws IOException {
        ProductTuple tuple = record.tuple();
        output.writeUTF(tuple.reactionId());
        output.writeInt(tuple.synthonIds().size());
        for (int i = 0; i < tuple.synthonIds().size(); i++) {
            output.writeUTF(tuple.synthonIds().get(i));
            output.writeInt(tuple.synthonOrdinals().get(i));
        }
        for (float value : record.embedding()) output.writeFloat(value);
        recordCount++;
    }

    public long recordCount() { return recordCount; }

    public void flush() throws IOException { output.flush(); }

    public static void writeManifest(Path path, ProductFingerprintIndexManifest manifest) throws IOException {
        manifest.validate();
        new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(path.toFile(), manifest);
    }

    @Override public void close() throws IOException { output.close(); }
}
