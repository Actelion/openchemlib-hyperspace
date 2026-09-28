package com.idorsia.research.chem.hyperspace3d.mining;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;

/** Versioned fixed-width query-shard reader/writer shared with Python. */
public final class PheSAQueryShardIO {
    public static final byte[] MAGIC = {'H', '3', 'D', 'Q', 'P', 'R', '0', '1'};
    public static final int FORMAT_VERSION = 1;
    public static final int HEADER_BYTES = 64;
    public static final int RECORD_BYTES = 72;
    private static final ByteOrder ORDER = ByteOrder.LITTLE_ENDIAN;

    private PheSAQueryShardIO() {}

    public static void write(Path path, long queryId, long screenedCount,
            int miningRound, List<PheSAQueryPairRecord> records) throws IOException {
        if (queryId < 0 || screenedCount < records.size() || miningRound < 0) {
            throw new IllegalArgumentException("invalid query shard metadata");
        }
        Path parent = path.toAbsolutePath().normalize().getParent();
        if (parent != null) Files.createDirectories(parent);
        try (FileChannel out = FileChannel.open(path, StandardOpenOption.CREATE_NEW,
                StandardOpenOption.WRITE)) {
            ByteBuffer header = ByteBuffer.allocate(HEADER_BYTES).order(ORDER);
            header.put(MAGIC).putInt(FORMAT_VERSION).putInt(RECORD_BYTES)
                    .putLong(records.size()).putLong(screenedCount).putLong(queryId)
                    .putInt(miningRound).putInt(0);
            header.position(HEADER_BYTES).flip();
            writeFully(out, header);
            ByteBuffer buffer = ByteBuffer.allocate(RECORD_BYTES * Math.min(4096,
                    Math.max(1, records.size()))).order(ORDER);
            for (PheSAQueryPairRecord value : records) {
                if (buffer.remaining() < RECORD_BYTES) {
                    buffer.flip(); writeFully(out, buffer); buffer.clear();
                }
                put(buffer, value);
            }
            buffer.flip(); writeFully(out, buffer);
        }
    }

    public static Shard read(Path path) throws IOException {
        try (FileChannel input = FileChannel.open(path)) {
            ByteBuffer header = ByteBuffer.allocate(HEADER_BYTES).order(ORDER);
            readFully(input, header); header.flip();
            byte[] magic = new byte[MAGIC.length]; header.get(magic);
            if (!java.util.Arrays.equals(magic, MAGIC)
                    || header.getInt() != FORMAT_VERSION || header.getInt() != RECORD_BYTES) {
                throw new IOException("unsupported PheSA query shard: " + path);
            }
            long count = header.getLong();
            long screened = header.getLong();
            long query = header.getLong();
            int round = header.getInt();
            header.getInt();
            if (count < 0 || count > Integer.MAX_VALUE || screened < count
                    || query < 0 || round < 0
                    || input.size() != HEADER_BYTES + count * RECORD_BYTES) {
                throw new IOException("inconsistent PheSA query shard metadata");
            }
            List<PheSAQueryPairRecord> values = new ArrayList<>((int) count);
            ByteBuffer rows = ByteBuffer.allocate(RECORD_BYTES * Math.min(4096,
                    Math.max(1, (int) count))).order(ORDER);
            while (values.size() < count) {
                rows.clear();
                rows.limit(RECORD_BYTES * Math.min(4096, (int) count - values.size()));
                readFully(input, rows); rows.flip();
                while (rows.remaining() >= RECORD_BYTES) values.add(get(rows));
            }
            return new Shard(query, screened, round, List.copyOf(values));
        }
    }

    private static void put(ByteBuffer out, PheSAQueryPairRecord value) {
        out.putLong(value.candidateId());
        out.putFloat(value.predictedTotal()).putFloat(value.predictedShape())
                .putFloat(value.predictedPharmacophore());
        out.putLong(value.rankTotal()).putLong(value.rankShape())
                .putLong(value.rankPharmacophore());
        out.putFloat(value.exactTotal()).putFloat(value.exactShape())
                .putFloat(value.exactPharmacophore());
        out.putLong(value.selectionMask()).putInt(value.miningRound())
                .putInt(value.exactStatus());
    }

    private static PheSAQueryPairRecord get(ByteBuffer input) {
        return new PheSAQueryPairRecord(input.getLong(), input.getFloat(), input.getFloat(),
                input.getFloat(), input.getLong(), input.getLong(), input.getLong(),
                input.getFloat(), input.getFloat(), input.getFloat(), input.getLong(),
                input.getInt(), input.getInt());
    }

    private static void writeFully(FileChannel channel, ByteBuffer value) throws IOException {
        while (value.hasRemaining()) channel.write(value);
    }
    private static void readFully(FileChannel channel, ByteBuffer value) throws IOException {
        while (value.hasRemaining()) {
            if (channel.read(value) < 0) throw new IOException("truncated PheSA query shard");
        }
    }

    public record Shard(long queryId, long screenedCount, int miningRound,
            List<PheSAQueryPairRecord> records) {}
}
