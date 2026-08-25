package com.idorsia.research.chem.hyperspace3d.index;

import java.io.IOException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/** Bounded background FP16 decoding pipeline for sequential multi-shard scans. */
public final class MoleculeVectorBatchPrefetcher implements AutoCloseable {
    private final BlockingQueue<Item> queue;
    private final Thread producer;
    private volatile boolean closed;

    public MoleculeVectorBatchPrefetcher(MoleculeFingerprintDataSource source,
            MoleculeFingerprintColumn column, int batchSize, int prefetchDepth) {
        if (batchSize < 1 || prefetchDepth < 1) {
            throw new IllegalArgumentException("batch size and prefetch depth must be positive");
        }
        queue = new ArrayBlockingQueue<>(prefetchDepth);
        producer = Thread.ofPlatform().name("hyperspace3d-vector-prefetch").daemon(true)
                .start(() -> produce(source, column, batchSize));
    }

    /** Returns null after the last shard. */
    public MoleculeVectorBatch take() throws IOException {
        try {
            Item item = queue.take();
            if (item.error != null) throw item.error;
            return item.completed ? null : item.batch;
        } catch (InterruptedException error) {
            Thread.currentThread().interrupt();
            throw new IOException("interrupted while waiting for fingerprint vectors", error);
        }
    }

    private void produce(MoleculeFingerprintDataSource source,
            MoleculeFingerprintColumn column, int batchSize) {
        try {
            for (int shard = 0; shard < source.shardCount() && !closed; shard++) {
                try (var reader = source.openShard(shard, column)) {
                    while (!closed) {
                        MoleculeVectorBatch batch = reader.readBatch(batchSize);
                        if (batch.size() == 0) break;
                        queue.put(Item.batch(batch));
                    }
                }
            }
            if (!closed) queue.put(Item.endMarker());
        } catch (InterruptedException interrupted) {
            Thread.currentThread().interrupt();
        } catch (IOException | RuntimeException error) {
            if (!closed) {
                try { queue.put(Item.error(error instanceof IOException io
                        ? io : new IOException("fingerprint prefetch failed", error))); }
                catch (InterruptedException interrupted) { Thread.currentThread().interrupt(); }
            }
        }
    }

    @Override public void close() {
        closed = true;
        producer.interrupt();
        try { producer.join(); }
        catch (InterruptedException error) { Thread.currentThread().interrupt(); }
    }

    private record Item(MoleculeVectorBatch batch, IOException error, boolean completed) {
        static Item batch(MoleculeVectorBatch value) { return new Item(value, null, false); }
        static Item error(IOException value) { return new Item(null, value, false); }
        static Item endMarker() { return new Item(null, null, true); }
    }
}
