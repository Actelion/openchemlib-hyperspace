package com.idorsia.research.chem.hyperspace3d.benchmark;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSources;
import java.nio.file.Path;
import java.util.Locale;

/** Sequential FP16 read/decode benchmark for native or Python molecule caches. */
public final class MoleculeCacheReadBenchmark {
    private MoleculeCacheReadBenchmark() {}

    public static void main(String[] args) throws Exception {
        if (args.length < 1 || args.length > 4) {
            throw new IllegalArgumentException(
                    "usage: MoleculeCacheReadBenchmark CACHE [BASE_128|COMPACT_16] [SHARDS] [BATCH]");
        }
        MoleculeFingerprintColumn column = args.length >= 2
                ? MoleculeFingerprintColumn.valueOf(args[1].toUpperCase(Locale.ROOT))
                : MoleculeFingerprintColumn.BASE_128;
        int requestedShards = args.length >= 3 ? Integer.parseInt(args[2]) : 1;
        int batchSize = args.length >= 4 ? Integer.parseInt(args[3]) : 32768;
        long records = 0;
        double checksum = 0;
        try (var source = MoleculeFingerprintDataSources.open(Path.of(args[0]))) {
            int shards = Math.min(requestedShards, source.shardCount());
            long started = System.nanoTime();
            for (int shard = 0; shard < shards; shard++) {
                try (var reader = source.openShard(shard, column)) {
                    while (true) {
                        var batch = reader.readBatch(batchSize);
                        if (batch.size() == 0) break;
                        for (int row = 0; row < batch.size(); row++) {
                            checksum += batch.values()[row * column.dimension()];
                        }
                        records += batch.size();
                    }
                }
            }
            double seconds = (System.nanoTime() - started) / 1.0e9;
            double gib = records * column.dimension() * 2.0 / (1024 * 1024 * 1024);
            System.out.printf(Locale.ROOT, "artifact=%s column=%s shards=%d records=%d%n",
                    source.artifactType(), column, shards, records);
            System.out.printf(Locale.ROOT,
                    "elapsed_seconds=%.6f molecules_per_second=%.0f fp16_gib_per_second=%.3f checksum=%.6f%n",
                    seconds, records / seconds, gib / seconds, checksum);
        }
    }
}
