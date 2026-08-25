package com.idorsia.research.chem.hyperspace3d.index;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

/** Direct reader for caches emitted by cache_supplier_fingerprints.py. */
public final class DeepSpace7SupplierCacheDataSource implements MoleculeFingerprintDataSource {
    public static final String ARTIFACT_TYPE = "deepspace7_streamed_supplier_fingerprints_v1";
    private static final String SHARD_ARTIFACT_TYPE =
            "deepspace7_streamed_supplier_fingerprint_shard_v1";
    private static final ObjectMapper MAPPER = new ObjectMapper();
    private final Path directory;
    private final long recordCount;
    private final MoleculeFingerprintProvenance provenance;
    private final List<Shard> shards;

    public DeepSpace7SupplierCacheDataSource(Path directory) throws IOException {
        this.directory = directory.toAbsolutePath().normalize();
        JsonNode root = readJson(required(this.directory.resolve("manifest.json")));
        require(ARTIFACT_TYPE.equals(text(root, "artifact_type")),
                "unsupported Deepspace7 supplier-cache artifact");
        require(root.path("complete_source").asBoolean(false), "supplier cache is incomplete");
        requireFingerprintSummary(root.path("fingerprints"));
        long sourceRows = nonnegative(root, "source_rows");
        this.recordCount = nonnegative(root, "valid_molecules");
        long rejected = nonnegative(root, "rejected_molecules");
        require(Math.addExact(recordCount, rejected) == sourceRows,
                "supplier-cache record totals are inconsistent");
        int shardCount = integer(root, "shards");
        require(shardCount >= 0, "invalid supplier-cache shard count");
        this.provenance = provenance(root.path("provenance"));

        List<Shard> loaded = new ArrayList<>(shardCount);
        long expectedSourceStart = 0;
        long validTotal = 0;
        long rejectedTotal = 0;
        for (int index = 0; index < shardCount; index++) {
            Path shardDirectory = this.directory.resolve(String.format("shard_%06d", index));
            required(shardDirectory.resolve(".complete"));
            JsonNode value = readJson(required(shardDirectory.resolve("manifest.json")));
            require(SHARD_ARTIFACT_TYPE.equals(text(value, "artifact_type")),
                    "unsupported supplier-cache shard artifact at " + index);
            require(integer(value, "shard_index") == index, "non-contiguous shard index");
            long start = nonnegative(value, "source_row_start");
            long stop = nonnegative(value, "source_row_stop");
            long sourceCount = nonnegative(value, "source_rows");
            long valid = nonnegative(value, "valid_molecules");
            long shardRejected = nonnegative(value, "rejected_molecules");
            require(start == expectedSourceStart && stop >= start && stop - start == sourceCount,
                    "non-contiguous supplier-cache source rows at shard " + index);
            require(Math.addExact(valid, shardRejected) == sourceCount,
                    "inconsistent counts in supplier-cache shard " + index);
            require(provenance.equals(provenance(value.path("provenance"))),
                    "checkpoint provenance differs in supplier-cache shard " + index);
            Path base = vectorPath(shardDirectory, value, "128d", valid, 128);
            Path compact = vectorPath(shardDirectory, value, "16d", valid, 16);
            JsonNode records = value.path("records");
            Path metadata = inside(shardDirectory, text(records, "path"));
            required(metadata);
            loaded.add(new Shard(index, start, stop, valid, base, compact, metadata));
            expectedSourceStart = stop;
            validTotal = Math.addExact(validTotal, valid);
            rejectedTotal = Math.addExact(rejectedTotal, shardRejected);
        }
        require(expectedSourceStart == sourceRows && validTotal == recordCount
                        && rejectedTotal == rejected,
                "supplier-cache shard totals do not match the top-level manifest");
        this.shards = List.copyOf(loaded);
    }

    @Override public String artifactType() { return ARTIFACT_TYPE; }
    @Override public long recordCount() { return recordCount; }
    @Override public int shardCount() { return shards.size(); }
    @Override public MoleculeFingerprintProvenance provenance() { return provenance; }

    @Override public MoleculeFingerprintShardReader openShard(int shardIndex,
            MoleculeFingerprintColumn column) throws IOException {
        if (shardIndex < 0 || shardIndex >= shards.size()) {
            throw new IllegalArgumentException("invalid shard index");
        }
        Shard shard = shards.get(shardIndex);
        Path vector = column == MoleculeFingerprintColumn.BASE_128 ? shard.base : shard.compact;
        return new RawFp16ShardReader(vector, shardIndex, shard.records,
                column.dimension());
    }

    @Override public Map<MoleculeVectorReference, MoleculeFingerprintMetadata> resolveMetadata(
            Collection<MoleculeVectorReference> references) throws IOException {
        Map<Integer, TreeSet<Long>> byShard = new HashMap<>();
        for (var reference : references) {
            requireReference(reference);
            byShard.computeIfAbsent(reference.shardIndex(), ignored -> new TreeSet<>())
                    .add(reference.localRow());
        }
        Map<MoleculeVectorReference, MoleculeFingerprintMetadata> result = new HashMap<>();
        for (var entry : byShard.entrySet()) {
            resolveShardMetadata(shards.get(entry.getKey()), entry.getValue(), result);
        }
        return result;
    }

    @Override public void validateCompatibility(DeepSpaceModelBundle model,
            CompactSkelSpheresModelBundle compact, boolean compactRequired) {
        require(provenance.foundationCheckpointSha256().equals(
                        model.manifest().foundationCheckpointSha256)
                        && provenance.predictorCheckpointSha256().equals(
                        model.manifest().predictorCheckpointSha256),
                "supplier cache and primary model use different checkpoints");
        if (compactRequired) {
            require(compact != null, "compact model bundle is required");
            require(provenance.foundationCheckpointSha256().equals(
                            compact.manifest().sourceFoundationCheckpointSha256)
                            && provenance.predictorCheckpointSha256().equals(
                            compact.manifest().sourcePredictorCheckpointSha256)
                            && provenance.projectionCheckpointSha256().equals(
                            compact.manifest().projectionCheckpointSha256),
                    "supplier cache and compact model use different checkpoints");
        }
    }

    private void resolveShardMetadata(Shard shard, TreeSet<Long> requested,
            Map<MoleculeVectorReference, MoleculeFingerprintMetadata> output) throws IOException {
        if (requested.isEmpty()) return;
        long maximum = requested.last();
        try (var input = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(Files.newInputStream(shard.metadata)), StandardCharsets.UTF_8),
                1024 * 1024)) {
            String header = input.readLine();
            require("local_index\tsource_row\tmolecule_id\tsmiles\tcanonical_smiles\tnum_atoms"
                    .equals(header), "unsupported supplier-cache metadata columns");
            long expected = 0;
            for (String line; expected <= maximum && (line = input.readLine()) != null; expected++) {
                String[] fields = line.split("\t", -1);
                require(fields.length == 6, "malformed supplier-cache metadata row");
                long local = parseLong(fields[0], "local_index");
                require(local == expected, "non-contiguous supplier-cache local_index");
                if (requested.contains(local)) {
                    output.put(new MoleculeVectorReference(shard.index, local),
                            new MoleculeFingerprintMetadata(parseLong(fields[1], "source_row"),
                                    parseInt(fields[5], "num_atoms"), fields[2], fields[3], fields[4]));
                }
            }
        }
        for (long local : requested) {
            require(output.containsKey(new MoleculeVectorReference(shard.index, local)),
                    "supplier-cache metadata ended before requested row " + local);
        }
    }

    private void requireReference(MoleculeVectorReference reference) {
        require(reference.shardIndex() >= 0 && reference.shardIndex() < shards.size(),
                "molecule vector reference has invalid shard");
        require(reference.localRow() >= 0
                        && reference.localRow() < shards.get(reference.shardIndex()).records,
                "molecule vector reference is outside its shard");
    }

    private static Path vectorPath(Path directory, JsonNode shard, String name,
            long records, int dimension) throws IOException {
        JsonNode value = shard.path("fingerprints").path(name);
        require("float16".equals(text(value, "dtype")), "unsupported fingerprint dtype");
        JsonNode shape = value.path("shape");
        require(shape.isArray() && shape.size() == 2 && shape.get(0).asLong(-1) == records
                        && shape.get(1).asInt(-1) == dimension,
                "unsupported fingerprint shape for " + name);
        Path path = inside(directory, text(value, "path"));
        required(path);
        require(Files.size(path) == Math.multiplyExact(records, dimension * 2L),
                "fingerprint file size does not match manifest: " + path);
        return path;
    }

    private static void requireFingerprintSummary(JsonNode value) {
        require("float16".equals(text(value.path("128d"), "dtype"))
                        && value.path("128d").path("bytes_per_molecule").asInt(-1) == 256,
                "unsupported top-level 128D fingerprint representation");
        require("float16".equals(text(value.path("16d"), "dtype"))
                        && value.path("16d").path("bytes_per_molecule").asInt(-1) == 32,
                "unsupported top-level 16D fingerprint representation");
    }

    private static MoleculeFingerprintProvenance provenance(JsonNode value) {
        String foundation = text(value, "integrated_checkpoint_sha256");
        String predictor = text(value, "predictor_checkpoint_sha256");
        String projection = text(value, "projection_checkpoint_sha256");
        require(hash(foundation) && hash(predictor) && hash(projection),
                "supplier-cache checkpoint provenance is incomplete");
        return new MoleculeFingerprintProvenance(foundation, predictor, projection);
    }

    private static Path inside(Path directory, String relative) throws IOException {
        require(relative != null && !relative.isBlank(), "missing relative cache path");
        Path path = directory.resolve(relative).normalize().toAbsolutePath();
        if (!path.startsWith(directory.toAbsolutePath().normalize())) {
            throw new IOException("cache path escapes its shard directory: " + relative);
        }
        return path;
    }

    private static JsonNode readJson(Path path) throws IOException {
        return MAPPER.readTree(path.toFile());
    }
    private static Path required(Path path) throws IOException {
        if (!Files.isRegularFile(path)) throw new IOException("missing cache file: " + path);
        return path;
    }
    private static String text(JsonNode value, String field) {
        JsonNode child = value.path(field);
        return child.isTextual() ? child.textValue() : null;
    }
    private static long nonnegative(JsonNode value, String field) {
        JsonNode child = value.path(field);
        require(child.canConvertToLong() && child.asLong() >= 0, "invalid " + field);
        return child.asLong();
    }
    private static int integer(JsonNode value, String field) {
        JsonNode child = value.path(field);
        require(child.canConvertToInt(), "invalid " + field);
        return child.asInt();
    }
    private static long parseLong(String value, String label) throws IOException {
        try { return Long.parseLong(value); }
        catch (NumberFormatException error) { throw new IOException("invalid " + label, error); }
    }
    private static int parseInt(String value, String label) throws IOException {
        try { return Integer.parseInt(value); }
        catch (NumberFormatException error) { throw new IOException("invalid " + label, error); }
    }
    private static boolean hash(String value) { return value != null && value.matches("[0-9a-f]{64}"); }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    private record Shard(int index, long sourceStart, long sourceStop, long records,
            Path base, Path compact, Path metadata) {}

    private static final class RawFp16ShardReader implements MoleculeFingerprintShardReader {
        private final FileChannel channel;
        private final ByteBuffer data;
        private final int shardIndex;
        private final long records;
        private final int dimension;
        private long row;

        private RawFp16ShardReader(Path path, int shardIndex, long records, int dimension)
                throws IOException {
            this.channel = FileChannel.open(path);
            if (channel.size() > Integer.MAX_VALUE) {
                channel.close();
                throw new IOException("supplier-cache shard exceeds Java mapping limit: " + path);
            }
            this.data = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size())
                    .order(ByteOrder.LITTLE_ENDIAN);
            this.shardIndex = shardIndex;
            this.records = records;
            this.dimension = dimension;
        }

        @Override public MoleculeVectorBatch readBatch(int maximum) {
            if (maximum < 1) throw new IllegalArgumentException("maximum must be positive");
            int count = (int) Math.min(maximum, records - row);
            float[] values = new float[Math.multiplyExact(count, dimension)];
            List<MoleculeVectorReference> references = new ArrayList<>(count);
            for (int batchRow = 0; batchRow < count; batchRow++) {
                long local = row++;
                copy(local, values, batchRow * dimension);
                references.add(new MoleculeVectorReference(shardIndex, local));
            }
            return new MoleculeVectorBatch(values, dimension, references);
        }

        @Override public float[] readVector(long localRow) {
            if (localRow < 0 || localRow >= records) throw new IndexOutOfBoundsException();
            float[] result = new float[dimension];
            copy(localRow, result, 0);
            return result;
        }

        private void copy(long localRow, float[] output, int outputOffset) {
            int offset = Math.toIntExact(localRow * dimension * 2L);
            for (int index = 0; index < dimension; index++) {
                output[outputOffset + index] = Float.float16ToFloat(data.getShort(offset + index * 2));
            }
        }

        @Override public long recordCount() { return records; }
        @Override public void close() throws IOException { channel.close(); }
    }
}
