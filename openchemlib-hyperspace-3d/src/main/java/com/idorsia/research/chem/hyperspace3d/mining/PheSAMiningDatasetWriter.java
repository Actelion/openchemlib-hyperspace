package com.idorsia.research.chem.hyperspace3d.mining;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

/** Atomic finalizer for the dependency-free Java/Python mining cache. */
public final class PheSAMiningDatasetWriter {
    private final ObjectMapper mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);

    public Path write(Path output, Provenance provenance, List<MoleculeEntry> molecules,
            List<QueryEntry> queries) throws IOException {
        Path directory = output.toAbsolutePath().normalize();
        if (Files.exists(directory)) throw new IOException("output already exists: " + directory);
        Path partial = directory.resolveSibling(directory.getFileName() + ".partial");
        if (Files.exists(partial)) throw new IOException("partial output already exists: " + partial);
        Files.createDirectories(partial.resolve("queries"));
        try {
            validate(molecules, queries);
            Path moleculePath = partial.resolve("molecules.tsv.gz");
            Path descriptorPath = partial.resolve("descriptors.f16");
            writeMolecules(moleculePath, molecules);
            writeDescriptors(descriptorPath, molecules);
            List<QueryIndexEntry> index = writeQueries(partial, queries);
            Path indexPath = partial.resolve("query-index.tsv");
            writeQueryIndex(indexPath, index);

            long pairs = queries.stream().mapToLong(value -> value.records().size()).sum();
            var manifest = new PheSAMiningManifest();
            manifest.phesaPpWeight = provenance.phesaPpWeight();
            manifest.maxConformers = provenance.maxConformers();
            manifest.miningRound = provenance.miningRound();
            manifest.samplingSeed = provenance.samplingSeed();
            manifest.fingerprintModelHash = provenance.fingerprintModelHash();
            manifest.sourceIndexHash = provenance.sourceIndexHash();
            manifest.candidateUniverseHash = provenance.candidateUniverseHash();
            manifest.queryRegistryHash = provenance.queryRegistryHash();
            manifest.oclVersion = provenance.oclVersion();
            manifest.moleculeCount = molecules.size();
            manifest.queryCount = queries.size();
            manifest.pairCount = pairs;
            manifest.rankBoundaries = doubles(PredictedRankMiner.DEFAULT_BOUNDARIES);
            manifest.rankQuotas = integers(PredictedRankMiner.DEFAULT_QUOTAS);
            manifest.randomQuota = PredictedRankMiner.DEFAULT_RANDOM_QUOTA;
            manifest.chemistryStratumQuotas = Map.of(
                    "FLEXIBLE", MiningStratumSelector.DEFAULT_FLEXIBLE_QUOTA,
                    "HIGH_SP3", MiningStratumSelector.DEFAULT_HIGH_SP3_QUOTA,
                    "STEREOCHEMICAL", MiningStratumSelector.DEFAULT_STEREOCHEMICAL_QUOTA,
                    "SIZE_EXTREME_PER_TAIL", MiningStratumSelector.DEFAULT_SIZE_EXTREME_QUOTA_PER_TAIL);
            manifest.selectionMask = MiningSelection.manifestNames();
            manifest.files = new LinkedHashMap<>();
            addFile(manifest.files, "molecules", partial, moleculePath);
            addFile(manifest.files, "descriptors", partial, descriptorPath);
            addFile(manifest.files, "query_index", partial, indexPath);
            manifest.validate();
            mapper.writeValue(partial.resolve("manifest.json").toFile(), manifest);
            Files.createFile(partial.resolve(".complete"));
            try {
                return Files.move(partial, directory, StandardCopyOption.ATOMIC_MOVE);
            } catch (IOException error) {
                throw new IOException("cannot atomically finalize mining cache " + directory, error);
            }
        } catch (Throwable error) {
            if (error instanceof IOException value) throw value;
            if (error instanceof RuntimeException value) throw value;
            throw new IOException(error);
        }
    }

    private static void validate(List<MoleculeEntry> molecules, List<QueryEntry> queries) {
        if (molecules == null || queries == null || molecules.isEmpty() || queries.isEmpty()) {
            throw new IllegalArgumentException("molecules and queries must not be empty");
        }
        Map<Long, MoleculeEntry> byId = new LinkedHashMap<>();
        for (int index = 0; index < molecules.size(); index++) {
            MoleculeEntry value = molecules.get(index);
            if (value.moleculeId() < 0 || value.descriptorOffset() != index
                    || value.fingerprint() == null || value.fingerprint().length != 128
                    || value.moleculeUid() == null || value.moleculeUid().isBlank()
                    || value.canonicalSmiles() == null || value.canonicalSmiles().isBlank()
                    || byId.put(value.moleculeId(), value) != null) {
                throw new IllegalArgumentException("invalid or duplicate molecule entry at " + index);
            }
            for (float component : value.fingerprint()) {
                if (!Float.isFinite(component)) throw new IllegalArgumentException("non-finite fingerprint");
            }
        }
        for (QueryEntry query : queries) {
            if (!byId.containsKey(query.queryId()) || query.screenedCount() < query.records().size()
                    || query.split() == null || query.split().isBlank()) {
                throw new IllegalArgumentException("invalid query entry: " + query.queryId());
            }
            for (PheSAQueryPairRecord record : query.records()) {
                if (!byId.containsKey(record.candidateId())) {
                    throw new IllegalArgumentException("query references unknown candidate");
                }
            }
        }
    }

    private static void writeMolecules(Path path, List<MoleculeEntry> values) throws IOException {
        try (var out = new BufferedWriter(new OutputStreamWriter(
                new GZIPOutputStream(Files.newOutputStream(path), 1024 * 1024),
                StandardCharsets.UTF_8), 1024 * 1024)) {
            out.write("molecule_id\tmolecule_uid\tsplit\trole\tdescriptor_offset\tcanonical_smiles"
                    + "\tsource_id\tsource_shard\tsource_row\theavy_atoms\trotatable_bonds"
                    + "\tsp3_fraction\tstereo_centers\n");
            for (MoleculeEntry value : values) {
                out.write(value.moleculeId() + "\t" + text(value.moleculeUid()) + "\t"
                        + text(value.split()) + "\t" + text(value.role()) + "\t"
                        + value.descriptorOffset() + "\t" + text(value.canonicalSmiles()) + "\t"
                        + text(value.sourceId()) + "\t" + value.sourceShard() + "\t"
                        + value.sourceRow() + "\t" + value.heavyAtoms() + "\t"
                        + value.rotatableBonds() + "\t" + value.sp3Fraction() + "\t"
                        + value.stereoCenters() + "\n");
            }
        }
    }

    private static void writeDescriptors(Path path, List<MoleculeEntry> values) throws IOException {
        try (FileChannel out = FileChannel.open(path, StandardOpenOption.CREATE_NEW,
                StandardOpenOption.WRITE)) {
            ByteBuffer buffer = ByteBuffer.allocate(128 * 2 * 4096).order(ByteOrder.LITTLE_ENDIAN);
            for (MoleculeEntry value : values) {
                if (buffer.remaining() < 256) { buffer.flip(); writeFully(out, buffer); buffer.clear(); }
                for (float component : value.fingerprint()) buffer.putShort(Float.floatToFloat16(component));
            }
            buffer.flip(); writeFully(out, buffer);
        }
    }

    private static List<QueryIndexEntry> writeQueries(Path directory, List<QueryEntry> queries)
            throws IOException {
        List<QueryIndexEntry> result = new ArrayList<>(queries.size());
        for (int index = 0; index < queries.size(); index++) {
            QueryEntry value = queries.get(index);
            String relative = String.format("queries/q%06d.bin", index);
            Path path = directory.resolve(relative);
            PheSAQueryShardIO.write(path, value.queryId(), value.screenedCount(),
                    value.miningRound(), value.records());
            result.add(new QueryIndexEntry(value.queryId(), value.queryUid(), value.split(),
                    relative, value.records().size(), value.screenedCount(), sha256(path)));
        }
        return result;
    }

    private static void writeQueryIndex(Path path, List<QueryIndexEntry> values) throws IOException {
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("query_id\tquery_uid\tsplit\tpath\trow_count\tscreened_count\tsha256\n");
            for (QueryIndexEntry value : values) {
                out.write(value.queryId + "\t" + text(value.queryUid) + "\t" + text(value.split)
                        + "\t" + value.path + "\t" + value.rowCount + "\t"
                        + value.screenedCount + "\t" + value.sha256 + "\n");
            }
        }
    }

    private static void addFile(Map<String, PheSAMiningManifest.FileEntry> files, String name,
            Path directory, Path path) throws IOException {
        files.put(name, new PheSAMiningManifest.FileEntry(
                directory.relativize(path).toString(), sha256(path), Files.size(path)));
    }
    private static String text(String value) {
        return value == null ? "" : value.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
    }
    private static List<Double> doubles(double[] values) {
        List<Double> result = new ArrayList<>(values.length);
        for (double value : values) result.add(value);
        return result;
    }
    private static List<Integer> integers(int[] values) {
        List<Integer> result = new ArrayList<>(values.length);
        for (int value : values) result.add(value);
        return result;
    }
    private static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (var input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int read; (read = input.read(buffer)) >= 0;) digest.update(buffer, 0, read);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (java.security.NoSuchAlgorithmException error) {
            throw new AssertionError(error);
        }
    }
    private static void writeFully(FileChannel channel, ByteBuffer value) throws IOException {
        while (value.hasRemaining()) channel.write(value);
    }

    public record Provenance(String fingerprintModelHash, String sourceIndexHash,
            String candidateUniverseHash, String queryRegistryHash, String oclVersion, double phesaPpWeight,
            int maxConformers, int miningRound, long samplingSeed) {}
    public record MoleculeEntry(long moleculeId, String moleculeUid, String split, String role,
            int descriptorOffset, String canonicalSmiles, String sourceId, int sourceShard,
            long sourceRow, int heavyAtoms, int rotatableBonds, double sp3Fraction,
            int stereoCenters, float[] fingerprint) {}
    public record QueryEntry(long queryId, String queryUid, String split, long screenedCount,
            int miningRound, List<PheSAQueryPairRecord> records) {}
    private record QueryIndexEntry(long queryId, String queryUid, String split, String path,
            long rowCount, long screenedCount, String sha256) {}
}
