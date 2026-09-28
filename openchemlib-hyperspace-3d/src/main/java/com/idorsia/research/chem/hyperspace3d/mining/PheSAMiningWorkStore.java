package com.idorsia.research.chem.hyperspace3d.mining;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Base64;
import java.util.HexFormat;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/** Restartable per-query work state. Query metadata is committed last. */
public final class PheSAMiningWorkStore {
    private static final String ARTIFACT_TYPE = "hyperspace-phesa-mining-work";
    private final Path directory;
    private final ObjectMapper mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);

    public PheSAMiningWorkStore(Path directory, String identitySha256, int queryCount)
            throws IOException {
        this.directory = directory.toAbsolutePath().normalize();
        Path manifestPath = this.directory.resolve("manifest.json");
        if (Files.exists(this.directory)) {
            Manifest manifest = mapper.readValue(manifestPath.toFile(), Manifest.class);
            if (!ARTIFACT_TYPE.equals(manifest.artifactType) || manifest.formatVersion != 1
                    || !identitySha256.equals(manifest.identitySha256)
                    || manifest.queryCount != queryCount) {
                throw new IOException("mining work directory belongs to a different run");
            }
        } else {
            Files.createDirectories(this.directory.resolve("queries"));
            Files.createDirectories(this.directory.resolve("descriptors"));
            Manifest manifest = new Manifest();
            manifest.identitySha256 = identitySha256; manifest.queryCount = queryCount;
            mapper.writeValue(manifestPath.toFile(), manifest);
        }
    }

    public CompletedQuery load(int index) throws IOException {
        Path metadataPath = metadata(index);
        if (!Files.isRegularFile(metadataPath)) {
            Files.deleteIfExists(query(index)); Files.deleteIfExists(descriptors(index));
            Files.deleteIfExists(temporary(query(index)));
            Files.deleteIfExists(temporary(descriptors(index)));
            Files.deleteIfExists(temporary(metadata(index)));
            return null;
        }
        QueryMetadata metadata = mapper.readValue(metadataPath.toFile(), QueryMetadata.class);
        Path queryPath = query(index); Path descriptorPath = descriptors(index);
        verify(queryPath, metadata.querySha256); verify(descriptorPath, metadata.descriptorsSha256);
        PheSAQueryShardIO.Shard shard = PheSAQueryShardIO.read(queryPath);
        List<DescriptorEntry> descriptorEntries = readDescriptors(descriptorPath);
        if (shard.queryId() != metadata.queryId || shard.records().size() != metadata.recordCount
                || descriptorEntries.size() != metadata.newDescriptorCount) {
            throw new IOException("inconsistent completed mining query " + index);
        }
        return new CompletedQuery(shard, descriptorEntries);
    }

    public void save(int index, long queryId, long screenedCount, int miningRound,
            List<PheSAQueryPairRecord> records, List<DescriptorEntry> descriptorEntries)
            throws IOException {
        if (Files.exists(metadata(index))) throw new IOException("query work already exists: " + index);
        Path queryTemporary = temporary(query(index));
        Path descriptorTemporary = temporary(descriptors(index));
        Path metadataTemporary = temporary(metadata(index));
        Files.deleteIfExists(queryTemporary); Files.deleteIfExists(descriptorTemporary);
        Files.deleteIfExists(metadataTemporary);
        PheSAQueryShardIO.write(queryTemporary, queryId, screenedCount, miningRound, records);
        writeDescriptors(descriptorTemporary, descriptorEntries);
        QueryMetadata metadata = new QueryMetadata(); metadata.queryId = queryId;
        metadata.recordCount = records.size(); metadata.newDescriptorCount = descriptorEntries.size();
        metadata.querySha256 = sha256(queryTemporary);
        metadata.descriptorsSha256 = sha256(descriptorTemporary);
        mapper.writeValue(metadataTemporary.toFile(), metadata);
        move(queryTemporary, query(index)); move(descriptorTemporary, descriptors(index));
        move(metadataTemporary, metadata(index));
    }

    private static void writeDescriptors(Path path, List<DescriptorEntry> values) throws IOException {
        try (BufferedWriter output = new BufferedWriter(new OutputStreamWriter(
                new GZIPOutputStream(Files.newOutputStream(path)), StandardCharsets.UTF_8))) {
            output.write("candidate_ordinal\tcandidate_id\tdescriptor_base64\n");
            for (DescriptorEntry value : values) {
                String encoded = value.descriptor() == null ? "" : Base64.getEncoder().encodeToString(
                        value.descriptor().getBytes(StandardCharsets.UTF_8));
                output.write(value.candidateOrdinal() + "\t" + value.candidateId() + "\t"
                        + encoded + "\n");
            }
        }
    }

    private static List<DescriptorEntry> readDescriptors(Path path) throws IOException {
        List<DescriptorEntry> result = new ArrayList<>();
        try (BufferedReader input = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(Files.newInputStream(path)), StandardCharsets.UTF_8))) {
            if (!"candidate_ordinal\tcandidate_id\tdescriptor_base64".equals(input.readLine())) {
                throw new IOException("unsupported work descriptor table");
            }
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\t", -1);
                if (fields.length != 3) throw new IOException("malformed work descriptor row");
                String descriptor = fields[2].isEmpty() ? null : new String(
                        Base64.getDecoder().decode(fields[2]), StandardCharsets.UTF_8);
                result.add(new DescriptorEntry(Integer.parseInt(fields[0]),
                        Long.parseLong(fields[1]), descriptor));
            }
        }
        return result;
    }

    private Path query(int index) { return directory.resolve(String.format("queries/q%06d.bin", index)); }
    private Path descriptors(int index) { return directory.resolve(String.format("descriptors/q%06d.tsv.gz", index)); }
    private Path metadata(int index) { return directory.resolve(String.format("queries/q%06d.json", index)); }
    private static Path temporary(Path path) { return path.resolveSibling(path.getFileName() + ".partial"); }
    private static void move(Path source, Path target) throws IOException {
        Files.move(source, target, StandardCopyOption.ATOMIC_MOVE);
    }
    private static void verify(Path path, String expected) throws IOException {
        if (!Files.isRegularFile(path) || !sha256(path).equals(expected)) {
            throw new IOException("mining work checksum mismatch: " + path);
        }
    }
    public static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (var input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int read; (read = input.read(buffer)) >= 0;) digest.update(buffer, 0, read);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException error) { throw new AssertionError(error); }
    }

    public record DescriptorEntry(int candidateOrdinal, long candidateId, String descriptor) {}
    public record CompletedQuery(PheSAQueryShardIO.Shard shard,
            List<DescriptorEntry> descriptorEntries) {}
    public static final class Manifest {
        public String artifactType = ARTIFACT_TYPE; public int formatVersion = 1;
        public String identitySha256; public int queryCount;
    }
    public static final class QueryMetadata {
        public long queryId; public int recordCount; public int newDescriptorCount;
        public String querySha256; public String descriptorsSha256;
    }
}
