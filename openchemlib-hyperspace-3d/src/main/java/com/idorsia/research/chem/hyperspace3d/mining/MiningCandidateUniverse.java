package com.idorsia.research.chem.hyperspace3d.mining;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintColumn;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSource;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
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
import java.security.NoSuchAlgorithmException;
import java.util.HexFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/** Materialized deterministic subset of a larger molecule fingerprint source. */
public final class MiningCandidateUniverse implements AutoCloseable {
    public static final String ARTIFACT_TYPE = "hyperspace-phesa-mining-universe";
    private final Path directory;
    private final Manifest manifest;
    private final List<Candidate> candidates;
    private final FileChannel vectorChannel;
    private final ByteBuffer vectors;

    private MiningCandidateUniverse(Path directory, Manifest manifest, List<Candidate> candidates,
            FileChannel vectorChannel, ByteBuffer vectors) {
        this.directory = directory; this.manifest = manifest; this.candidates = candidates;
        this.vectorChannel = vectorChannel; this.vectors = vectors;
    }

    public static Path build(Path output, MoleculeFingerprintDataSource source, int count,
            int oversample, long seed, Set<String> excludedCanonicalSmiles,
            String sourceIndexHash, String fingerprintModelHash) throws IOException {
        if (count < 1 || oversample < 0 || count + oversample < 0) {
            throw new IllegalArgumentException("invalid universe size");
        }
        Path directory = output.toAbsolutePath().normalize();
        Path partial = directory.resolveSibling(directory.getFileName() + ".partial");
        if (Files.exists(directory) || Files.exists(partial)) {
            throw new IOException("candidate universe or partial output already exists");
        }
        var sampler = new CandidateReferenceSampler(Math.addExact(count, oversample), seed);
        for (int shard = 0; shard < source.shardCount(); shard++) {
            try (var reader = source.openShard(shard, MoleculeFingerprintColumn.BASE_128)) {
                for (long row = 0; row < reader.recordCount(); row++) sampler.offer(shard, row);
            }
        }
        List<CandidateReferenceSampler.SampledReference> sampled = sampler.results();
        List<MoleculeVectorReference> references = sampled.stream()
                .map(CandidateReferenceSampler.SampledReference::reference).toList();
        Map<MoleculeVectorReference, MoleculeFingerprintMetadata> metadata =
                source.resolveMetadata(references);
        List<BuildCandidate> selected = new ArrayList<>(count);
        for (var value : sampled) {
            var molecule = metadata.get(value.reference());
            String canonical = canonical(molecule);
            if (excludedCanonicalSmiles != null && excludedCanonicalSmiles.contains(canonical)) continue;
            selected.add(new BuildCandidate(selected.size(), value.reference(), molecule, canonical,
                    properties(canonical)));
            if (selected.size() == count) break;
        }
        if (selected.size() != count) throw new IOException("universe oversample was insufficient after exclusions");
        Files.createDirectories(partial);
        write(partial, source, selected, seed, sourceIndexHash, fingerprintModelHash);
        Files.move(partial, directory, StandardCopyOption.ATOMIC_MOVE);
        return directory;
    }

    private static void write(Path directory, MoleculeFingerprintDataSource source,
            List<BuildCandidate> selected, long seed, String sourceHash, String modelHash)
            throws IOException {
        Path vectorPath = directory.resolve("descriptors.f16");
        try (FileChannel out = FileChannel.open(vectorPath, StandardOpenOption.CREATE_NEW,
                StandardOpenOption.WRITE)) {
            out.position((long) selected.size() * 256 - 1); out.write(ByteBuffer.wrap(new byte[]{0}));
            Map<Integer, List<BuildCandidate>> byShard = new HashMap<>();
            for (BuildCandidate value : selected) byShard.computeIfAbsent(
                    value.reference.shardIndex(), ignored -> new ArrayList<>()).add(value);
            ByteBuffer encoded = ByteBuffer.allocate(256).order(ByteOrder.LITTLE_ENDIAN);
            for (var entry : byShard.entrySet()) {
                try (var reader = source.openShard(entry.getKey(), MoleculeFingerprintColumn.BASE_128)) {
                    for (BuildCandidate value : entry.getValue()) {
                        encoded.clear();
                        for (float component : reader.readVector(value.reference.localRow()))
                            encoded.putShort(Float.floatToFloat16(component));
                        encoded.flip(); out.position((long) value.id * 256);
                        while (encoded.hasRemaining()) out.write(encoded);
                    }
                }
            }
        }
        Path moleculePath = directory.resolve("molecules.tsv.gz");
        try (BufferedWriter out = new BufferedWriter(new OutputStreamWriter(
                new GZIPOutputStream(Files.newOutputStream(moleculePath)),
                StandardCharsets.UTF_8))) {
            out.write("candidate_id\tshard\tlocal_row\tsource_row\tmolecule_id\tcanonical_smiles\theavy_atoms\trotatable_bonds\tsp3_fraction\tstereo_centers\n");
            selected.sort(Comparator.comparingInt(BuildCandidate::id));
            for (BuildCandidate value : selected) {
                out.write(value.id + "\t" + value.reference.shardIndex() + "\t"
                        + value.reference.localRow() + "\t" + value.metadata.sourceRow() + "\t"
                        + text(value.metadata.moleculeId()) + "\t" + text(value.canonicalSmiles)
                        + "\t" + value.properties.heavyAtoms + "\t" + value.properties.rotatableBonds
                        + "\t" + value.properties.sp3Fraction + "\t" + value.properties.stereoCenters + "\n");
            }
        }
        var manifest = new Manifest(); manifest.sourceArtifactType = source.artifactType();
        manifest.sourceIndexHash = sourceHash; manifest.fingerprintModelHash = modelHash;
        manifest.recordCount = selected.size(); manifest.samplingSeed = seed;
        manifest.descriptorsSha256 = sha256(vectorPath);
        manifest.moleculesSha256 = sha256(moleculePath);
        new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
                .writeValue(directory.resolve("manifest.json").toFile(), manifest);
        Files.createFile(directory.resolve(".complete"));
    }

    public static MiningCandidateUniverse open(Path path) throws IOException {
        Path directory = path.toAbsolutePath().normalize();
        if (!Files.isRegularFile(directory.resolve(".complete"))) throw new IOException("universe is incomplete");
        Manifest manifest = new ObjectMapper().readValue(directory.resolve("manifest.json").toFile(), Manifest.class);
        manifest.validate();
        verify(directory.resolve("descriptors.f16"), manifest.descriptorsSha256);
        verify(directory.resolve("molecules.tsv.gz"), manifest.moleculesSha256);
        List<Candidate> candidates = new ArrayList<>((int) manifest.recordCount);
        try (BufferedReader input = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(Files.newInputStream(directory.resolve("molecules.tsv.gz"))),
                StandardCharsets.UTF_8))) {
            if (!("candidate_id\tshard\tlocal_row\tsource_row\tmolecule_id\tcanonical_smiles"
                    + "\theavy_atoms\trotatable_bonds\tsp3_fraction\tstereo_centers").equals(input.readLine()))
                throw new IOException("unsupported universe molecule table");
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\t", -1);
                int id = Integer.parseInt(fields[0]);
                if (id != candidates.size()) throw new IOException("non-contiguous candidate IDs");
                candidates.add(new Candidate(id, new MoleculeVectorReference(Integer.parseInt(fields[1]),
                        Long.parseLong(fields[2])), Long.parseLong(fields[3]), fields[4], fields[5],
                        Integer.parseInt(fields[6]), Integer.parseInt(fields[7]),
                        Double.parseDouble(fields[8]), Integer.parseInt(fields[9])));
            }
        }
        if (candidates.size() != manifest.recordCount) throw new IOException("universe count mismatch");
        FileChannel channel = FileChannel.open(directory.resolve("descriptors.f16"));
        if (channel.size() != manifest.recordCount * 256) { channel.close(); throw new IOException("universe vector size mismatch"); }
        ByteBuffer vectors = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size()).order(ByteOrder.LITTLE_ENDIAN);
        return new MiningCandidateUniverse(directory, manifest, List.copyOf(candidates), channel, vectors);
    }

    public float[] fingerprint(int id) {
        if (id < 0 || id >= candidates.size()) throw new IndexOutOfBoundsException();
        float[] value = new float[128]; int offset = id * 256;
        for (int i = 0; i < 128; i++) value[i] = Float.float16ToFloat(vectors.getShort(offset + i * 2));
        return value;
    }
    public Candidate candidate(int id) { return candidates.get(id); }
    public List<Candidate> candidates() { return candidates; }
    public Manifest manifest() { return manifest; }
    @Override public void close() throws IOException { vectorChannel.close(); }

    private static void verify(Path path, String expected) throws IOException {
        String actual = sha256(path);
        if (!actual.equals(expected)) throw new IOException("candidate-universe checksum mismatch: " + path);
    }
    private static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (var input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int read; (read = input.read(buffer)) >= 0;) digest.update(buffer, 0, read);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException error) {
            throw new AssertionError(error);
        }
    }
    private static String canonical(MoleculeFingerprintMetadata value) {
        return value.canonicalSmiles() == null || value.canonicalSmiles().isBlank()
                ? value.smiles() : value.canonicalSmiles();
    }
    private static Properties properties(String smiles) throws IOException {
        try {
            StereoMolecule molecule = new StereoMolecule();
            new SmilesParser().parse(molecule, smiles);
            molecule.ensureHelperArrays(Molecule.cHelperCIP);
            boolean[] rotatable = new boolean[molecule.getBonds()];
            TorsionDB.findRotatableBonds(molecule, true, rotatable);
            int rotors = 0, sp3 = 0, stereo = 0;
            for (boolean value : rotatable) if (value) rotors++;
            for (int atom = 0; atom < molecule.getAtoms(); atom++) {
                if (molecule.getAtomPi(atom) == 0) sp3++;
                if (molecule.getAtomParity(atom) != Molecule.cAtomParityNone) stereo++;
            }
            return new Properties(molecule.getAtoms(), rotors,
                    molecule.getAtoms() == 0 ? 0.0 : (double) sp3 / molecule.getAtoms(), stereo);
        } catch (Exception error) {
            throw new IOException("cannot parse sampled candidate: " + smiles, error);
        }
    }
    private static String text(String value) { return value == null ? "" : value.replace('\t', ' '); }
    private record BuildCandidate(int id, MoleculeVectorReference reference,
            MoleculeFingerprintMetadata metadata, String canonicalSmiles, Properties properties) {}
    private record Properties(int heavyAtoms, int rotatableBonds, double sp3Fraction,
            int stereoCenters) {}
    public record Candidate(int id, MoleculeVectorReference reference, long sourceRow,
            String sourceId, String canonicalSmiles, int heavyAtoms, int rotatableBonds,
            double sp3Fraction, int stereoCenters) {}

    public static final class Manifest {
        public String artifactType = ARTIFACT_TYPE; public int formatVersion = 1;
        public String sourceArtifactType; public String sourceIndexHash; public String fingerprintModelHash;
        public String vectorDtype = "float16"; public int embeddingDim = 128;
        public long recordCount; public long samplingSeed;
        public String descriptorsSha256; public String moleculesSha256;
        public void validate() {
            if (!ARTIFACT_TYPE.equals(artifactType) || formatVersion != 1 || embeddingDim != 128
                    || !"float16".equals(vectorDtype) || recordCount < 1
                    || sourceArtifactType == null || !hash(sourceIndexHash) || !hash(fingerprintModelHash)
                    || !hash(descriptorsSha256) || !hash(moleculesSha256))
                throw new IllegalArgumentException("invalid mining universe manifest");
        }
        private static boolean hash(String value) { return value != null && value.matches("[0-9a-f]{64}"); }
    }
}
