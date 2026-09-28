package com.idorsia.research.chem.hyperspace3d.cli;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.mining.ExactPheSALabeler;
import com.idorsia.research.chem.hyperspace3d.mining.MiningCandidateUniverse;
import com.idorsia.research.chem.hyperspace3d.mining.MiningStratumSelector;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAMiningDatasetWriter;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAMiningWorkStore;
import com.idorsia.research.chem.hyperspace3d.mining.PheSAQueryPairRecord;
import com.idorsia.research.chem.hyperspace3d.mining.PredictedRankMiner;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Comparator;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/** Screens learned rank channels and exact-labels their deduplicated union. */
public final class PheSAQueryPairMiningCLI {
    private PheSAQueryPairMiningCLI() {}

    public static void main(String[] args) throws Exception {
        Path configPath = config(args); ObjectMapper mapper = new ObjectMapper();
        var config = mapper.readValue(configPath.toFile(), PheSAQueryPairMiningConfig.class);
        config.validate(); var paths = config.resolve(configPath);
        var model = DeepSpaceModelBundle.load(paths.modelBundle());
        var runtime = new DeepSpaceOnnxEnvironment(config.device(), config.runtime.cudaDeviceId);
        List<Query> queries = readQueries(paths.queries());
        Map<Integer, Long> candidateIds = new LinkedHashMap<>();
        Map<Integer, String> descriptorCache = new HashMap<>();
        List<PheSAMiningDatasetWriter.QueryEntry> queryEntries = new ArrayList<>();

        try (var universe = MiningCandidateUniverse.open(paths.universe());
             var encoder = new DeepSpaceV1Encoder(runtime, model);
             var comparator = new DeepSpaceV1Comparator(runtime, model)) {
            String fingerprintHash = sha256(model.bundleHash());
            if (!fingerprintHash.equals(universe.manifest().fingerprintModelHash)) {
                throw new IllegalArgumentException("candidate universe and model bundle differ");
            }
            String universeHash = sha256(paths.universe().resolve("manifest.json"));
            String queriesHash = sha256(paths.queries());
            String workIdentity = sha256(fingerprintHash + ":" + universeHash + ":" + queriesHash
                    + ":" + config.mining.seed + ":" + config.mining.round + ":"
                    + config.mining.maxConformers + ":" + config.mining.phesaPpWeight);
            var work = new PheSAMiningWorkStore(paths.work(), workIdentity, queries.size());

            List<float[]> queryFingerprints = new ArrayList<>();
            List<StereoMolecule> queryMolecules = new ArrayList<>();
            var tensorBuilder = new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer());
            for (Query query : queries) {
                StereoMolecule molecule = parse(query.smiles);
                queryFingerprints.add(encoder.encode(tensorBuilder.build(List.of(molecule)))[0]);
                queryMolecules.add(molecule);
            }

            long[] nextMoleculeId = {queries.size()};
            int totalIndex = model.manifest().targetIndex("phesa_total");
            int shapeIndex = model.manifest().targetIndex("phesa_shape");
            int ppIndex = model.manifest().targetIndex("phesa_pharmacophore");
            var miner = PredictedRankMiner.defaults(config.mining.seed);
            var stratumSelections = MiningStratumSelector.defaults(
                    universe.candidates(), config.mining.seed);
            ExecutorService workers = Executors.newFixedThreadPool(config.runtime.exactLabelThreads);
            ThreadLocal<ExactPheSALabeler> labelers = ThreadLocal.withInitial(() ->
                    new ExactPheSALabeler(config.mining.maxConformers,
                            config.mining.phesaPpWeight));
            try {
                for (int queryIndex = 0; queryIndex < queries.size(); queryIndex++) {
                    var completed = work.load(queryIndex);
                    if (completed != null) {
                        restoreDescriptors(completed.descriptorEntries(), candidateIds,
                                descriptorCache, nextMoleculeId);
                        Query query = queries.get(queryIndex);
                        queryEntries.add(new PheSAMiningDatasetWriter.QueryEntry(queryIndex,
                                query.uid, query.split, completed.shard().screenedCount(),
                                completed.shard().miningRound(), completed.shard().records()));
                        System.out.printf("resumed query %d/%d: %s -> %d exact pairs%n",
                                queryIndex + 1, queries.size(), query.uid,
                                completed.shard().records().size());
                        continue;
                    }

                    float[][] scores = screen(universe, comparator,
                            queryFingerprints.get(queryIndex), config.runtime.comparatorBatchSize,
                            totalIndex, shapeIndex, ppIndex);
                    List<PredictedRankMiner.Selection> selected =
                            miner.select(scores, stratumSelections);
                    ExactPheSALabeler queryLabeler = new ExactPheSALabeler(
                            config.mining.maxConformers, config.mining.phesaPpWeight);
                    String queryDescriptor = queryLabeler.encode(
                            queryLabeler.describe(queryMolecules.get(queryIndex)));

                    List<Future<LabeledChoice>> futures = new ArrayList<>(selected.size());
                    for (var choice : selected) {
                        int ordinal = choice.candidateOrdinal();
                        long candidateId = candidateIds.computeIfAbsent(ordinal,
                                ignored -> nextMoleculeId[0]++);
                        boolean newDescriptor = !descriptorCache.containsKey(ordinal);
                        String cached = descriptorCache.get(ordinal);
                        futures.add(workers.submit(labelTask(universe, choice, candidateId,
                                queryDescriptor, cached, newDescriptor, labelers,
                                config.mining.round)));
                    }
                    List<PheSAQueryPairRecord> records = new ArrayList<>(futures.size());
                    List<PheSAMiningWorkStore.DescriptorEntry> descriptorEntries = new ArrayList<>();
                    for (Future<LabeledChoice> future : futures) {
                        LabeledChoice value = future.get(); records.add(value.record);
                        if (value.newDescriptor) {
                            descriptorCache.put(value.candidateOrdinal, value.descriptor);
                            descriptorEntries.add(new PheSAMiningWorkStore.DescriptorEntry(
                                    value.candidateOrdinal, value.candidateId, value.descriptor));
                        }
                    }
                    work.save(queryIndex, queryIndex, universe.candidates().size(),
                            config.mining.round, records, descriptorEntries);
                    Query query = queries.get(queryIndex);
                    queryEntries.add(new PheSAMiningDatasetWriter.QueryEntry(queryIndex,
                            query.uid, query.split, universe.candidates().size(),
                            config.mining.round, records));
                    System.out.printf("mined query %d/%d: %s -> %d exact pairs (%d new descriptors)%n",
                            queryIndex + 1, queries.size(), query.uid, records.size(),
                            descriptorEntries.size());
                }
            } finally {
                workers.shutdownNow();
            }

            List<PheSAMiningDatasetWriter.MoleculeEntry> molecules = new ArrayList<>();
            for (int index = 0; index < queries.size(); index++) {
                Query query = queries.get(index);
                molecules.add(entry(index, query.uid, query.split, "query", query.smiles,
                        -1, -1, queryMolecules.get(index), queryFingerprints.get(index)));
            }
            candidateIds.entrySet().stream().sorted(Map.Entry.comparingByValue()).forEach(value -> {
                try {
                    var candidate = universe.candidate(value.getKey());
                    StereoMolecule structure = parse(candidate.canonicalSmiles());
                    molecules.add(entry(value.getValue(), stableUid(candidate.canonicalSmiles()),
                            "CANDIDATE_LIBRARY", "candidate", candidate.canonicalSmiles(),
                            candidate.reference().shardIndex(), candidate.sourceRow(), structure,
                            universe.fingerprint(value.getKey())));
                } catch (Exception error) { throw new CandidateMaterializationException(error); }
            });
            molecules.sort(Comparator.comparingLong(PheSAMiningDatasetWriter.MoleculeEntry::moleculeId));
            var provenance = new PheSAMiningDatasetWriter.Provenance(fingerprintHash,
                    universe.manifest().sourceIndexHash, universeHash, queriesHash,
                    ExactPheSALabeler.class.getPackage().getImplementationVersion() == null
                            ? "OpenChemLib-PheSA" : ExactPheSALabeler.class.getPackage().getImplementationVersion(),
                    config.mining.phesaPpWeight, config.mining.maxConformers,
                    config.mining.round, config.mining.seed);
            new PheSAMiningDatasetWriter().write(paths.dataset(), provenance, molecules, queryEntries);
        } catch (CandidateMaterializationException error) {
            throw (Exception) error.getCause();
        }
    }

    private static Callable<LabeledChoice> labelTask(MiningCandidateUniverse universe,
            PredictedRankMiner.Selection choice, long candidateId, String queryDescriptor,
            String cachedDescriptor, boolean newDescriptor,
            ThreadLocal<ExactPheSALabeler> labelers, int miningRound) {
        return () -> {
            int ordinal = choice.candidateOrdinal(); ExactPheSALabeler labeler = labelers.get();
            String descriptor = cachedDescriptor;
            if (newDescriptor) {
                descriptor = labeler.encode(labeler.describe(parse(
                        universe.candidate(ordinal).canonicalSmiles())));
            }
            var exact = labeler.score(labeler.decode(queryDescriptor), labeler.decode(descriptor));
            float[] predicted = choice.scores();
            var record = new PheSAQueryPairRecord(candidateId, predicted[0], predicted[1],
                    predicted[2], choice.rankTotal(), choice.rankShape(),
                    choice.rankPharmacophore(), Float.NaN, Float.NaN, Float.NaN,
                    choice.selectionMask(), miningRound,
                    PheSAQueryPairRecord.EXACT_DESCRIPTOR_FAILED).withExact(exact);
            return new LabeledChoice(ordinal, candidateId, descriptor, newDescriptor, record);
        };
    }

    private static void restoreDescriptors(List<PheSAMiningWorkStore.DescriptorEntry> entries,
            Map<Integer, Long> candidateIds, Map<Integer, String> descriptors, long[] nextId) {
        for (var value : entries) {
            Long previous = candidateIds.putIfAbsent(value.candidateOrdinal(), value.candidateId());
            if (previous != null && previous != value.candidateId()) {
                throw new IllegalArgumentException("candidate ID changed in resumed work");
            }
            descriptors.put(value.candidateOrdinal(), value.descriptor());
            nextId[0] = Math.max(nextId[0], value.candidateId() + 1);
        }
    }

    private static float[][] screen(MiningCandidateUniverse universe,
            DeepSpaceV1Comparator comparator, float[] query, int batchSize,
            int total, int shape, int pp) {
        float[][] result = new float[universe.candidates().size()][3];
        for (int start = 0; start < result.length; start += batchSize) {
            int count = Math.min(batchSize, result.length - start);
            float[] vectors = new float[count * 128];
            for (int row = 0; row < count; row++) {
                System.arraycopy(universe.fingerprint(start + row), 0,
                        vectors, row * 128, 128);
            }
            float[][] scores = comparator.compareFlat(query, vectors, count);
            for (int row = 0; row < count; row++) {
                result[start + row] = new float[]{scores[row][total],
                        scores[row][shape], scores[row][pp]};
            }
        }
        return result;
    }

    private static PheSAMiningDatasetWriter.MoleculeEntry entry(long id, String uid,
            String split, String role, String smiles, int shard, long sourceRow,
            StereoMolecule molecule, float[] fingerprint) {
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        boolean[] rotatable = new boolean[molecule.getBonds()];
        TorsionDB.findRotatableBonds(molecule, true, rotatable); int rotors = 0;
        for (boolean value : rotatable) if (value) rotors++;
        int sp3 = 0, stereo = 0;
        for (int atom = 0; atom < molecule.getAtoms(); atom++) {
            if (molecule.getAtomPi(atom) == 0) sp3++;
            if (molecule.getAtomParity(atom) != Molecule.cAtomParityNone) stereo++;
        }
        return new PheSAMiningDatasetWriter.MoleculeEntry(id, uid, split, role, (int) id,
                smiles, uid, shard, sourceRow, molecule.getAtoms(), rotors,
                molecule.getAtoms() == 0 ? 0 : (double) sp3 / molecule.getAtoms(),
                stereo, fingerprint);
    }

    private static List<Query> readQueries(Path path) throws Exception {
        List<Query> result = new ArrayList<>();
        try (BufferedReader input = Files.newBufferedReader(path)) {
            if (!"query_uid\tsplit\tsmiles".equals(input.readLine())) {
                throw new IllegalArgumentException("query TSV header must be query_uid, split, smiles");
            }
            for (String line; (line = input.readLine()) != null;) {
                String[] fields = line.split("\t", -1);
                if (fields.length != 3) throw new IllegalArgumentException("malformed query TSV row");
                result.add(new Query(fields[0], fields[1], fields[2]));
            }
        }
        if (result.isEmpty()) throw new IllegalArgumentException("query TSV is empty");
        return result;
    }
    private static StereoMolecule parse(String smiles) throws Exception {
        StereoMolecule value = new StereoMolecule(); new SmilesParser().parse(value, smiles);
        return value;
    }
    private static String stableUid(String value) throws Exception {
        return sha256(value).substring(0, 24);
    }
    private static String sha256(Path path) throws Exception {
        return PheSAMiningWorkStore.sha256(path);
    }
    private static String sha256(String value) throws Exception {
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256")
                .digest(value.getBytes(StandardCharsets.UTF_8)));
    }
    private static Path config(String[] args) {
        if (args.length != 2 || !"--config".equals(args[0])) {
            throw new IllegalArgumentException("usage: --config FILE");
        }
        return Path.of(args[1]).toAbsolutePath();
    }

    private record Query(String uid, String split, String smiles) {}
    private record LabeledChoice(int candidateOrdinal, long candidateId, String descriptor,
            boolean newDescriptor, PheSAQueryPairRecord record) {}
    private static final class CandidateMaterializationException extends RuntimeException {
        private CandidateMaterializationException(Throwable cause) { super(cause); }
    }
}
