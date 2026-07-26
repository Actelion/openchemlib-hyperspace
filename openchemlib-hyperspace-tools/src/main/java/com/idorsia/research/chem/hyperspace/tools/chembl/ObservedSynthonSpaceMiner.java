package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonShredder;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.SynthonReactionValidator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

/**
 * Mines large, observed-pattern pseudo-reactions from ChEMBL-style molecule tables.
 */
public final class ObservedSynthonSpaceMiner {

    private static final String SOURCE_FORMAT = "chembl-observed-synthon-space-v1";

    public Result mine(Options options) throws Exception {
        Objects.requireNonNull(options, "options");

        Map<String, CandidateReaction> candidates = new LinkedHashMap<>();
        MiningStats stats = new MiningStats();
        IDCodeParser parser = new IDCodeParser();

        try (BufferedReader reader = newReader(options.input())) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Input file is empty: " + options.input());
            }
            HeaderIndex columns = HeaderIndex.parse(header);
            int idColumn = columns.required(options.idColumn());
            int structureColumn = columns.required(options.structureColumn());

            String line;
            while ((line = reader.readLine()) != null) {
                if (options.maxMolecules() > 0 && stats.processedMolecules >= options.maxMolecules()) {
                    break;
                }
                stats.inputRows++;
                String[] fields = splitTsv(line);
                if (idColumn >= fields.length || structureColumn >= fields.length) {
                    stats.skippedRows++;
                    continue;
                }
                String moleculeId = fields[idColumn].trim();
                String idcode = fields[structureColumn].trim();
                if (moleculeId.isEmpty() || idcode.isEmpty()) {
                    stats.skippedRows++;
                    continue;
                }

                StereoMolecule molecule = new StereoMolecule();
                try {
                    parser.parse(molecule, idcode);
                    molecule.ensureHelperArrays(Molecule.cHelperRings);
                } catch (RuntimeException ex) {
                    stats.parseFailures++;
                    continue;
                }

                stats.processedMolecules++;
                mineMolecule(options, candidates, stats, moleculeId, idcode, molecule);
            }
        }

        List<CandidateReaction> selected = candidates.values().stream()
                .filter(candidate -> candidate.passes(options))
                .sorted(CandidateReaction.ORDER)
                .limit(options.topReactions())
                .collect(Collectors.toList());

        RawSynthonSpace rawSpace = buildRawSpace(options, selected, stats);
        return new Result(rawSpace, selected.stream().map(CandidateReaction::toSummary).collect(Collectors.toList()), stats);
    }

    private void mineMolecule(Options options,
                              Map<String, CandidateReaction> candidates,
                              MiningStats stats,
                              String moleculeId,
                              String sourceIdcode,
                              StereoMolecule molecule) {
        List<Integer> cutBonds = candidateCutBonds(molecule);
        int emitted = 0;
        for (int bond : cutBonds) {
            if (emitted >= options.maxCutsetsPerMolecule()) {
                return;
            }
            if (addCutset(options, candidates, stats, moleculeId, sourceIdcode, molecule, new int[]{bond})) {
                emitted++;
            }
        }

        if (options.maxCuts() < 2) {
            return;
        }
        for (int i = 0; i < cutBonds.size(); i++) {
            for (int j = i + 1; j < cutBonds.size(); j++) {
                if (emitted >= options.maxCutsetsPerMolecule()) {
                    return;
                }
                int b1 = cutBonds.get(i);
                int b2 = cutBonds.get(j);
                if (shareAtom(molecule, b1, b2) && !options.allowAdjacentTwoCuts()) {
                    continue;
                }
                if (addCutset(options, candidates, stats, moleculeId, sourceIdcode, molecule, new int[]{b1, b2})) {
                    emitted++;
                }
            }
        }
    }

    private boolean addCutset(Options options,
                              Map<String, CandidateReaction> candidates,
                              MiningStats stats,
                              String moleculeId,
                              String sourceIdcode,
                              StereoMolecule molecule,
                              int[] cutset) {
        if (!options.splitFilter().accepts(cutset.length)) {
            return false;
        }
        SynthonShredder.SplitResult split = SynthonShredder.trySplit(molecule,
                Arrays.copyOf(cutset, cutset.length),
                cutset.length + 1);
        if (split == null || split.fragments.length != cutset.length + 1) {
            stats.rejectedSplits++;
            return false;
        }

        List<StereoMolecule[]> variants = split.getAllSplitsWithUniqueConnectors();
        SplitCandidate best = null;
        for (StereoMolecule[] fragments : variants) {
            if (!passesFragmentFilters(options, fragments)) {
                continue;
            }
            if (!roundTrips(sourceIdcode, fragments)) {
                continue;
            }
            SplitCandidate candidate = SplitCandidate.from(moleculeId, sourceIdcode, split.cutset_hash, cutset.length, fragments);
            if (best == null || candidate.groupKey.compareTo(best.groupKey) < 0) {
                best = candidate;
            }
        }
        if (best == null) {
            stats.rejectedSplits++;
            return false;
        }

        int splitCount = best.splitCount;
        String signature = best.signature;
        candidates.computeIfAbsent(best.groupKey, key -> new CandidateReaction(key, splitCount, signature))
                .add(best);
        stats.acceptedSplits++;
        return true;
    }

    private RawSynthonSpace buildRawSpace(Options options,
                                          List<CandidateReaction> selected,
                                          MiningStats stats) throws Exception {
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(options.spaceName())
                .version("1.0")
                .putMetadata(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT, SOURCE_FORMAT)
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_SHORT_NAME, options.descriptorShortName())
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_BITS, Integer.toString(options.descriptorBits()))
                .putMetadata("source.input", options.input().getFileName().toString())
                .putMetadata("source.idColumn", options.idColumn())
                .putMetadata("source.structureColumn", options.structureColumn())
                .putMetadata("mining.maxCuts", Integer.toString(options.maxCuts()))
                .putMetadata("mining.splitFilter", options.splitFilter().cliValue())
                .putMetadata("mining.maxMolecules", Integer.toString(options.maxMolecules()))
                .putMetadata("mining.maxCutsetsPerMolecule", Integer.toString(options.maxCutsetsPerMolecule()))
                .putMetadata("mining.topReactions", Integer.toString(options.topReactions()))
                .putMetadata("mining.seed", Long.toString(options.seed()))
                .putMetadata("stats.inputRows", Long.toString(stats.inputRows))
                .putMetadata("stats.processedMolecules", Long.toString(stats.processedMolecules))
                .putMetadata("stats.acceptedSplits", Long.toString(stats.acceptedSplits))
                .putMetadata("stats.rejectedSplits", Long.toString(stats.rejectedSplits));

        int reactionNumber = 0;
        for (CandidateReaction candidate : selected) {
            String reactionId = String.format(Locale.ROOT,
                    "chembl_observed_%dsplit_%03d",
                    candidate.splitCount,
                    reactionNumber++);
            candidate.writeTo(builder, reactionId, options);
        }
        return builder.build();
    }

    private static boolean passesFragmentFilters(Options options, StereoMolecule[] fragments) {
        for (StereoMolecule fragment : fragments) {
            fragment.ensureHelperArrays(Molecule.cHelperNeighbours);
            BitSet connectors = SynthonSpace.computeConnectorBitSet(fragment);
            if (connectors.isEmpty()) {
                return false;
            }
            for (int atom = 0; atom < fragment.getAtoms(); atom++) {
                int atomicNo = fragment.getAtomicNo(atom);
                if (atomicNo >= 92 && atomicNo <= 95 && fragment.getConnAtoms(atom) != 1) {
                    return false;
                }
            }
            if (heavyAtomCount(fragment, false) > options.maxFragmentHeavyAtoms()) {
                return false;
            }
        }
        return true;
    }

    private static boolean roundTrips(String sourceIdcode, StereoMolecule[] fragments) {
        try {
            List<StereoMolecule> parts = Arrays.stream(fragments)
                    .map(StereoMolecule::new)
                    .collect(Collectors.toList());
            StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
            return sourceIdcode.equals(assembled.getIDCode());
        } catch (RuntimeException ex) {
            return false;
        }
    }

    private static List<Integer> candidateCutBonds(StereoMolecule molecule) {
        List<Integer> bonds = new ArrayList<>();
        molecule.ensureHelperArrays(Molecule.cHelperRings);
        for (int bond = 0; bond < molecule.getBonds(); bond++) {
            int a1 = molecule.getBondAtom(0, bond);
            int a2 = molecule.getBondAtom(1, bond);
            if (molecule.getAtomicNo(a1) <= 1 || molecule.getAtomicNo(a2) <= 1) {
                continue;
            }
            if (molecule.isRingBond(bond) || molecule.isAromaticBond(bond)) {
                continue;
            }
            if (molecule.getBondOrder(bond) != 1) {
                continue;
            }
            bonds.add(bond);
        }
        return bonds;
    }

    private static boolean shareAtom(StereoMolecule molecule, int b1, int b2) {
        int a10 = molecule.getBondAtom(0, b1);
        int a11 = molecule.getBondAtom(1, b1);
        int a20 = molecule.getBondAtom(0, b2);
        int a21 = molecule.getBondAtom(1, b2);
        return a10 == a20 || a10 == a21 || a11 == a20 || a11 == a21;
    }

    private static int heavyAtomCount(StereoMolecule molecule, boolean includeConnectors) {
        int count = 0;
        for (int atom = 0; atom < molecule.getAtoms(); atom++) {
            int atomicNo = molecule.getAtomicNo(atom);
            if (atomicNo <= 1) {
                continue;
            }
            if (!includeConnectors && atomicNo >= 92 && atomicNo <= 95) {
                continue;
            }
            count++;
        }
        return count;
    }

    private static BufferedReader newReader(Path path) throws IOException {
        InputStream in = Files.newInputStream(path);
        if (path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz")) {
            in = new GZIPInputStream(in);
        }
        return new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
    }

    private static String[] splitTsv(String line) {
        return line.split("\t", -1);
    }

    private static String connectorSetKey(StereoMolecule fragment) {
        return SynthonSpace.computeConnectorSet(fragment).stream()
                .sorted()
                .map(connector -> Integer.toString(connector - 92))
                .collect(Collectors.joining(","));
    }

    private static String fragmentSignature(StereoMolecule fragment) {
        fragment.ensureHelperArrays(Molecule.cHelperRings);
        Map<Integer, Integer> connectorPositions = SynthonSpace.computeConnectorPositions(fragment);
        return connectorPositions.keySet().stream()
                .sorted()
                .map(connector -> (connector - 92) + ":" + connectorEndpointSignature(fragment, connectorPositions.get(connector)))
                .collect(Collectors.joining(","));
    }

    private static String connectorEndpointSignature(StereoMolecule fragment, int connectorAtom) {
        if (fragment.getConnAtoms(connectorAtom) != 1) {
            return "bad";
        }
        int neighbor = fragment.getConnAtom(connectorAtom, 0);
        int bond = fragment.getConnBond(connectorAtom, 0);
        return fragment.getAtomicNo(neighbor)
                + "|bo" + fragment.getBondOrder(bond)
                + "|bt" + fragment.getBondType(bond)
                + "|ar" + fragment.isAromaticAtom(neighbor)
                + "|ring" + fragment.isRingAtom(neighbor)
                + "|deg" + heavyNeighborCount(fragment, neighbor);
    }

    private static int heavyNeighborCount(StereoMolecule fragment, int atom) {
        int count = 0;
        for (int i = 0; i < fragment.getConnAtoms(atom); i++) {
            int neighbor = fragment.getConnAtom(atom, i);
            int atomicNo = fragment.getAtomicNo(neighbor);
            if (atomicNo > 1 && (atomicNo < 92 || atomicNo > 95)) {
                count++;
            }
        }
        return count;
    }

    private static final class HeaderIndex {
        private final Map<String, Integer> index;

        private HeaderIndex(Map<String, Integer> index) {
            this.index = index;
        }

        static HeaderIndex parse(String header) {
            String[] columns = splitTsv(header);
            Map<String, Integer> index = new HashMap<>();
            for (int i = 0; i < columns.length; i++) {
                index.put(columns[i], i);
            }
            return new HeaderIndex(index);
        }

        int required(String column) {
            Integer idx = index.get(column);
            if (idx == null) {
                throw new IllegalArgumentException("Missing required input column: " + column);
            }
            return idx;
        }
    }

    private static final class SplitCandidate {
        private final String moleculeId;
        private final String sourceIdcode;
        private final String cutsetHash;
        private final int splitCount;
        private final String groupKey;
        private final String signature;
        private final List<StereoMolecule> fragments;

        private SplitCandidate(String moleculeId,
                               String sourceIdcode,
                               String cutsetHash,
                               int splitCount,
                               String groupKey,
                               String signature,
                               List<StereoMolecule> fragments) {
            this.moleculeId = moleculeId;
            this.sourceIdcode = sourceIdcode;
            this.cutsetHash = cutsetHash;
            this.splitCount = splitCount;
            this.groupKey = groupKey;
            this.signature = signature;
            this.fragments = fragments;
        }

        static SplitCandidate from(String moleculeId,
                                   String sourceIdcode,
                                   String cutsetHash,
                                   int splitCount,
                                   StereoMolecule[] fragments) {
            List<StereoMolecule> copied = Arrays.stream(fragments)
                    .map(StereoMolecule::new)
                    .collect(Collectors.toList());
            if (splitCount == 2) {
                copied = canonicalizeTwoSplitFragmentOrder(copied);
            }
            String topology = copied.stream()
                    .map(ObservedSynthonSpaceMiner::connectorSetKey)
                    .collect(Collectors.joining("|"));
            String signature = copied.stream()
                    .map(ObservedSynthonSpaceMiner::fragmentSignature)
                    .collect(Collectors.joining("|"));
            String groupKey = splitCount + "split;" + topology + ";" + signature;
            return new SplitCandidate(moleculeId, sourceIdcode, cutsetHash, splitCount, groupKey, signature, copied);
        }

        private static List<StereoMolecule> canonicalizeTwoSplitFragmentOrder(List<StereoMolecule> fragments) {
            return fragments.stream()
                    .map(fragment -> new OrderedFragment(fragmentOrderKey(fragment), fragment))
                    .sorted(Comparator.comparing(OrderedFragment::orderKey))
                    .map(OrderedFragment::fragment)
                    .collect(Collectors.toList());
        }

        private static String fragmentOrderKey(StereoMolecule fragment) {
            return connectorSetKey(fragment) + ";" + fragmentSignature(fragment) + ";" + fragment.getIDCode();
        }
    }

    private record OrderedFragment(String orderKey, StereoMolecule fragment) {
    }

    private static final class CandidateReaction {
        private static final Comparator<CandidateReaction> ORDER =
                Comparator.comparingLong(CandidateReaction::estimatedProductCount).reversed()
                        .thenComparingInt(CandidateReaction::sourceMoleculeCount).reversed()
                        .thenComparing(candidate -> candidate.groupKey);

        private final String groupKey;
        private final int splitCount;
        private final String signature;
        private final List<Map<String, SynthonEntry>> fragmentSets;
        private final Set<String> sourceMolecules = new LinkedHashSet<>();
        private final List<String> representativeCompounds = new ArrayList<>();

        private CandidateReaction(String groupKey, int splitCount, String signature) {
            this.groupKey = groupKey;
            this.splitCount = splitCount;
            this.signature = signature;
            this.fragmentSets = new ArrayList<>();
        }

        void add(SplitCandidate candidate) {
            while (fragmentSets.size() < candidate.fragments.size()) {
                fragmentSets.add(new LinkedHashMap<>());
            }
            sourceMolecules.add(candidate.moleculeId);
            if (representativeCompounds.size() < 100) {
                representativeCompounds.add(candidate.sourceIdcode);
            }
            for (int i = 0; i < candidate.fragments.size(); i++) {
                StereoMolecule fragment = candidate.fragments.get(i);
                String idcode = fragment.getIDCode();
                fragmentSets.get(i).putIfAbsent(idcode,
                        new SynthonEntry(candidate.moleculeId + "|" + candidate.cutsetHash + "|frag_" + i,
                                fragment,
                                candidate.moleculeId));
            }
        }

        boolean passes(Options options) {
            if (sourceMolecules.size() < options.minSourceMolecules()) {
                return false;
            }
            for (Map<String, SynthonEntry> set : fragmentSets) {
                if (set.size() < options.minSetSize()) {
                    return false;
                }
            }
            return true;
        }

        long estimatedProductCount() {
            long count = 1L;
            for (Map<String, SynthonEntry> set : fragmentSets) {
                count = saturatedMultiply(count, set.size());
            }
            return count;
        }

        int sourceMoleculeCount() {
            return sourceMolecules.size();
        }

        Summary toSummary() {
            return new Summary(groupKey,
                    splitCount,
                    signature,
                    fragmentSets.stream().map(Map::size).collect(Collectors.toList()),
                    sourceMolecules.size(),
                    estimatedProductCount());
        }

        void writeTo(RawSynthonSpace.Builder builder, String reactionId, Options options) throws Exception {
            Map<Integer, List<Object>> validationMap = new TreeMap<>();
            for (int fragIdx = 0; fragIdx < fragmentSets.size(); fragIdx++) {
                List<RawSynthon> rawSynthons = new ArrayList<>();
                for (SynthonEntry entry : fragmentSets.get(fragIdx).values()) {
                    String fragmentId = reactionId + "|" + entry.fragmentId();
                    RawSynthon synthon = RawSynthon.fromMolecule(reactionId, fragIdx, fragmentId, entry.fragment());
                    rawSynthons.add(synthon);
                    builder.addFragmentAttribute(reactionId, fragmentId, "source.compoundId", entry.sourceMoleculeId());
                }
                builder.addRawFragments(reactionId, fragIdx, rawSynthons);
                validationMap.put(fragIdx, rawSynthons.stream().map(RawSynthon::getIdcode).collect(Collectors.toList()));
            }
            SynthonReactionValidator.validate(validationMap);

            builder.addReactionMetadata(reactionId, "mining.groupKey", groupKey);
            builder.addReactionMetadata(reactionId, "mining.signature", signature);
            builder.addReactionMetadata(reactionId, "mining.splitType", splitCount + "-split");
            builder.addReactionMetadata(reactionId, "mining.sourceMoleculeCount", Integer.toString(sourceMolecules.size()));
            builder.addReactionMetadata(reactionId, "mining.estimatedProductCount", Long.toString(estimatedProductCount()));
            builder.addReactionMetadata(reactionId, "mining.synthonCounts", fragmentSets.stream()
                    .map(set -> Integer.toString(set.size()))
                    .collect(Collectors.joining(",")));
            builder.addRepresentativeCompounds(reactionId, sampleProducts(options));
            builder.addExampleScaffolds(reactionId, representativeCompounds.stream()
                    .limit(options.sampleProducts())
                    .collect(Collectors.toList()));
        }

        private List<String> sampleProducts(Options options) {
            if (options.sampleProducts() <= 0) {
                return Collections.emptyList();
            }
            Random random = new Random(options.seed() ^ groupKey.hashCode());
            IDCodeParser parser = new IDCodeParser();
            List<List<SynthonEntry>> sets = fragmentSets.stream()
                    .map(set -> new ArrayList<>(set.values()))
                    .collect(Collectors.toList());
            List<String> products = new ArrayList<>();
            Set<String> seen = new HashSet<>();
            int attempts = Math.max(options.sampleProducts() * 10, 20);
            for (int attempt = 0; attempt < attempts && products.size() < options.sampleProducts(); attempt++) {
                List<StereoMolecule> parts = new ArrayList<>();
                for (List<SynthonEntry> set : sets) {
                    SynthonEntry entry = set.get(random.nextInt(set.size()));
                    StereoMolecule part = new StereoMolecule();
                    parser.parse(part, entry.fragment().getIDCode());
                    parts.add(part);
                }
                try {
                    String idcode = SynthonAssembler.assembleSynthons_faster(parts).getIDCode();
                    if (seen.add(idcode)) {
                        products.add(idcode);
                    }
                } catch (RuntimeException ignored) {
                    // Invalid samples are not fatal; validator covers connector consistency.
                }
            }
            return products;
        }
    }

    private record SynthonEntry(String fragmentId, StereoMolecule fragment, String sourceMoleculeId) {
    }

    private static long saturatedMultiply(long a, long b) {
        if (a == 0 || b == 0) {
            return 0;
        }
        if (a > Long.MAX_VALUE / b) {
            return Long.MAX_VALUE;
        }
        return a * b;
    }

    public record Summary(String groupKey,
                          int splitCount,
                          String signature,
                          List<Integer> synthonCounts,
                          int sourceMoleculeCount,
                          long estimatedProductCount) {
    }

    public record Result(RawSynthonSpace rawSpace, List<Summary> reactions, MiningStats stats) {
    }

    public static final class MiningStats {
        private long inputRows;
        private long processedMolecules;
        private long skippedRows;
        private long parseFailures;
        private long acceptedSplits;
        private long rejectedSplits;

        public long inputRows() {
            return inputRows;
        }

        public long processedMolecules() {
            return processedMolecules;
        }

        public long skippedRows() {
            return skippedRows;
        }

        public long parseFailures() {
            return parseFailures;
        }

        public long acceptedSplits() {
            return acceptedSplits;
        }

        public long rejectedSplits() {
            return rejectedSplits;
        }
    }

    public record Options(Path input,
                          String spaceName,
                          String idColumn,
                          String structureColumn,
                          int maxMolecules,
                          int maxCuts,
                          int topReactions,
                          int minSetSize,
                          int minSourceMolecules,
                          int maxFragmentHeavyAtoms,
                          int maxCutsetsPerMolecule,
                          int sampleProducts,
                          long seed,
                          boolean allowAdjacentTwoCuts,
                          SplitFilter splitFilter,
                          String descriptorShortName,
                          int descriptorBits) {

        public Options {
            Objects.requireNonNull(input, "input");
            Objects.requireNonNull(spaceName, "spaceName");
            Objects.requireNonNull(idColumn, "idColumn");
            Objects.requireNonNull(structureColumn, "structureColumn");
            Objects.requireNonNull(splitFilter, "splitFilter");
            Objects.requireNonNull(descriptorShortName, "descriptorShortName");
            if (maxCuts < 1 || maxCuts > 2) {
                throw new IllegalArgumentException("maxCuts must be 1 or 2");
            }
            if (topReactions < 1 || minSetSize < 1 || minSourceMolecules < 1) {
                throw new IllegalArgumentException("reaction selection limits must be positive");
            }
            if (maxFragmentHeavyAtoms < 1 || maxCutsetsPerMolecule < 1 || sampleProducts < 0) {
                throw new IllegalArgumentException("fragment and sample limits are invalid");
            }
            if (descriptorBits < 1) {
                throw new IllegalArgumentException("descriptorBits must be positive");
            }
        }

        public static Builder builder() {
            return new Builder();
        }
    }

    public static final class Builder {
        private Path input;
        private String spaceName = "chembl_observed_synthon_space";
        private String idColumn = "chembl_id";
        private String structureColumn = "idcode";
        private int maxMolecules = 0;
        private int maxCuts = 2;
        private int topReactions = 20;
        private int minSetSize = 100;
        private int minSourceMolecules = 100;
        private int maxFragmentHeavyAtoms = 40;
        private int maxCutsetsPerMolecule = 500;
        private int sampleProducts = 25;
        private long seed = 7L;
        private boolean allowAdjacentTwoCuts = false;
        private SplitFilter splitFilter = SplitFilter.BOTH;
        private String descriptorShortName = "FragFp";
        private int descriptorBits = 1024;

        public Builder input(Path input) {
            this.input = input;
            return this;
        }

        public Builder spaceName(String spaceName) {
            this.spaceName = spaceName;
            return this;
        }

        public Builder idColumn(String idColumn) {
            this.idColumn = idColumn;
            return this;
        }

        public Builder structureColumn(String structureColumn) {
            this.structureColumn = structureColumn;
            return this;
        }

        public Builder maxMolecules(int maxMolecules) {
            this.maxMolecules = maxMolecules;
            return this;
        }

        public Builder maxCuts(int maxCuts) {
            this.maxCuts = maxCuts;
            return this;
        }

        public Builder topReactions(int topReactions) {
            this.topReactions = topReactions;
            return this;
        }

        public Builder minSetSize(int minSetSize) {
            this.minSetSize = minSetSize;
            return this;
        }

        public Builder minSourceMolecules(int minSourceMolecules) {
            this.minSourceMolecules = minSourceMolecules;
            return this;
        }

        public Builder maxFragmentHeavyAtoms(int maxFragmentHeavyAtoms) {
            this.maxFragmentHeavyAtoms = maxFragmentHeavyAtoms;
            return this;
        }

        public Builder maxCutsetsPerMolecule(int maxCutsetsPerMolecule) {
            this.maxCutsetsPerMolecule = maxCutsetsPerMolecule;
            return this;
        }

        public Builder sampleProducts(int sampleProducts) {
            this.sampleProducts = sampleProducts;
            return this;
        }

        public Builder seed(long seed) {
            this.seed = seed;
            return this;
        }

        public Builder allowAdjacentTwoCuts(boolean allowAdjacentTwoCuts) {
            this.allowAdjacentTwoCuts = allowAdjacentTwoCuts;
            return this;
        }

        public Builder splitFilter(SplitFilter splitFilter) {
            this.splitFilter = splitFilter;
            return this;
        }

        public Builder descriptorShortName(String descriptorShortName) {
            this.descriptorShortName = descriptorShortName;
            return this;
        }

        public Builder descriptorBits(int descriptorBits) {
            this.descriptorBits = descriptorBits;
            return this;
        }

        public Options build() {
            return new Options(input,
                    spaceName,
                    idColumn,
                    structureColumn,
                    maxMolecules,
                    maxCuts,
                    topReactions,
                    minSetSize,
                    minSourceMolecules,
                    maxFragmentHeavyAtoms,
                    maxCutsetsPerMolecule,
                    sampleProducts,
                    seed,
                    allowAdjacentTwoCuts,
                    splitFilter,
                    descriptorShortName,
                    descriptorBits);
        }
    }

    public enum SplitFilter {
        ONE("1"),
        TWO("2"),
        BOTH("both");

        private final String cliValue;

        SplitFilter(String cliValue) {
            this.cliValue = cliValue;
        }

        public boolean accepts(int cutCount) {
            return this == BOTH || (this == ONE && cutCount == 1) || (this == TWO && cutCount == 2);
        }

        public String cliValue() {
            return cliValue;
        }

        public static SplitFilter parse(String raw) {
            if (raw == null || raw.isBlank() || "both".equalsIgnoreCase(raw)) {
                return BOTH;
            }
            if ("1".equals(raw) || "one".equalsIgnoreCase(raw) || "1split".equalsIgnoreCase(raw)) {
                return ONE;
            }
            if ("2".equals(raw) || "two".equalsIgnoreCase(raw) || "2split".equalsIgnoreCase(raw)) {
                return TWO;
            }
            throw new IllegalArgumentException("Unsupported splitFilter: " + raw);
        }
    }
}
