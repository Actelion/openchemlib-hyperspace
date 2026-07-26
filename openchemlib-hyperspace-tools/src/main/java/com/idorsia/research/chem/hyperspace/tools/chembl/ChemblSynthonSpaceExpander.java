package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.LSHProvider;
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
import java.util.Base64;
import java.util.BitSet;
import java.util.Collection;
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
 * Expands observed ChEMBL pseudo-reaction rawspaces with connector-compatible synthons
 * mined from an additional molecule table.
 */
public final class ChemblSynthonSpaceExpander {

    private static final String SOURCE_FORMAT = "chembl-observed-synthon-space-expanded-v1";
    private static final int MAX_ANCHORS_PER_SET = 256;

    public Result expand(RawSynthonSpace seed, Options options) throws Exception {
        Objects.requireNonNull(seed, "seed");
        Objects.requireNonNull(options, "options");

        CachedDescriptorProvider descriptorProvider = new CachedDescriptorProvider(options.descriptorShortName());
        Map<SetKey, SeedSetProfile> profiles = buildProfiles(seed, options, descriptorProvider);
        Catalog catalog = mineCatalog(options, descriptorProvider);
        RawSynthonSpace expanded = buildExpandedRawSpace(seed, profiles, catalog, options);
        return new Result(expanded, catalog.stats, summarize(seed, expanded, profiles));
    }

    private Map<SetKey, SeedSetProfile> buildProfiles(RawSynthonSpace seed,
                                                       Options options,
                                                       CachedDescriptorProvider descriptorProvider) {
        Map<SetKey, SeedSetProfile> profiles = new LinkedHashMap<>();
        IDCodeParser parser = new IDCodeParser();
        seed.getReactions().forEach((reactionId, reactionData) ->
                reactionData.getRawFragmentSets().forEach((fragIdx, synthons) -> {
                    if (synthons == null || synthons.isEmpty()) {
                        return;
                    }
                    SeedSetProfile profile = new SeedSetProfile(reactionId, fragIdx);
                    for (RawSynthon synthon : synthons) {
                        StereoMolecule molecule = new StereoMolecule();
                        parser.parse(molecule, synthon.getIdcode());
                        molecule.ensureHelperArrays(Molecule.cHelperRings);
                        FragmentFeatures features = FragmentFeatures.from(molecule,
                                options.connectorRegionSize(),
                                descriptorProvider);
                        if (features.connectorCount != 1 && features.connectorCount != 2) {
                            continue;
                        }
                        profile.addSeed(synthon, molecule, features);
                    }
                    if (profile.connectorCount == 1 || profile.connectorCount == 2) {
                        profile.finishAnchors(options.seed());
                        profiles.put(new SetKey(reactionId, fragIdx), profile);
                    }
                }));
        return profiles;
    }

    private Catalog mineCatalog(Options options, CachedDescriptorProvider descriptorProvider) throws IOException {
        Catalog catalog = new Catalog();
        IDCodeParser parser = new IDCodeParser();
        long startNanos = System.nanoTime();
        try (BufferedReader reader = newReader(options.chemblInput())) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Input file is empty: " + options.chemblInput());
            }
            HeaderIndex columns = HeaderIndex.parse(header);
            int idColumn = columns.required(options.idColumn());
            int structureColumn = columns.required(options.structureColumn());

            String line;
            while ((line = reader.readLine()) != null) {
                if (options.catalogMaxMolecules() > 0
                        && catalog.stats.processedMolecules >= options.catalogMaxMolecules()) {
                    break;
                }
                catalog.stats.inputRows++;
                String[] fields = splitTsv(line);
                if (idColumn >= fields.length || structureColumn >= fields.length) {
                    catalog.stats.skippedRows++;
                    continue;
                }
                String moleculeId = fields[idColumn].trim();
                String idcode = fields[structureColumn].trim();
                if (moleculeId.isEmpty() || idcode.isEmpty()) {
                    catalog.stats.skippedRows++;
                    continue;
                }

                StereoMolecule molecule = new StereoMolecule();
                try {
                    parser.parse(molecule, idcode);
                    molecule.ensureHelperArrays(Molecule.cHelperRings);
                } catch (RuntimeException ex) {
                    catalog.stats.parseFailures++;
                    continue;
                }

                catalog.stats.processedMolecules++;
                mineMolecule(options, descriptorProvider, catalog, moleculeId, molecule);
                logCatalogProgress(options, catalog, startNanos);
            }
        }
        catalog.finish();
        return catalog;
    }

    private void mineMolecule(Options options,
                              CachedDescriptorProvider descriptorProvider,
                              Catalog catalog,
                              String moleculeId,
                              StereoMolecule molecule) {
        List<Integer> cutBonds = candidateCutBonds(molecule);
        int emitted = 0;
        if (options.catalogSplitFilter().accepts(1)) {
            for (int bond : cutBonds) {
                if (emitted >= options.maxCutsetsPerMolecule()) {
                    return;
                }
                emitted += addCutset(options, descriptorProvider, catalog, moleculeId, molecule, new int[]{bond}) ? 1 : 0;
            }
        }
        if (!options.catalogSplitFilter().accepts(2)) {
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
                emitted += addCutset(options, descriptorProvider, catalog, moleculeId, molecule, new int[]{b1, b2}) ? 1 : 0;
            }
        }
    }

    private boolean addCutset(Options options,
                              CachedDescriptorProvider descriptorProvider,
                              Catalog catalog,
                              String moleculeId,
                              StereoMolecule molecule,
                              int[] cutset) {
        SynthonShredder.SplitResult split = SynthonShredder.trySplit(molecule,
                Arrays.copyOf(cutset, cutset.length),
                cutset.length + 1);
        if (split == null || split.fragments.length != cutset.length + 1) {
            catalog.stats.rejectedSplits++;
            return false;
        }
        boolean acceptedAny = false;
        for (StereoMolecule[] variant : split.getAllSplitsWithUniqueConnectors()) {
            for (StereoMolecule fragment : variant) {
                StereoMolecule copy = new StereoMolecule(fragment);
                copy.ensureHelperArrays(Molecule.cHelperRings);
                FragmentFeatures features = FragmentFeatures.from(copy,
                        options.connectorRegionSize(),
                        descriptorProvider);
                if (!passesCatalogFragment(options, copy, features)) {
                    continue;
                }
                catalog.add(new CatalogEntry(moleculeId, copy.getIDCode(), features));
                acceptedAny = true;
            }
        }
        if (acceptedAny) {
            catalog.stats.acceptedSplits++;
        } else {
            catalog.stats.rejectedSplits++;
        }
        return acceptedAny;
    }

    private RawSynthonSpace buildExpandedRawSpace(RawSynthonSpace seed,
                                                  Map<SetKey, SeedSetProfile> profiles,
                                                  Catalog catalog,
                                                  Options options) throws Exception {
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(seed.getName() + "_expanded")
                .version(seed.getVersion());
        seed.getMetadata().forEach(builder::putMetadata);
        builder.putMetadata(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT, SOURCE_FORMAT);
        builder.putMetadata("expansion.source.format", "chembl-synthon-catalog");
        builder.putMetadata("expansion.source.input", options.chemblInput().getFileName().toString());
        builder.putMetadata("expansion.catalog.maxMolecules", Integer.toString(options.catalogMaxMolecules()));
        builder.putMetadata("expansion.catalog.splitFilter", options.catalogSplitFilter().cliValue());
        builder.putMetadata("expansion.targetConnectorCount", options.targetConnectorCount().cliValue());
        builder.putMetadata("expansion.maxAddedPerSet", Integer.toString(options.maxAddedPerSet()));
        builder.putMetadata("expansion.maxCandidateHeavyAtoms", Integer.toString(options.maxCandidateHeavyAtoms()));
        builder.putMetadata("expansion.minFragFpSimilarity", Double.toString(options.minFragFpSimilarity()));
        builder.putMetadata("expansion.connectorRegionSize", Integer.toString(options.connectorRegionSize()));
        builder.putMetadata("expansion.mode", options.mode().cliValue());
        builder.putMetadata("expansion.seed", Long.toString(options.seed()));
        builder.putMetadata("expansion.scoring.similarityWeight", Double.toString(options.similarityWeight()));
        builder.putMetadata("expansion.scoring.smallnessWeight", Double.toString(options.smallnessWeight()));
        builder.putMetadata("expansion.scoring.occurrenceWeight", Double.toString(options.occurrenceWeight()));
        builder.putMetadata("expansion.scoring.sizeBiasPower", Double.toString(options.sizeBiasPower()));

        for (Map.Entry<String, RawSynthonSpace.ReactionData> reactionEntry : seed.getReactions().entrySet()) {
            String reactionId = reactionEntry.getKey();
            RawSynthonSpace.ReactionData reactionData = reactionEntry.getValue();
            int addedForReaction = 0;
            List<Integer> beforeCounts = new ArrayList<>();
            List<Integer> afterCounts = new ArrayList<>();
            Map<Integer, List<Object>> validationMap = new TreeMap<>();

            for (Integer fragIdx : reactionData.getRawFragmentSets().keySet().stream().sorted().collect(Collectors.toList())) {
                List<RawSynthon> synthons = reactionData.getRawFragmentSets().get(fragIdx);
                beforeCounts.add(synthons.size());
                SeedSetProfile profile = profiles.get(new SetKey(reactionId, fragIdx));
                List<ScoredCandidate> additions = profile == null
                        ? Collections.emptyList()
                        : selectAdditions(profile, catalog, options);
                List<RawSynthon> expandedSet = new ArrayList<>(synthons);
                for (ScoredCandidate scored : additions) {
                    RawSynthon synthon = createExpandedSynthon(reactionId, fragIdx, profile, scored.entry);
                    expandedSet.add(synthon);
                    addExpansionAttributes(builder, reactionId, synthon.getFragmentId(), scored);
                }
                afterCounts.add(expandedSet.size());
                builder.addRawFragments(reactionId, fragIdx, expandedSet);
                validationMap.put(fragIdx, expandedSet.stream().map(RawSynthon::getIdcode).collect(Collectors.toList()));
                System.out.println("[Expansion] " + reactionId + ":" + fragIdx
                        + " connectors=" + (profile == null ? "n/a" : profile.connectorCount)
                        + " " + synthons.size() + " -> " + expandedSet.size()
                        + " added=" + additions.size());
            }

            reactionData.getExampleScaffolds().forEach(idcode ->
                    builder.addExampleScaffolds(reactionId, Collections.singletonList(idcode)));
            reactionData.getPartialAssemblies().forEach((idx, assemblies) ->
                    builder.addPartialAssemblies(reactionId, idx, assemblies));
            builder.addRepresentativeCompounds(reactionId, reactionData.getRepresentativeCompounds());
            reactionData.getDescriptors().forEach((key, value) ->
                    builder.addReactionDescriptor(reactionId, key, value));
            reactionData.getReactionMetadata().forEach((key, value) ->
                    builder.addReactionMetadata(reactionId, key, value));
            reactionData.getFragmentAttributes().forEach((fragmentId, attributes) -> {
                if (attributes != null) {
                    attributes.forEach((key, value) ->
                            builder.addFragmentAttribute(reactionId, fragmentId, key, value));
                }
            });

            for (int i = 0; i < beforeCounts.size(); i++) {
                addedForReaction += afterCounts.get(i) - beforeCounts.get(i);
            }
            SynthonReactionValidator.validate(validationMap);
            builder.addReactionMetadata(reactionId, "expansion.beforeSynthonCounts", joinInts(beforeCounts));
            builder.addReactionMetadata(reactionId, "expansion.afterSynthonCounts", joinInts(afterCounts));
            builder.addReactionMetadata(reactionId, "expansion.addedSynthons", Integer.toString(addedForReaction));
            builder.addReactionMetadata(reactionId, "expansion.estimatedProductCount",
                    Long.toString(estimatedProductCount(afterCounts)));
            System.out.println("[Expansion] " + reactionId
                    + " " + beforeCounts + " -> " + afterCounts
                    + " added=" + addedForReaction
                    + " products~" + estimatedProductCount(afterCounts));
        }
        return builder.build();
    }

    private List<ScoredCandidate> selectAdditions(SeedSetProfile profile,
                                                  Catalog catalog,
                                                  Options options) {
        if (!options.targetConnectorCount().accepts(profile.connectorCount)) {
            return Collections.emptyList();
        }
        Collection<CatalogEntry> bucket = profile.connectorCount == 1
                ? catalog.getOneConnector(profile.endpointSignatures)
                : catalog.getTwoConnector(profile.pairKeys);
        if (bucket.isEmpty()) {
            return Collections.emptyList();
        }
        int maxOccurrence = bucket.stream().mapToInt(entry -> entry.occurrence).max().orElse(1);
        List<ScoredCandidate> scored = new ArrayList<>();
        for (CatalogEntry entry : bucket) {
            if (profile.existingIdcodes.contains(entry.idcode)) {
                continue;
            }
            if (entry.features.heavyAtoms > options.maxCandidateHeavyAtoms()) {
                continue;
            }
            if (!profile.accepts(entry.features)) {
                continue;
            }
            double similarity = bestSimilarity(entry.features.fullFp, profile.anchorFps);
            if (similarity < options.minFragFpSimilarity()) {
                continue;
            }
            double smallness = 1.0 - ((double) entry.features.heavyAtoms / options.maxCandidateHeavyAtoms());
            smallness = Math.max(0.0, Math.min(1.0, smallness));
            smallness = Math.pow(smallness, options.sizeBiasPower());
            double occurrence = Math.log1p(entry.occurrence) / Math.log1p(maxOccurrence);
            double score = score(options, similarity, smallness, occurrence);
            scored.add(new ScoredCandidate(entry, score, similarity, smallness, occurrence));
        }
        scored.sort(Comparator.comparingDouble(ScoredCandidate::score).reversed()
                .thenComparing(scoredCandidate -> scoredCandidate.entry.features.heavyAtoms)
                .thenComparing(scoredCandidate -> scoredCandidate.entry.idcode));
        if (scored.size() > options.maxAddedPerSet()) {
            return new ArrayList<>(scored.subList(0, options.maxAddedPerSet()));
        }
        return scored;
    }

    private RawSynthon createExpandedSynthon(String reactionId,
                                             int fragIdx,
                                             SeedSetProfile profile,
                                             CatalogEntry entry) {
        StereoMolecule molecule = parse(entry.idcode);
        relabelConnectors(molecule, profile);
        String fragmentId = reactionId + "|expanded|" + fragIdx + "|" + stableFragmentSuffix(entry.idcode);
        return RawSynthon.fromMolecule(reactionId, fragIdx, fragmentId, molecule);
    }

    private void addExpansionAttributes(RawSynthonSpace.Builder builder,
                                        String reactionId,
                                        String fragmentId,
                                        ScoredCandidate scored) {
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.sourceCompoundId", scored.entry.sourceMoleculeId);
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.score",
                String.format(Locale.ROOT, "%.6f", scored.score));
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.fragFpSimilarity",
                String.format(Locale.ROOT, "%.6f", scored.fragFpSimilarity));
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.smallnessScore",
                String.format(Locale.ROOT, "%.6f", scored.smallnessScore));
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.occurrenceScore",
                String.format(Locale.ROOT, "%.6f", scored.occurrenceScore));
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.heavyAtoms",
                Integer.toString(scored.entry.features.heavyAtoms));
        builder.addFragmentAttribute(reactionId, fragmentId, "expansion.occurrence",
                Integer.toString(scored.entry.occurrence));
    }

    private List<ReactionSummary> summarize(RawSynthonSpace seed,
                                            RawSynthonSpace expanded,
                                            Map<SetKey, SeedSetProfile> profiles) {
        List<ReactionSummary> summaries = new ArrayList<>();
        expanded.getReactions().forEach((reactionId, data) -> {
            RawSynthonSpace.ReactionData before = seed.getReactions().get(reactionId);
            List<Integer> beforeCounts = before == null ? Collections.emptyList()
                    : before.getRawFragmentSets().keySet().stream().sorted()
                    .map(idx -> before.getRawFragmentSets().get(idx).size())
                    .collect(Collectors.toList());
            List<Integer> afterCounts = data.getRawFragmentSets().keySet().stream().sorted()
                    .map(idx -> data.getRawFragmentSets().get(idx).size())
                    .collect(Collectors.toList());
            int added = 0;
            for (int i = 0; i < Math.min(beforeCounts.size(), afterCounts.size()); i++) {
                added += afterCounts.get(i) - beforeCounts.get(i);
            }
            summaries.add(new ReactionSummary(reactionId, beforeCounts, afterCounts, added,
                    estimatedProductCount(afterCounts)));
        });
        return summaries;
    }

    private static boolean passesCatalogFragment(Options options,
                                                 StereoMolecule fragment,
                                                 FragmentFeatures features) {
        if (features.connectorCount != 1 && features.connectorCount != 2) {
            return false;
        }
        if (features.heavyAtoms > options.maxCandidateHeavyAtoms()) {
            return false;
        }
        for (int atom = 0; atom < fragment.getAtoms(); atom++) {
            int atomicNo = fragment.getAtomicNo(atom);
            if (atomicNo >= 92 && atomicNo <= 95 && fragment.getConnAtoms(atom) != 1) {
                return false;
            }
        }
        return true;
    }

    public static double score(Options options, double similarity, double smallness, double occurrence) {
        double weightSum = options.similarityWeight() + options.smallnessWeight() + options.occurrenceWeight();
        return ((options.similarityWeight() * similarity)
                + (options.smallnessWeight() * smallness)
                + (options.occurrenceWeight() * occurrence)) / weightSum;
    }

    private static void logCatalogProgress(Options options, Catalog catalog, long startNanos) {
        int interval = options.progressInterval();
        if (interval <= 0 || catalog.stats.processedMolecules == 0
                || catalog.stats.processedMolecules % interval != 0) {
            return;
        }
        long elapsedMillis = (System.nanoTime() - startNanos) / 1_000_000L;
        System.out.println("[Catalog] processed=" + catalog.stats.processedMolecules
                + " acceptedSplits=" + catalog.stats.acceptedSplits
                + " rejectedSplits=" + catalog.stats.rejectedSplits
                + " uniqueFragments=" + catalog.entriesByIdentity.size()
                + " elapsed=" + formatElapsed(elapsedMillis));
    }

    private static String formatElapsed(long millis) {
        long seconds = millis / 1_000L;
        long minutes = seconds / 60L;
        long remainingSeconds = seconds % 60L;
        return String.format(Locale.ROOT, "%d:%02d", minutes, remainingSeconds);
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

    private static String connectorEndpointSignature(StereoMolecule fragment, int connectorAtom) {
        fragment.ensureHelperArrays(Molecule.cHelperRings);
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

    private static int heavyAtomCount(StereoMolecule molecule) {
        int count = 0;
        for (int atom = 0; atom < molecule.getAtoms(); atom++) {
            int atomicNo = molecule.getAtomicNo(atom);
            if (atomicNo > 1 && (atomicNo < 92 || atomicNo > 95)) {
                count++;
            }
        }
        return count;
    }

    private static String canonicalPair(List<String> signatures) {
        List<String> sorted = new ArrayList<>(signatures);
        Collections.sort(sorted);
        return String.join("||", sorted);
    }

    private static String canonicalIdentity(StereoMolecule molecule, FragmentFeatures features) {
        StereoMolecule copy = new StereoMolecule(molecule);
        for (int atom = 0; atom < copy.getAtoms(); atom++) {
            int atomicNo = copy.getAtomicNo(atom);
            if (atomicNo >= 92 && atomicNo <= 95) {
                copy.setAtomicNo(atom, 92);
            }
        }
        copy.ensureHelperArrays(Molecule.cHelperCIP);
        return features.connectorCount + "|" + features.pairKey + "|" + copy.getIDCode();
    }

    private static double bestSimilarity(BitSet candidate, List<BitSet> anchors) {
        double best = 0.0;
        for (BitSet anchor : anchors) {
            best = Math.max(best, LSHProvider.tanimoto_similarity(candidate, anchor));
        }
        return best;
    }

    private static void relabelConnectors(StereoMolecule molecule, SeedSetProfile profile) {
        Map<Integer, Integer> positions = SynthonSpace.computeConnectorPositions(molecule);
        if (positions.size() != profile.targetConnectors.size()) {
            throw new IllegalArgumentException("Connector count mismatch for expansion candidate");
        }
        if (positions.size() == 1) {
            int atom = positions.values().iterator().next();
            molecule.setAtomicNo(atom, profile.targetConnectors.get(0));
            molecule.ensureHelperArrays(Molecule.cHelperCIP);
            return;
        }

        List<Integer> sourceConnectors = new ArrayList<>(positions.keySet());
        Collections.sort(sourceConnectors);
        List<Integer> targetLabels = resolveTwoConnectorAssignment(molecule, sourceConnectors, profile);
        for (int i = 0; i < sourceConnectors.size(); i++) {
            molecule.setAtomicNo(positions.get(sourceConnectors.get(i)), targetLabels.get(i));
        }
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
    }

    private static List<Integer> resolveTwoConnectorAssignment(StereoMolecule molecule,
                                                               List<Integer> sourceConnectors,
                                                               SeedSetProfile profile) {
        String sig0 = connectorEndpointSignature(molecule, SynthonSpace.computeConnectorPositions(molecule).get(sourceConnectors.get(0)));
        String sig1 = connectorEndpointSignature(molecule, SynthonSpace.computeConnectorPositions(molecule).get(sourceConnectors.get(1)));
        List<Integer> targets = profile.targetConnectors;
        boolean direct = profile.targetLabelSignatures.getOrDefault(targets.get(0), Collections.emptySet()).contains(sig0)
                && profile.targetLabelSignatures.getOrDefault(targets.get(1), Collections.emptySet()).contains(sig1);
        if (direct) {
            return targets;
        }
        boolean swapped = profile.targetLabelSignatures.getOrDefault(targets.get(0), Collections.emptySet()).contains(sig1)
                && profile.targetLabelSignatures.getOrDefault(targets.get(1), Collections.emptySet()).contains(sig0);
        if (swapped) {
            return Arrays.asList(targets.get(1), targets.get(0));
        }
        return targets;
    }

    private static StereoMolecule parse(String idcode) {
        StereoMolecule molecule = new StereoMolecule();
        new IDCodeParser().parse(molecule, idcode);
        molecule.ensureHelperArrays(Molecule.cHelperRings);
        return molecule;
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

    private static String joinInts(List<Integer> values) {
        return values.stream().map(Object::toString).collect(Collectors.joining(","));
    }

    private static long estimatedProductCount(List<Integer> counts) {
        long result = 1L;
        for (int count : counts) {
            if (result > Long.MAX_VALUE / Math.max(1, count)) {
                return Long.MAX_VALUE;
            }
            result *= count;
        }
        return result;
    }

    private static String bitSetKey(BitSet bitSet) {
        return Base64.getEncoder().encodeToString(bitSet.toByteArray());
    }

    private static String stableFragmentSuffix(String value) {
        return Integer.toUnsignedString(value.hashCode(), 36);
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

    private static final class FragmentFeatures {
        private final int connectorCount;
        private final int heavyAtoms;
        private final Map<Integer, String> signatureByConnector;
        private final Set<String> endpointSignatures;
        private final String pairKey;
        private final String connectorProximalKey;
        private final BitSet fullFp;

        private FragmentFeatures(int connectorCount,
                                 int heavyAtoms,
                                 Map<Integer, String> signatureByConnector,
                                 Set<String> endpointSignatures,
                                 String pairKey,
                                 String connectorProximalKey,
                                 BitSet fullFp) {
            this.connectorCount = connectorCount;
            this.heavyAtoms = heavyAtoms;
            this.signatureByConnector = signatureByConnector;
            this.endpointSignatures = endpointSignatures;
            this.pairKey = pairKey;
            this.connectorProximalKey = connectorProximalKey;
            this.fullFp = fullFp;
        }

        static FragmentFeatures from(StereoMolecule molecule,
                                     int connectorRegionSize,
                                     CachedDescriptorProvider descriptorProvider) {
            molecule.ensureHelperArrays(Molecule.cHelperRings);
            Map<Integer, Integer> connectorPositions = SynthonSpace.computeConnectorPositions(molecule);
            Map<Integer, String> signatureByConnector = new TreeMap<>();
            for (Map.Entry<Integer, Integer> entry : connectorPositions.entrySet()) {
                signatureByConnector.put(entry.getKey(), connectorEndpointSignature(molecule, entry.getValue()));
            }
            Set<String> endpointSignatures = new LinkedHashSet<>(signatureByConnector.values());
            String pairKey = canonicalPair(new ArrayList<>(signatureByConnector.values()));
            StereoMolecule connectorProximal = SynthonSpace.createConnectorProximalFragment(molecule, connectorRegionSize);
            BitSet connectorProximalFp = descriptorProvider.getFP_cached(connectorProximal);
            BitSet fullFp = descriptorProvider.getFP_cached(molecule);
            return new FragmentFeatures(connectorPositions.size(),
                    heavyAtomCount(molecule),
                    signatureByConnector,
                    endpointSignatures,
                    pairKey,
                    bitSetKey(connectorProximalFp),
                    fullFp);
        }
    }

    private record SetKey(String reactionId, int fragIdx) {
    }

    private static final class SeedSetProfile {
        private final String reactionId;
        private final int fragIdx;
        private int connectorCount = -1;
        private final Set<String> existingIdcodes = new HashSet<>();
        private final Set<String> endpointSignatures = new LinkedHashSet<>();
        private final Set<String> pairKeys = new LinkedHashSet<>();
        private final Set<String> connectorProximalKeys = new LinkedHashSet<>();
        private final Map<Integer, Set<String>> targetLabelSignatures = new TreeMap<>();
        private final List<Integer> targetConnectors = new ArrayList<>();
        private final List<Anchor> anchors = new ArrayList<>();
        private List<BitSet> anchorFps = Collections.emptyList();

        private SeedSetProfile(String reactionId, int fragIdx) {
            this.reactionId = reactionId;
            this.fragIdx = fragIdx;
        }

        void addSeed(RawSynthon synthon, StereoMolecule molecule, FragmentFeatures features) {
            if (connectorCount < 0) {
                connectorCount = features.connectorCount;
                targetConnectors.addAll(features.signatureByConnector.keySet());
                Collections.sort(targetConnectors);
            }
            if (features.connectorCount != connectorCount) {
                return;
            }
            existingIdcodes.add(synthon.getIdcode());
            endpointSignatures.addAll(features.endpointSignatures);
            pairKeys.add(features.pairKey);
            connectorProximalKeys.add(features.connectorProximalKey);
            features.signatureByConnector.forEach((connector, signature) ->
                    targetLabelSignatures.computeIfAbsent(connector, key -> new LinkedHashSet<>()).add(signature));
            anchors.add(new Anchor(synthon.getIdcode(), features.fullFp));
        }

        void finishAnchors(long seed) {
            anchors.sort(Comparator.comparing(Anchor::idcode));
            if (anchors.size() > MAX_ANCHORS_PER_SET) {
                Collections.shuffle(anchors, new Random(seed ^ reactionId.hashCode() ^ fragIdx));
                anchors.subList(MAX_ANCHORS_PER_SET, anchors.size()).clear();
            }
            anchorFps = anchors.stream().map(Anchor::fingerprint).collect(Collectors.toList());
        }

        boolean accepts(FragmentFeatures features) {
            if (features.connectorCount != connectorCount) {
                return false;
            }
            if (connectorCount == 1) {
                return !Collections.disjoint(endpointSignatures, features.endpointSignatures);
            }
            return pairKeys.contains(features.pairKey);
        }
    }

    private record Anchor(String idcode, BitSet fingerprint) {
    }

    private static final class CatalogEntry {
        private final String sourceMoleculeId;
        private final String idcode;
        private final FragmentFeatures features;
        private int occurrence = 1;

        private CatalogEntry(String sourceMoleculeId, String idcode, FragmentFeatures features) {
            this.sourceMoleculeId = sourceMoleculeId;
            this.idcode = idcode;
            this.features = features;
        }
    }

    private static final class Catalog {
        private final Map<String, CatalogEntry> entriesByIdentity = new LinkedHashMap<>();
        private final Map<String, List<CatalogEntry>> oneByEndpoint = new HashMap<>();
        private final Map<String, List<CatalogEntry>> twoByPair = new HashMap<>();
        private final CatalogStats stats = new CatalogStats();

        void add(CatalogEntry entry) {
            StereoMolecule molecule = parse(entry.idcode);
            String identity = canonicalIdentity(molecule, entry.features);
            CatalogEntry existing = entriesByIdentity.get(identity);
            if (existing != null) {
                existing.occurrence++;
                return;
            }
            entriesByIdentity.put(identity, entry);
        }

        void finish() {
            for (CatalogEntry entry : entriesByIdentity.values()) {
                if (entry.features.connectorCount == 1) {
                    for (String signature : entry.features.endpointSignatures) {
                        oneByEndpoint.computeIfAbsent(signature, key -> new ArrayList<>()).add(entry);
                    }
                } else if (entry.features.connectorCount == 2) {
                    twoByPair.computeIfAbsent(entry.features.pairKey, key -> new ArrayList<>()).add(entry);
                }
            }
            stats.catalogFragments = entriesByIdentity.size();
            stats.oneConnectorFragments = oneByEndpoint.values().stream().mapToLong(List::size).sum();
            stats.twoConnectorFragments = twoByPair.values().stream().mapToLong(List::size).sum();
        }

        Collection<CatalogEntry> getOneConnector(Set<String> endpointSignatures) {
            Map<String, CatalogEntry> result = new LinkedHashMap<>();
            for (String signature : endpointSignatures) {
                for (CatalogEntry entry : oneByEndpoint.getOrDefault(signature, Collections.emptyList())) {
                    result.put(entry.idcode, entry);
                }
            }
            return result.values();
        }

        Collection<CatalogEntry> getTwoConnector(Set<String> pairKeys) {
            Map<String, CatalogEntry> result = new LinkedHashMap<>();
            for (String pairKey : pairKeys) {
                for (CatalogEntry entry : twoByPair.getOrDefault(pairKey, Collections.emptyList())) {
                    result.put(entry.idcode, entry);
                }
            }
            return result.values();
        }
    }

    private record ScoredCandidate(CatalogEntry entry,
                                   double score,
                                   double fragFpSimilarity,
                                   double smallnessScore,
                                   double occurrenceScore) {
    }

    public record Result(RawSynthonSpace rawSpace,
                         CatalogStats catalogStats,
                         List<ReactionSummary> reactions) {
    }

    public static final class CatalogStats {
        private long inputRows;
        private long processedMolecules;
        private long skippedRows;
        private long parseFailures;
        private long acceptedSplits;
        private long rejectedSplits;
        private long catalogFragments;
        private long oneConnectorFragments;
        private long twoConnectorFragments;

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

        public long catalogFragments() {
            return catalogFragments;
        }

        public long oneConnectorFragments() {
            return oneConnectorFragments;
        }

        public long twoConnectorFragments() {
            return twoConnectorFragments;
        }
    }

    public record ReactionSummary(String reactionId,
                                  List<Integer> beforeCounts,
                                  List<Integer> afterCounts,
                                  int addedSynthons,
                                  long estimatedProductCount) {
    }

    public record Options(Path chemblInput,
                          String idColumn,
                          String structureColumn,
                          int catalogMaxMolecules,
                          ObservedSynthonSpaceMiner.SplitFilter catalogSplitFilter,
                          TargetConnectorCount targetConnectorCount,
                          int maxCutsetsPerMolecule,
                          int maxAddedPerSet,
                          int maxCandidateHeavyAtoms,
                          double minFragFpSimilarity,
                          int connectorRegionSize,
                          long seed,
                          boolean allowAdjacentTwoCuts,
                          ExpansionMode mode,
                          String descriptorShortName,
                          double similarityWeight,
                          double smallnessWeight,
                          double occurrenceWeight,
                          double sizeBiasPower,
                          int progressInterval) {

        public Options {
            Objects.requireNonNull(chemblInput, "chemblInput");
            Objects.requireNonNull(idColumn, "idColumn");
            Objects.requireNonNull(structureColumn, "structureColumn");
            Objects.requireNonNull(catalogSplitFilter, "catalogSplitFilter");
            Objects.requireNonNull(targetConnectorCount, "targetConnectorCount");
            Objects.requireNonNull(mode, "mode");
            Objects.requireNonNull(descriptorShortName, "descriptorShortName");
            if (catalogMaxMolecules < 0 || maxCutsetsPerMolecule < 1 || maxAddedPerSet < 0 || progressInterval < 0) {
                throw new IllegalArgumentException("catalog and expansion limits are invalid");
            }
            if (maxCandidateHeavyAtoms < 1 || connectorRegionSize < 1) {
                throw new IllegalArgumentException("fragment limits are invalid");
            }
            if (minFragFpSimilarity < 0.0 || minFragFpSimilarity > 1.0) {
                throw new IllegalArgumentException("minFragFpSimilarity must be in [0,1]");
            }
            if (similarityWeight < 0.0 || smallnessWeight < 0.0 || occurrenceWeight < 0.0
                    || similarityWeight + smallnessWeight + occurrenceWeight <= 0.0) {
                throw new IllegalArgumentException("scoring weights must be non-negative and not all zero");
            }
            if (sizeBiasPower <= 0.0) {
                throw new IllegalArgumentException("sizeBiasPower must be positive");
            }
        }

        public static Builder builder() {
            return new Builder();
        }
    }

    public static final class Builder {
        private Path chemblInput;
        private String idColumn = "chembl_id";
        private String structureColumn = "idcode";
        private int catalogMaxMolecules = 10_000;
        private ObservedSynthonSpaceMiner.SplitFilter catalogSplitFilter = ObservedSynthonSpaceMiner.SplitFilter.BOTH;
        private TargetConnectorCount targetConnectorCount = TargetConnectorCount.BOTH;
        private int maxCutsetsPerMolecule = 500;
        private int maxAddedPerSet = 1_000;
        private int maxCandidateHeavyAtoms = 24;
        private double minFragFpSimilarity = 0.25;
        private int connectorRegionSize = SynthonSpace.CONNECTOR_REGION_SIZE;
        private long seed = 7L;
        private boolean allowAdjacentTwoCuts = false;
        private ExpansionMode mode = ExpansionMode.STRICT;
        private String descriptorShortName = "FragFp";
        private double similarityWeight = 0.55;
        private double smallnessWeight = 0.35;
        private double occurrenceWeight = 0.10;
        private double sizeBiasPower = 1.0;
        private int progressInterval = 1_000;

        public Builder chemblInput(Path chemblInput) {
            this.chemblInput = chemblInput;
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

        public Builder catalogMaxMolecules(int catalogMaxMolecules) {
            this.catalogMaxMolecules = catalogMaxMolecules;
            return this;
        }

        public Builder catalogSplitFilter(ObservedSynthonSpaceMiner.SplitFilter catalogSplitFilter) {
            this.catalogSplitFilter = catalogSplitFilter;
            return this;
        }

        public Builder targetConnectorCount(TargetConnectorCount targetConnectorCount) {
            this.targetConnectorCount = targetConnectorCount;
            return this;
        }

        public Builder maxCutsetsPerMolecule(int maxCutsetsPerMolecule) {
            this.maxCutsetsPerMolecule = maxCutsetsPerMolecule;
            return this;
        }

        public Builder maxAddedPerSet(int maxAddedPerSet) {
            this.maxAddedPerSet = maxAddedPerSet;
            return this;
        }

        public Builder maxCandidateHeavyAtoms(int maxCandidateHeavyAtoms) {
            this.maxCandidateHeavyAtoms = maxCandidateHeavyAtoms;
            return this;
        }

        public Builder minFragFpSimilarity(double minFragFpSimilarity) {
            this.minFragFpSimilarity = minFragFpSimilarity;
            return this;
        }

        public Builder connectorRegionSize(int connectorRegionSize) {
            this.connectorRegionSize = connectorRegionSize;
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

        public Builder mode(ExpansionMode mode) {
            this.mode = mode;
            return this;
        }

        public Builder descriptorShortName(String descriptorShortName) {
            this.descriptorShortName = descriptorShortName;
            return this;
        }

        public Builder similarityWeight(double similarityWeight) {
            this.similarityWeight = similarityWeight;
            return this;
        }

        public Builder smallnessWeight(double smallnessWeight) {
            this.smallnessWeight = smallnessWeight;
            return this;
        }

        public Builder occurrenceWeight(double occurrenceWeight) {
            this.occurrenceWeight = occurrenceWeight;
            return this;
        }

        public Builder sizeBiasPower(double sizeBiasPower) {
            this.sizeBiasPower = sizeBiasPower;
            return this;
        }

        public Builder progressInterval(int progressInterval) {
            this.progressInterval = progressInterval;
            return this;
        }

        public Options build() {
            return new Options(chemblInput,
                    idColumn,
                    structureColumn,
                    catalogMaxMolecules,
                    catalogSplitFilter,
                    targetConnectorCount,
                    maxCutsetsPerMolecule,
                    maxAddedPerSet,
                    maxCandidateHeavyAtoms,
                    minFragFpSimilarity,
                    connectorRegionSize,
                    seed,
                    allowAdjacentTwoCuts,
                    mode,
                    descriptorShortName,
                    similarityWeight,
                    smallnessWeight,
                    occurrenceWeight,
                    sizeBiasPower,
                    progressInterval);
        }
    }

    public enum TargetConnectorCount {
        ONE("1"),
        TWO("2"),
        BOTH("both");

        private final String cliValue;

        TargetConnectorCount(String cliValue) {
            this.cliValue = cliValue;
        }

        public boolean accepts(int connectorCount) {
            return this == BOTH || (this == ONE && connectorCount == 1) || (this == TWO && connectorCount == 2);
        }

        public String cliValue() {
            return cliValue;
        }

        public static TargetConnectorCount parse(String raw) {
            if (raw == null || raw.isBlank() || "both".equalsIgnoreCase(raw)) {
                return BOTH;
            }
            if ("1".equals(raw) || "one".equalsIgnoreCase(raw) || "1connector".equalsIgnoreCase(raw)) {
                return ONE;
            }
            if ("2".equals(raw) || "two".equalsIgnoreCase(raw) || "2connector".equalsIgnoreCase(raw)) {
                return TWO;
            }
            throw new IllegalArgumentException("Unsupported targetConnectorCount: " + raw);
        }
    }

    public enum ExpansionMode {
        STRICT("strict");

        private final String cliValue;

        ExpansionMode(String cliValue) {
            this.cliValue = cliValue;
        }

        public String cliValue() {
            return cliValue;
        }

        public static ExpansionMode parse(String raw) {
            if (raw == null || raw.isBlank() || "strict".equalsIgnoreCase(raw)) {
                return STRICT;
            }
            throw new IllegalArgumentException("Unsupported expansion mode: " + raw);
        }
    }
}
