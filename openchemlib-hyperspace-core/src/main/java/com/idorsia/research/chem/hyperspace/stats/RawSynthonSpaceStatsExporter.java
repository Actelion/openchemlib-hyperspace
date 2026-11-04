package com.idorsia.research.chem.hyperspace.stats;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public final class RawSynthonSpaceStatsExporter {
    private RawSynthonSpaceStatsExporter() {
    }

    public static RawSynthonSpaceStatsReport analyze(RawSynthonSpace rawSpace,
                                                     RawSynthonSpaceStatsOptions options) {
        Objects.requireNonNull(rawSpace, "rawSpace");
        RawSynthonSpaceStatsOptions effectiveOptions = options == null
                ? RawSynthonSpaceStatsOptions.builder().build()
                : options;

        List<Map.Entry<String, RawSynthonSpace.ReactionData>> reactions = new ArrayList<>(rawSpace.getReactions().entrySet());
        reactions.sort(Map.Entry.comparingByKey());
        if (effectiveOptions.getMaxReactions() > 0 && reactions.size() > effectiveOptions.getMaxReactions()) {
            reactions = new ArrayList<>(reactions.subList(0, effectiveOptions.getMaxReactions()));
        }

        ConcurrentMap<String, Integer> atomCountCache = new ConcurrentHashMap<>();
        List<ReactionReport> reactionReports = analyzeReactions(reactions, effectiveOptions, atomCountCache);
        reactionReports.sort(Comparator.comparing(report -> report.stats.reactionId()));

        List<RawSynthonSpaceStatsReport.ReactionStats> reactionStats = reactionReports.stream()
                .map(report -> report.stats)
                .collect(Collectors.toList());
        List<RawSynthonSpaceStatsReport.SynthonSetStats> synthonSetStats = reactionReports.stream()
                .flatMap(report -> report.synthonSetStats.stream())
                .collect(Collectors.toList());
        List<RawSynthonSpaceStatsReport.ExampleProduct> examples = reactionReports.stream()
                .flatMap(report -> report.examples.stream())
                .collect(Collectors.toList());
        List<RawSynthonSpaceStatsReport.SourceStats> sourceStats = buildSourceStats(reactionReports);
        RawSynthonSpaceStatsReport.SpaceSummary summary = buildSpaceSummary(rawSpace, reactionStats,
                synthonSetStats, reactionReports, effectiveOptions);
        return new RawSynthonSpaceStatsReport(summary, reactionStats, sourceStats, synthonSetStats, examples);
    }

    private static List<ReactionReport> analyzeReactions(List<Map.Entry<String, RawSynthonSpace.ReactionData>> reactions,
                                                         RawSynthonSpaceStatsOptions options,
                                                         ConcurrentMap<String, Integer> atomCountCache) {
        if (reactions.isEmpty()) {
            return Collections.emptyList();
        }
        ExecutorService executor = Executors.newFixedThreadPool(options.getThreads(), new StatsThreadFactory());
        try {
            List<Future<ReactionReport>> futures = new ArrayList<>();
            for (Map.Entry<String, RawSynthonSpace.ReactionData> reaction : reactions) {
                futures.add(executor.submit(() -> analyzeReaction(reaction.getKey(), reaction.getValue(), options, atomCountCache)));
            }

            List<ReactionReport> results = new ArrayList<>(futures.size());
            int completed = 0;
            for (Future<ReactionReport> future : futures) {
                try {
                    ReactionReport report = future.get();
                    results.add(report);
                    completed++;
                    options.getProgressListener().onReactionCompleted(completed, futures.size(), report.stats.reactionId());
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new IllegalStateException("Rawspace stats export interrupted", e);
                } catch (ExecutionException e) {
                    throw new IllegalStateException("Rawspace stats export failed", e.getCause());
                }
            }
            return results;
        } finally {
            executor.shutdownNow();
        }
    }

    private static ReactionReport analyzeReaction(String reactionId,
                                                  RawSynthonSpace.ReactionData reactionData,
                                                  RawSynthonSpaceStatsOptions options,
                                                  ConcurrentMap<String, Integer> atomCountCache) {
        Map<Integer, List<RawSynthon>> sets = nonEmptySets(reactionData.getRawFragmentSets());
        String sourceSpace = findMetadataValue(reactionData.getReactionMetadata(), "source.spaceName", ".spaceName");
        String sourcePath = findMetadataValue(reactionData.getReactionMetadata(), "source.spacePath", ".spacePath");
        String originalReactionId = findMetadataValue(reactionData.getReactionMetadata(), "source.originalReactionId", ".originalReactionId");

        List<Integer> setIndexes = new ArrayList<>(sets.keySet());
        Collections.sort(setIndexes);
        long totalSynthons = 0L;
        Set<String> uniqueIdcodes = new HashSet<>();
        IntStats synthonAtomStats = new IntStats();
        List<RawSynthonSpaceStatsReport.SynthonSetStats> setStats = new ArrayList<>();
        BigInteger productCount = BigInteger.ONE;
        List<String> setSizeTokens = new ArrayList<>();

        for (Integer setIndex : setIndexes) {
            List<RawSynthon> synthons = sets.get(setIndex);
            productCount = productCount.multiply(BigInteger.valueOf(synthons.size()));
            setSizeTokens.add(setIndex + ":" + synthons.size());
            totalSynthons += synthons.size();
            Set<String> uniqueSetIdcodes = new HashSet<>();
            IntStats setAtomStats = new IntStats();
            for (RawSynthon synthon : synthons) {
                uniqueIdcodes.add(synthon.getIdcode());
                uniqueSetIdcodes.add(synthon.getIdcode());
                int atoms = atomCount(synthon.getIdcode(), atomCountCache);
                synthonAtomStats.add(atoms);
                setAtomStats.add(atoms);
            }
            setStats.add(new RawSynthonSpaceStatsReport.SynthonSetStats(
                    reactionId,
                    sourceSpace,
                    setIndex,
                    synthons.size(),
                    uniqueSetIdcodes.size(),
                    setAtomStats.min(),
                    setAtomStats.mean(),
                    setAtomStats.max()));
        }
        if (sets.isEmpty()) {
            productCount = BigInteger.ZERO;
        }

        SamplingResult sampling = sampleProducts(reactionId, sourceSpace, sets, options);
        RawSynthonSpaceStatsReport.ReactionStats stats = new RawSynthonSpaceStatsReport.ReactionStats(
                reactionId,
                sourceSpace,
                sourcePath,
                originalReactionId,
                sets.size(),
                String.join(",", setSizeTokens),
                productCount,
                totalSynthons,
                uniqueIdcodes.size(),
                synthonAtomStats.min(),
                synthonAtomStats.mean(),
                synthonAtomStats.max(),
                sampling.productStats.count(),
                sampling.failures,
                sampling.productStats.atomStats.min(),
                sampling.productStats.atomStats.mean(),
                sampling.productStats.atomStats.max(),
                sampling.productStats.rotatableStats.min(),
                sampling.productStats.rotatableStats.mean(),
                sampling.productStats.rotatableStats.max());
        return new ReactionReport(stats, setStats, sampling.examples, uniqueIdcodes);
    }

    private static Map<Integer, List<RawSynthon>> nonEmptySets(Map<Integer, List<RawSynthon>> rawSets) {
        Map<Integer, List<RawSynthon>> result = new LinkedHashMap<>();
        rawSets.forEach((setIndex, synthons) -> {
            if (synthons != null && !synthons.isEmpty()) {
                result.put(setIndex, synthons);
            }
        });
        return result;
    }

    private static SamplingResult sampleProducts(String reactionId,
                                                 String sourceSpace,
                                                 Map<Integer, List<RawSynthon>> sets,
                                                 RawSynthonSpaceStatsOptions options) {
        int productSamples = options.getProductSamplesPerReaction();
        int examplesRequested = options.getExamplesPerReaction();
        if (sets.size() < 2 || (productSamples <= 0 && examplesRequested <= 0)) {
            return new SamplingResult(new ProductStats(), 0, Collections.emptyList());
        }

        List<Integer> setIndexes = new ArrayList<>(sets.keySet());
        Collections.sort(setIndexes);
        Random random = new Random(options.getSeed() ^ reactionId.hashCode());
        IDCodeParser parser = new IDCodeParser();
        ProductStats productStats = new ProductStats();
        List<RawSynthonSpaceStatsReport.ExampleProduct> examples = new ArrayList<>();
        int failures = 0;
        int attempts = Math.max(productSamples, examplesRequested);

        for (int attempt = 0; attempt < attempts; attempt++) {
            try {
                List<StereoMolecule> parts = new ArrayList<>(setIndexes.size());
                List<String> fragmentIds = new ArrayList<>(setIndexes.size());
                for (int setIndex : setIndexes) {
                    List<RawSynthon> pool = sets.get(setIndex);
                    RawSynthon synthon = pool.get(random.nextInt(pool.size()));
                    StereoMolecule part = new StereoMolecule();
                    parser.parse(part, synthon.getIdcode());
                    part.ensureHelperArrays(Molecule.cHelperCIP);
                    parts.add(part);
                    fragmentIds.add(synthon.getFragmentId());
                }
                StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
                assembled.ensureHelperArrays(Molecule.cHelperCIP);
                int atoms = assembled.getAtoms();
                int rotatable = countRotatableBonds(assembled);
                if (attempt < productSamples) {
                    productStats.add(atoms, rotatable);
                }
                if (examples.size() < examplesRequested) {
                    examples.add(new RawSynthonSpaceStatsReport.ExampleProduct(
                            reactionId,
                            sourceSpace,
                            examples.size() + 1,
                            String.join(",", fragmentIds),
                            assembled.getIDCode(),
                            atoms,
                            rotatable));
                }
            } catch (RuntimeException ex) {
                failures++;
            }
        }
        return new SamplingResult(productStats, failures, examples);
    }

    private static int atomCount(String idcode, ConcurrentMap<String, Integer> atomCountCache) {
        return atomCountCache.computeIfAbsent(idcode, key -> {
            try {
                StereoMolecule mol = new StereoMolecule();
                new IDCodeParser().parse(mol, key);
                return mol.getAtoms();
            } catch (RuntimeException ex) {
                return -1;
            }
        });
    }

    private static int countRotatableBonds(StereoMolecule molecule) {
        molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
        boolean[] rotatable = new boolean[molecule.getBonds()];
        TorsionDB.findRotatableBonds(molecule, true, rotatable);
        int count = 0;
        for (boolean value : rotatable) {
            if (value) {
                count++;
            }
        }
        return count;
    }

    private static List<RawSynthonSpaceStatsReport.SourceStats> buildSourceStats(List<ReactionReport> reactionReports) {
        Map<String, SourceAccumulator> accumulators = new TreeMap<>();
        for (ReactionReport report : reactionReports) {
            RawSynthonSpaceStatsReport.ReactionStats stats = report.stats;
            String source = isBlank(stats.sourceSpace()) ? "" : stats.sourceSpace();
            SourceAccumulator accumulator = accumulators.computeIfAbsent(source, SourceAccumulator::new);
            accumulator.reactionCount++;
            accumulator.synthonSetCount += stats.synthonSetCount();
            accumulator.totalSynthons += stats.totalSynthons();
            accumulator.totalProductCount = accumulator.totalProductCount.add(stats.productCount());
            accumulator.uniqueIdcodes.addAll(report.uniqueIdcodes);
            if (stats.synthonSetCount() == 2) {
                accumulator.twoSetReactionCount++;
            } else if (stats.synthonSetCount() == 3) {
                accumulator.threeSetReactionCount++;
            } else {
                accumulator.otherReactionCount++;
            }
        }
        return accumulators.values().stream()
                .map(SourceAccumulator::toStats)
                .collect(Collectors.toList());
    }

    private static RawSynthonSpaceStatsReport.SpaceSummary buildSpaceSummary(RawSynthonSpace rawSpace,
                                                                             List<RawSynthonSpaceStatsReport.ReactionStats> reactionStats,
                                                                             List<RawSynthonSpaceStatsReport.SynthonSetStats> synthonSetStats,
                                                                             List<ReactionReport> reactionReports,
                                                                             RawSynthonSpaceStatsOptions options) {
        int twoSet = 0;
        int threeSet = 0;
        int other = 0;
        long totalSynthons = 0L;
        BigInteger totalProductCount = BigInteger.ZERO;
        Set<String> uniqueIdcodes = new HashSet<>();
        for (ReactionReport report : reactionReports) {
            RawSynthonSpaceStatsReport.ReactionStats stats = report.stats;
            if (stats.synthonSetCount() == 2) {
                twoSet++;
            } else if (stats.synthonSetCount() == 3) {
                threeSet++;
            } else {
                other++;
            }
            totalSynthons += stats.totalSynthons();
            totalProductCount = totalProductCount.add(stats.productCount());
            uniqueIdcodes.addAll(report.uniqueIdcodes);
        }
        return new RawSynthonSpaceStatsReport.SpaceSummary(
                rawSpace.getName(),
                rawSpace.getVersion(),
                rawSpace.getMetadata(),
                reactionStats.size(),
                twoSet,
                threeSet,
                other,
                synthonSetStats.size(),
                totalSynthons,
                uniqueIdcodes.size(),
                totalProductCount,
                options.getExamplesPerReaction(),
                options.getProductSamplesPerReaction(),
                options.getSeed());
    }

    private static String findMetadataValue(Map<String, String> metadata, String preferredKey, String fallbackSuffix) {
        String preferred = metadata.get(preferredKey);
        if (!isBlank(preferred)) {
            return preferred;
        }
        for (Map.Entry<String, String> entry : metadata.entrySet()) {
            if (entry.getKey() != null && entry.getKey().endsWith(fallbackSuffix) && !isBlank(entry.getValue())) {
                return entry.getValue();
            }
        }
        return "";
    }

    private static boolean isBlank(String value) {
        return value == null || value.isBlank();
    }

    private static final class IntStats {
        private int count;
        private long sum;
        private int min = Integer.MAX_VALUE;
        private int max = Integer.MIN_VALUE;

        private void add(int value) {
            if (value < 0) {
                return;
            }
            count++;
            sum += value;
            min = Math.min(min, value);
            max = Math.max(max, value);
        }

        private int count() {
            return count;
        }

        private int min() {
            return count == 0 ? -1 : min;
        }

        private int max() {
            return count == 0 ? -1 : max;
        }

        private double mean() {
            return count == 0 ? Double.NaN : (double) sum / count;
        }
    }

    private static final class ProductStats {
        private final IntStats atomStats = new IntStats();
        private final IntStats rotatableStats = new IntStats();

        private void add(int atoms, int rotatable) {
            atomStats.add(atoms);
            rotatableStats.add(rotatable);
        }

        private int count() {
            return atomStats.count();
        }
    }

    private record SamplingResult(ProductStats productStats,
                                  int failures,
                                  List<RawSynthonSpaceStatsReport.ExampleProduct> examples) {
    }

    private record ReactionReport(RawSynthonSpaceStatsReport.ReactionStats stats,
                                  List<RawSynthonSpaceStatsReport.SynthonSetStats> synthonSetStats,
                                  List<RawSynthonSpaceStatsReport.ExampleProduct> examples,
                                  Set<String> uniqueIdcodes) {
    }

    private static final class SourceAccumulator {
        private final String source;
        private int reactionCount;
        private int synthonSetCount;
        private long totalSynthons;
        private BigInteger totalProductCount = BigInteger.ZERO;
        private int twoSetReactionCount;
        private int threeSetReactionCount;
        private int otherReactionCount;
        private final Set<String> uniqueIdcodes = new HashSet<>();

        private SourceAccumulator(String source) {
            this.source = source;
        }

        private RawSynthonSpaceStatsReport.SourceStats toStats() {
            return new RawSynthonSpaceStatsReport.SourceStats(
                    source,
                    reactionCount,
                    synthonSetCount,
                    totalSynthons,
                    uniqueIdcodes.size(),
                    totalProductCount,
                    twoSetReactionCount,
                    threeSetReactionCount,
                    otherReactionCount);
        }
    }

    private static final class StatsThreadFactory implements ThreadFactory {
        private final AtomicInteger counter = new AtomicInteger();

        @Override
        public Thread newThread(Runnable runnable) {
            Thread thread = new Thread(runnable, String.format(Locale.ROOT, "RawspaceStats-%d", counter.incrementAndGet()));
            thread.setDaemon(true);
            return thread;
        }
    }
}
