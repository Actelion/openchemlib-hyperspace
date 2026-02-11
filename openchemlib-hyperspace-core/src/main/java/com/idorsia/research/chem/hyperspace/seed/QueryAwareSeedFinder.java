package com.idorsia.research.chem.hyperspace.seed;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.seed.QueryAwareSeedFinderResult.ReactionStats;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Performs the random query-aware seed discovery on downsampled synthon pools.
 */
public class QueryAwareSeedFinder {

    private final DownsampledSynthonSpace downsampledSpace;

    public QueryAwareSeedFinder(DownsampledSynthonSpace downsampledSpace) {
        this.downsampledSpace = downsampledSpace;
    }

    public QueryAwareSeedFinderResult run(PheSAMolecule query,
                                          QueryAwareSeedFinderRequest request) throws IOException {
        List<ReactionStats> stats = new ArrayList<>();
        Map<String, Map<Integer, List<SynthonSpace.FragId>>> rxnMap = downsampledSpace.getDownsampledSets();
        if (rxnMap.isEmpty()) {
            return new QueryAwareSeedFinderResult(Collections.emptyList());
        }

        ThreadLocal<IDCodeParser> parser = ThreadLocal.withInitial(IDCodeParser::new);
        ThreadLocal<DescriptorHandlerShape> descriptorHandlers = ThreadLocal.withInitial(DescriptorHandlerShape::new);
        ThreadLocal<ConformerSetGenerator> conformerGenerators = ThreadLocal.withInitial(() -> new ConformerSetGenerator(1));
        Path outputPath = request.getOutputPath();

        try (SeedHitWriter writer = new SeedHitWriter(outputPath, request.getMaxHitsPerFile())) {
            ExecutorService executor = Executors.newFixedThreadPool(request.getNumThreads(), new SeedFinderThreadFactory());
            try {
                List<Future<ReactionStats>> futures = new ArrayList<>();
                for (Map.Entry<String, Map<Integer, List<SynthonSpace.FragId>>> entry : rxnMap.entrySet()) {
                    futures.add(executor.submit(new ReactionTask(entry.getKey(),
                            entry.getValue(),
                            query,
                            request,
                            writer,
                            parser,
                            descriptorHandlers,
                            conformerGenerators)));
                }
                for (Future<ReactionStats> future : futures) {
                    stats.add(future.get());
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IOException("Seed finder interrupted", e);
            } catch (ExecutionException e) {
                throw new IOException("Seed finder failed", e.getCause());
            } finally {
                executor.shutdownNow();
            }
        }
        return new QueryAwareSeedFinderResult(stats);
    }

    private static final class ReactionTask implements Callable<ReactionStats> {
        private final String reactionId;
        private final Map<Integer, List<SynthonSpace.FragId>> synthonSets;
        private final PheSAMolecule query;
        private final QueryAwareSeedFinderRequest request;
        private final SeedHitWriter writer;
        private final ThreadLocal<IDCodeParser> parserSupplier;
        private final ThreadLocal<DescriptorHandlerShape> descriptorSupplier;
        private final ThreadLocal<ConformerSetGenerator> conformerGeneratorSupplier;

        private ReactionTask(String reactionId,
                             Map<Integer, List<SynthonSpace.FragId>> synthonSets,
                             PheSAMolecule query,
                             QueryAwareSeedFinderRequest request,
                             SeedHitWriter writer,
                             ThreadLocal<IDCodeParser> parserSupplier,
                             ThreadLocal<DescriptorHandlerShape> descriptorSupplier,
                             ThreadLocal<ConformerSetGenerator> conformerGeneratorSupplier) {
            this.reactionId = reactionId;
            this.synthonSets = synthonSets;
            this.query = query;
            this.request = request;
            this.writer = writer;
            this.parserSupplier = parserSupplier;
            this.descriptorSupplier = descriptorSupplier;
            this.conformerGeneratorSupplier = conformerGeneratorSupplier;
        }

        @Override
        public ReactionStats call() {
            ReactionStats stats = new ReactionStats(reactionId);
            if (synthonSets.size() < 2) {
                return stats;
            }

            List<Integer> fragOrder = new ArrayList<>(synthonSets.keySet());
            fragOrder.sort(Comparator.naturalOrder());
            Random rng = new Random(request.getRandomSeed() ^ reactionId.hashCode());
            List<StereoMolecule> parts = new ArrayList<>(fragOrder.size());
            List<SynthonSpace.FragId> chosen = new ArrayList<>(fragOrder.size());

            for (long attempt = 0; attempt < request.getAttemptsPerReaction(); attempt++) {
                stats.incrementAttempt();
                parts.clear();
                chosen.clear();

                boolean missingSynthons = false;
                for (int fragIdx : fragOrder) {
                    List<SynthonSpace.FragId> pool = synthonSets.get(fragIdx);
                    if (pool == null || pool.isEmpty()) {
                        missingSynthons = true;
                        break;
                    }
                    SynthonSpace.FragId frag = pool.get(rng.nextInt(pool.size()));
                    chosen.add(frag);
                    parts.add(parseFrag(frag));
                }
                if (missingSynthons) {
                    break;
                }

                StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
                stats.incrementAssemblies();
                int atomCount = assembled.getAtoms();
                if (atomCount < request.getMinAtoms()) {
                    continue;
                }
                if (request.getMaxAtoms() > 0 && atomCount > request.getMaxAtoms()) {
                    continue;
                }

                assembled.ensureHelperArrays(Molecule.cHelperCIP);
                int rotatable = countRotatableBonds(assembled);
                if (rotatable > request.getMaxRotatableBonds()) {
                    continue;
                }
                stats.incrementThresholdPass();

                ConformerSetGenerator generator = conformerGeneratorSupplier.get();
                ConformerSet conformers = generator.generateConformerSet(assembled);
                if (conformers.isEmpty()) {
                    continue;
                }
                DescriptorHandlerShape descriptorHandler = descriptorSupplier.get();
                PheSAMolecule candidate = descriptorHandler.createDescriptor(conformers);
                if (descriptorHandler.calculationFailed(candidate)) {
                    continue;
                }
                stats.incrementPhesaEvaluation();
                double similarity = descriptorHandler.getSimilarity(query, candidate);
                if (similarity < request.getMinPhesaSimilarity()) {
                    continue;
                }
                stats.incrementHit();
                writer.writeHit(new SeedHitWriter.HitRecord(
                        reactionId,
                        extractFragmentIds(chosen),
                        assembled.getIDCode(),
                        atomCount,
                        rotatable,
                        similarity,
                        stats.getAttemptCount()));
            }
            return stats;
        }

        private StereoMolecule parseFrag(SynthonSpace.FragId frag) {
            StereoMolecule molecule = new StereoMolecule();
            parserSupplier.get().parse(molecule, frag.idcode);
            molecule.ensureHelperArrays(Molecule.cHelperCIP);
            return molecule;
        }

        private static List<String> extractFragmentIds(List<SynthonSpace.FragId> fragIds) {
            List<String> ids = new ArrayList<>(fragIds.size());
            for (SynthonSpace.FragId fragId : fragIds) {
                ids.add(fragId.fragment_id);
            }
            return ids;
        }

        private int countRotatableBonds(StereoMolecule molecule) {
            molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
            boolean[] rotatable = new boolean[molecule.getBonds()];
            TorsionDB.findRotatableBonds(molecule, true, rotatable);
            int count = 0;
            for (boolean b : rotatable) {
                if (b) {
                    count++;
                }
            }
            return count;
        }
    }

    private static final class SeedFinderThreadFactory implements ThreadFactory {
        private final AtomicInteger counter = new AtomicInteger();

        @Override
        public Thread newThread(Runnable r) {
            Thread t = new Thread(r, "SeedFinder-" + counter.incrementAndGet());
            t.setDaemon(true);
            return t;
        }
    }
}
