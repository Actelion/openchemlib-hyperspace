package com.idorsia.research.chem.hyperspace.screening;

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;

/**
 * Produces random candidates from a downsampled synthon pool.
 */
public final class CandidateSampler {

    private final Map<String, ReactionSampler> samplers = new ConcurrentHashMap<>();
    private final PheSAMolecule query;
    private final Config config;
    private final LongAdder comparisonCounter;
    private final ReactionSamplingTracker samplingTracker;

    public CandidateSampler(DownsampledSynthonSpace downsampledSpace,
                            PheSAMolecule query,
                            Config config) {
        this(downsampledSpace, query, config, null, null);
    }

    public CandidateSampler(DownsampledSynthonSpace downsampledSpace,
                            PheSAMolecule query,
                            Config config,
                            LongAdder comparisonCounter) {
        this(downsampledSpace, query, config, comparisonCounter, null);
    }

    public CandidateSampler(DownsampledSynthonSpace downsampledSpace,
                            PheSAMolecule query,
                            Config config,
                            LongAdder comparisonCounter,
                            ReactionSamplingTracker samplingTracker) {
        this.query = query;
        this.config = config;
        this.comparisonCounter = comparisonCounter;
        this.samplingTracker = samplingTracker;
        downsampledSpace.getDownsampledSets().forEach((rxnId, sets) ->
                samplers.put(rxnId, new ReactionSampler(rxnId, sets, config, comparisonCounter, samplingTracker)));
    }

    public ScreeningCandidate sample(String reactionId, Random rng) {
        ReactionSampler sampler = samplers.get(reactionId);
        if (sampler == null) {
            return null;
        }
        return sampler.sample(rng, query);
    }

    public static final class Config {
        private final long attemptsPerReaction;
        private final int minAtoms;
        private final int maxAtoms;
        private final int maxRotatableBonds;
        private final double minSimilarity;

        public Config(long attemptsPerReaction,
                      int minAtoms,
                      int maxAtoms,
                      int maxRotatableBonds,
                      double minSimilarity) {
            this.attemptsPerReaction = Math.max(1, attemptsPerReaction);
            this.minAtoms = Math.max(1, minAtoms);
            this.maxAtoms = maxAtoms;
            this.maxRotatableBonds = Math.max(1, maxRotatableBonds);
            this.minSimilarity = minSimilarity;
        }
    }

    private static final class ReactionSampler {
        private final String reactionId;
        private final Map<Integer, List<SynthonSpace.FragId>> synthonSets;
        private final List<Integer> fragOrder;
        private final ThreadLocal<IDCodeParser> parser = ThreadLocal.withInitial(IDCodeParser::new);
        private final ThreadLocal<DescriptorHandlerShape> descriptorHandler =
                ThreadLocal.withInitial(DescriptorHandlerShape::new);
        private final ThreadLocal<ConformerSetGenerator> conformerGenerator =
                ThreadLocal.withInitial(() -> new ConformerSetGenerator(1));
        private final Config config;
        private final LongAdder comparisonCounter;
        private final ReactionSamplingTracker samplingTracker;

        private ReactionSampler(String reactionId,
                                Map<Integer, List<SynthonSpace.FragId>> synthonSets,
                                Config config,
                                LongAdder comparisonCounter,
                                ReactionSamplingTracker samplingTracker) {
            this.reactionId = reactionId;
            this.synthonSets = synthonSets;
            this.config = config;
            this.comparisonCounter = comparisonCounter;
            this.samplingTracker = samplingTracker;
            List<Integer> order = new ArrayList<>(synthonSets.keySet());
            order.sort(Comparator.naturalOrder());
            this.fragOrder = Collections.unmodifiableList(order);
        }

        private ScreeningCandidate sample(Random rng, PheSAMolecule query) {
            if (fragOrder.size() < 2) {
                return null;
            }
            List<StereoMolecule> parts = new ArrayList<>(fragOrder.size());
            List<SynthonSpace.FragId> chosen = new ArrayList<>(fragOrder.size());
            for (long attempt = 0; attempt < config.attemptsPerReaction; attempt++) {
                if (samplingTracker != null) {
                    samplingTracker.recordAttempt(reactionId);
                }
                parts.clear();
                chosen.clear();
                boolean missing = false;
                for (int fragIdx : fragOrder) {
                    List<SynthonSpace.FragId> pool = synthonSets.get(fragIdx);
                    if (pool == null || pool.isEmpty()) {
                        missing = true;
                        break;
                    }
                    SynthonSpace.FragId frag = pool.get(rng.nextInt(pool.size()));
                    chosen.add(frag);
                    parts.add(parse(frag));
                }
                if (missing) {
                    return null;
                }

                StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
                int atoms = assembled.getAtoms();
                if (atoms < config.minAtoms) {
                    continue;
                }
                if (config.maxAtoms > 0 && atoms > config.maxAtoms) {
                    continue;
                }
                assembled.ensureHelperArrays(Molecule.cHelperCIP);
                int rotatable = countRotatableBonds(assembled);
                if (rotatable > config.maxRotatableBonds) {
                    continue;
                }

                ConformerSetGenerator generator = conformerGenerator.get();
                ConformerSet conformers = generator.generateConformerSet(assembled);
                if (conformers.isEmpty()) {
                    continue;
                }
                DescriptorHandlerShape handler = descriptorHandler.get();
                PheSAMolecule descriptor = handler.createDescriptor(conformers);
                if (handler.calculationFailed(descriptor)) {
                    continue;
                }
                double similarity = handler.getSimilarity(query, descriptor);
                if (comparisonCounter != null) {
                    comparisonCounter.increment();
                }
                if (similarity < config.minSimilarity) {
                    continue;
                }
                List<String> fragmentIds = new ArrayList<>(chosen.size());
                for (SynthonSpace.FragId fragId : chosen) {
                    fragmentIds.add(fragId.fragment_id);
                }
                if (samplingTracker != null) {
                    samplingTracker.recordSuccess(reactionId);
                }
                return new ScreeningCandidate(reactionId,
                        fragmentIds,
                        assembled.getIDCode(),
                        similarity,
                        atoms,
                        rotatable);
            }
            return null;
        }

        private StereoMolecule parse(SynthonSpace.FragId frag) {
            StereoMolecule molecule = new StereoMolecule();
            parser.get().parse(molecule, frag.idcode);
            molecule.ensureHelperArrays(Molecule.cHelperCIP);
            return molecule;
        }

        private static int countRotatableBonds(StereoMolecule molecule) {
            molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
            boolean[] rotatable = new boolean[molecule.getBonds()];
            TorsionDB.findRotatableBonds(molecule, true, rotatable);
            int count = 0;
            for (boolean r : rotatable) {
                if (r) {
                    count++;
                }
            }
            return count;
        }
    }

    public interface ReactionSamplingTracker {
        void recordAttempt(String reactionId);

        void recordSuccess(String reactionId);
    }
}
