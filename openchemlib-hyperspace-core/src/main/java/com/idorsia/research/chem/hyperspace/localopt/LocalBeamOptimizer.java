package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

public class LocalBeamOptimizer {

    private final SynthonSetAccessor provider;
    private final AssemblyScorer scorer;
    private final NeighborSampler neighborSampler;

    public LocalBeamOptimizer(SynthonSetAccessor provider, AssemblyScorer scorer) {
        this(provider, scorer, new SkelSpheresNeighborSampler());
    }

    public LocalBeamOptimizer(SynthonSetAccessor provider,
                              AssemblyScorer scorer,
                              NeighborSampler neighborSampler) {
        this.provider = provider;
        this.scorer = scorer;
        this.neighborSampler = neighborSampler;
    }

    public LocalOptimizationResult optimize(SeedAssembly seed, LocalOptimizationRequest request) {
        Map<Integer, List<SynthonSpace.FragId> > synthonSets = provider.getSynthonSets(seed.getReactionId());
        if (synthonSets.isEmpty()) {
            throw new IllegalArgumentException("No synthon sets for reaction " + seed.getReactionId());
        }
        List<Integer> positions = new ArrayList<>(synthonSets.keySet());
        Collections.sort(positions);
        List<SynthonSpace.FragId> seedFragments = mapSeedFragments(seed, positions);
        LocalOptimizationLogger logger = new LocalOptimizationLogger(request.getLogLevel(), seed.getFragmentIds());

        Map<String, LocalOptimizationResult.BeamEntry> scoreCache = new HashMap<>();
        LocalOptimizationResult.BeamEntry seedEntry = scoreCandidate(seed.getReactionId(), seedFragments, scoreCache, 0, logger);
        if (seedEntry == null) {
            throw new IllegalStateException("Unable to score seed assembly for reaction " + seed.getReactionId());
        }
        List<LocalOptimizationResult.BeamEntry> beam = new ArrayList<>();
        beam.add(seedEntry);
        double bestScore = seedEntry.getScore();
        int noImprovementRounds = 0;
        Random rng = new Random(request.getRandomSeed() ^ seedKey(seed));

        for (int round = 1; round <= request.getMaxRounds(); round++) {
            for (int posIndex = 0; posIndex < positions.size(); posIndex++) {
                int position = positions.get(posIndex);
                beam = expandPosition(seed.getReactionId(), position, posIndex, beam, synthonSets,
                        request, scoreCache, rng, logger);
                if (beam.isEmpty()) {
                    break;
                }
            }
            if (beam.isEmpty()) {
                break;
            }
            double roundBest = beam.get(0).getScore();
            if (roundBest > bestScore + request.getImprovementTolerance()) {
                bestScore = roundBest;
                noImprovementRounds = 0;
                logger.logImprovement(round, roundBest, beam.get(0).getIdcode());
            } else {
                noImprovementRounds++;
                if (noImprovementRounds >= request.getPatience()) {
                    break;
                }
            }
        }

        return new LocalOptimizationResult(seed.getReactionId(), beam, seed.getFragmentIds());
    }

    private List<LocalOptimizationResult.BeamEntry> expandPosition(String reactionId,
                                                                   int positionKey,
                                                                   int positionIndex,
                                                                   List<LocalOptimizationResult.BeamEntry> beam,
                                                                   Map<Integer, List<SynthonSpace.FragId> > synthonSets,
                                                                   LocalOptimizationRequest request,
                                                                   Map<String, LocalOptimizationResult.BeamEntry> scoreCache,
                                                                   Random rng,
                                                                   LocalOptimizationLogger logger) {
        Map<String, LocalOptimizationResult.BeamEntry> candidates = new LinkedHashMap<>();
        for (LocalOptimizationResult.BeamEntry entry : beam) {
            candidates.put(buildKey(reactionId, entry.getFragments()), entry);
            List<SynthonSpace.FragId> sampled = neighborSampler.sampleNeighbors(reactionId,
                    positionKey,
                    entry.getFragments().get(positionIndex),
                    synthonSets,
                    request,
                    rng);
            for (SynthonSpace.FragId neighbor : sampled) {
                if (neighbor.fragment_id.equals(entry.getFragments().get(positionIndex).fragment_id)) {
                    continue;
                }
                List<SynthonSpace.FragId> updated = new ArrayList<>(entry.getFragments());
                updated.set(positionIndex, neighbor);
                String key = buildKey(reactionId, updated);
                if (candidates.containsKey(key)) {
                    continue;
                }
                LocalOptimizationResult.BeamEntry scored = scoreCandidate(reactionId, updated, scoreCache, positionIndex + 1, logger);
                if (scored != null) {
                    candidates.put(key, scored);
                }
            }
        }
        if (candidates.isEmpty()) {
            return Collections.emptyList();
        }
        List<LocalOptimizationResult.BeamEntry> sorted = new ArrayList<>(candidates.values());
        sorted.sort(Comparator.comparingDouble(LocalOptimizationResult.BeamEntry::getScore).reversed());
        return selectBeam(sorted, request.getBeamSize(), request.getPerPositionCap());
    }

    private List<LocalOptimizationResult.BeamEntry> selectBeam(List<LocalOptimizationResult.BeamEntry> candidates,
                                                               int beamSize,
                                                               int perPositionCap) {
        List<LocalOptimizationResult.BeamEntry> result = new ArrayList<>();
        Map<Integer, Map<String, Integer> > usage = new HashMap<Integer, Map<String, Integer> >();
        for (LocalOptimizationResult.BeamEntry candidate : candidates) {
            if (result.size() >= beamSize) {
                break;
            }
            if (!violatesCaps(candidate, usage, perPositionCap)) {
                addCandidate(candidate, usage);
                result.add(candidate);
            }
        }
        if (result.size() < beamSize) {
            for (LocalOptimizationResult.BeamEntry candidate : candidates) {
                if (result.contains(candidate)) {
                    continue;
                }
                result.add(candidate);
                if (result.size() >= beamSize) {
                    break;
                }
            }
        }
        return result;
    }

    private boolean violatesCaps(LocalOptimizationResult.BeamEntry entry,
                                 Map<Integer, Map<String, Integer> > usage,
                                 int cap) {
        if (cap <= 0) {
            return false;
        }
        for (int idx = 0; idx < entry.getFragments().size(); idx++) {
            String fragmentId = entry.getFragments().get(idx).fragment_id;
            Map<String, Integer> counts = usage.computeIfAbsent(idx, k -> new HashMap<>());
            int current = counts.getOrDefault(fragmentId, 0);
            if (current >= cap) {
                return true;
            }
        }
        return false;
    }

    private void addCandidate(LocalOptimizationResult.BeamEntry entry,
                              Map<Integer, Map<String, Integer> > usage) {
        for (int idx = 0; idx < entry.getFragments().size(); idx++) {
            Map<String, Integer> counts = usage.computeIfAbsent(idx, k -> new HashMap<>());
            String fragmentId = entry.getFragments().get(idx).fragment_id;
            counts.put(fragmentId, counts.getOrDefault(fragmentId, 0) + 1);
        }
    }

    private List<SynthonSpace.FragId> mapSeedFragments(SeedAssembly seed, List<Integer> positions) {
        if (seed.getFragmentIds().size() != positions.size()) {
            throw new IllegalArgumentException("Seed fragment count does not match reaction definition for " + seed.getReactionId());
        }
        List<SynthonSpace.FragId> fragments = new ArrayList<>();
        for (int i = 0; i < positions.size(); i++) {
            String fragmentId = seed.getFragmentIds().get(i);
            SynthonSpace.FragId frag = provider.findFragment(seed.getReactionId(), fragmentId);
            if (frag == null) {
                throw new IllegalArgumentException("Fragment " + fragmentId + " not found for reaction " + seed.getReactionId());
            }
            fragments.add(frag);
        }
        return fragments;
    }

    private LocalOptimizationResult.BeamEntry scoreCandidate(String reactionId,
                                                             List<SynthonSpace.FragId> fragments,
                                                             Map<String, LocalOptimizationResult.BeamEntry> cache,
                                                             int originatingRound,
                                                             LocalOptimizationLogger logger) {
        String key = buildKey(reactionId, fragments);
        LocalOptimizationResult.BeamEntry cached = cache.get(key);
        if (cached != null) {
            return cached;
        }
        LocalOptimizationResult.BeamEntry entry = scorer.score(reactionId, fragments, originatingRound, logger);
        if (entry != null) {
            cache.put(key, entry);
        }
        return entry;
    }

    private String buildKey(String reactionId, List<SynthonSpace.FragId> fragments) {
        return reactionId + ":" + fragments.stream().map(f -> f.fragment_id).collect(Collectors.joining("|"));
    }

    private long seedKey(SeedAssembly seed) {
        return (seed.getReactionId() + ":" + String.join("|", seed.getFragmentIds())).hashCode();
    }
}
