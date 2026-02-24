package com.idorsia.research.chem.hyperspace.localopt;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDescriptorCache;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

/**
 * NeighborSampler implementation that ranks synthons by SkelSpheres similarity per fragment.
 */
public class SkelSpheresNeighborSampler implements NeighborSampler {

    private final SynthonDescriptorCache<int[]> descriptorCache;
    private final ThreadLocal<DescriptorHandlerBinarySkelSpheres> descriptorHandlers;
    private final Map<NeighborKey, List<SynthonSpace.FragId>> cache = new ConcurrentHashMap<>();

    public SkelSpheresNeighborSampler() {
        this.descriptorHandlers = ThreadLocal.withInitial(DescriptorHandlerBinarySkelSpheres::new);
        this.descriptorCache = new SynthonDescriptorCache<>(this::computeDescriptor);
    }

    @Override
    public List<SynthonSpace.FragId> sampleNeighbors(String reactionId,
                                                     int fragIdx,
                                                     SynthonSpace.FragId center,
                                                     Map<Integer, List<SynthonSpace.FragId>> synthonSets,
                                                     LocalOptimizationRequest request,
                                                     Random rng) {
        if (request.getSampledNeighbors() <= 0) {
            return Collections.emptyList();
        }
        List<SynthonSpace.FragId> pool = synthonSets.get(fragIdx);
        if (pool == null || pool.isEmpty()) {
            return Collections.emptyList();
        }
        int neighborPoolSize = request.getNeighborPoolSize();
        NeighborKey key = new NeighborKey(reactionId, fragIdx, center.fragment_id, neighborPoolSize);
        List<SynthonSpace.FragId> topNeighbors = cache.computeIfAbsent(key,
                ignored -> computeNeighbors(center, pool, neighborPoolSize));
        if (topNeighbors.isEmpty()) {
            return topNeighbors;
        }
        int sampleSize = Math.min(request.getSampledNeighbors(), topNeighbors.size());
        List<SynthonSpace.FragId> copy = new ArrayList<>(topNeighbors);
        Collections.shuffle(copy, rng);
        return Collections.unmodifiableList(new ArrayList<>(copy.subList(0, sampleSize)));
    }

    private List<SynthonSpace.FragId> computeNeighbors(SynthonSpace.FragId center,
                                                       List<SynthonSpace.FragId> candidates,
                                                       int topL) {
        if (candidates.isEmpty()) {
            return Collections.emptyList();
        }
        int[] centerDesc = descriptorCache.getOrCompute(center);
        if (centerDesc == null) {
            return Collections.emptyList();
        }
        DescriptorHandlerBinarySkelSpheres handler = descriptorHandlers.get();
        List<SynthonSpace.FragId> sorted = new ArrayList<>(candidates);
        sorted.sort(Comparator.comparingDouble((SynthonSpace.FragId f) ->
                handler.getSimilarity(centerDesc, descriptorCache.getOrCompute(f))).reversed());
        if (topL > 0 && sorted.size() > topL) {
            sorted = new ArrayList<>(sorted.subList(0, topL));
        }
        return Collections.unmodifiableList(sorted);
    }

    private int[] computeDescriptor(SynthonSpace.FragId fragId) {
        try {
            StereoMolecule molecule = new StereoMolecule();
            SynchronizedIDCodeParser.parse(molecule, fragId.idcode);
            molecule.ensureHelperArrays(StereoMolecule.cHelperCIP);
            return descriptorHandlers.get().createDescriptor(molecule);
        } catch (Exception e) {
            return null;
        }
    }

    private static final class NeighborKey {
        private final String reactionId;
        private final int fragIdx;
        private final String fragmentId;
        private final int poolSize;

        private NeighborKey(String reactionId, int fragIdx, String fragmentId, int poolSize) {
            this.reactionId = reactionId;
            this.fragIdx = fragIdx;
            this.fragmentId = fragmentId;
            this.poolSize = poolSize;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof NeighborKey)) return false;
            NeighborKey that = (NeighborKey) o;
            return fragIdx == that.fragIdx &&
                    poolSize == that.poolSize &&
                    reactionId.equals(that.reactionId) &&
                    fragmentId.equals(that.fragmentId);
        }

        @Override
        public int hashCode() {
            int result = reactionId.hashCode();
            result = 31 * result + fragIdx;
            result = 31 * result + fragmentId.hashCode();
            result = 31 * result + poolSize;
            return result;
        }
    }
}
