package com.idorsia.research.chem.hyperspace.localopt;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDescriptorCache;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

class NeighborCache {

    private final SynthonDescriptorCache<int[]> descriptorCache;
    private final ThreadLocal<DescriptorHandlerBinarySkelSpheres> descriptorHandlers;
    private final Map<NeighborKey, List<SynthonSpace.FragId> > cache = new ConcurrentHashMap<>();
    private final IDCodeParser parser = new IDCodeParser();

    NeighborCache() {
        this.descriptorHandlers = ThreadLocal.withInitial(DescriptorHandlerBinarySkelSpheres::new);
        this.descriptorCache = new SynthonDescriptorCache<>(this::computeDescriptor);
    }

    @SuppressWarnings("unchecked")
    List<SynthonSpace.FragId> getNeighbors(String reactionId,
                                           int fragIdx,
                                           SynthonSpace.FragId center,
                                           Map pool,
                                           int topL) {
        Map<Integer, List<SynthonSpace.FragId> > typedPool = (Map<Integer, List<SynthonSpace.FragId> >) pool;
        NeighborKey key = new NeighborKey(reactionId, fragIdx, center.fragment_id);
        return cache.computeIfAbsent(key, k -> computeNeighbors(center, typedPool.get(fragIdx), topL));
    }

    private List<SynthonSpace.FragId> computeNeighbors(SynthonSpace.FragId center,
                                                       List<SynthonSpace.FragId> candidates,
                                                       int topL) {
        if (candidates == null || candidates.isEmpty()) {
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
        if (sorted.size() > topL) {
            return Collections.unmodifiableList(new ArrayList<>(sorted.subList(0, topL)));
        }
        return Collections.unmodifiableList(sorted);
    }

    private int[] computeDescriptor(SynthonSpace.FragId fragId) {
        try {
            StereoMolecule molecule = new StereoMolecule();
            parser.parse(molecule, fragId.idcode);
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

        private NeighborKey(String reactionId, int fragIdx, String fragmentId) {
            this.reactionId = reactionId;
            this.fragIdx = fragIdx;
            this.fragmentId = fragmentId;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof NeighborKey)) return false;
            NeighborKey that = (NeighborKey) o;
            return fragIdx == that.fragIdx &&
                    reactionId.equals(that.reactionId) &&
                    fragmentId.equals(that.fragmentId);
        }

        @Override
        public int hashCode() {
            int result = reactionId.hashCode();
            result = 31 * result + fragIdx;
            result = 31 * result + fragmentId.hashCode();
            return result;
        }
    }
}
