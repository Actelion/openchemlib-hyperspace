package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class SynthonSetAccessor {

    private final DownsampledSynthonSpace downsampledSpace;
    private final SynthonSpace synthonSpace;
    private final Map<String, Map<Integer, List<SynthonSpace.FragId> > > cache = new ConcurrentHashMap<>();
    private final Map<String, Map<String, SynthonSpace.FragId> > fragmentLookup = new ConcurrentHashMap<>();

    public SynthonSetAccessor(DownsampledSynthonSpace downsampledSpace) {
        this.downsampledSpace = downsampledSpace;
        this.synthonSpace = null;
    }

    public SynthonSetAccessor(SynthonSpace synthonSpace) {
        this.downsampledSpace = null;
        this.synthonSpace = synthonSpace;
    }

    @SuppressWarnings("unchecked")
    public Map<Integer, List<SynthonSpace.FragId> > getSynthonSets(String reactionId) {
        if (downsampledSpace != null) {
            Map<Integer, List<SynthonSpace.FragId> > map = downsampledSpace.getDownsampledSets().get(reactionId);
            return map == null ? Collections.emptyMap() : map;
        }
        return cache.computeIfAbsent(reactionId, this::loadFromSpace);
    }

    public SynthonSpace.FragId findFragment(String reactionId, String fragmentId) {
        Map<String, SynthonSpace.FragId> lookup = fragmentLookup.computeIfAbsent(reactionId, this::buildLookup);
        return lookup.get(fragmentId);
    }

    private Map<Integer, List<SynthonSpace.FragId> > loadFromSpace(String reactionId) {
        if (synthonSpace == null) {
            return Collections.emptyMap();
        }
        Map<Integer, SynthonSpace.FragType> fragTypes = synthonSpace.getFragTypes(reactionId);
        if (fragTypes == null) {
            return Collections.emptyMap();
        }
        Map<Integer, List<SynthonSpace.FragId> > result = new HashMap<>();
        for (int fragIndex : fragTypes.keySet()) {
            List<SynthonSpace.FragId> list = synthonSpace.getSynthonSet(reactionId, fragIndex);
            if (list != null) {
                result.put(fragIndex, Collections.unmodifiableList(new ArrayList<>(list)));
            }
        }
        return result;
    }

    private Map<String, SynthonSpace.FragId> buildLookup(String reactionId) {
        Map<String, SynthonSpace.FragId> map = new HashMap<>();
        Map<Integer, List<SynthonSpace.FragId> > sets = getSynthonSets(reactionId);
        for (List<SynthonSpace.FragId> list : sets.values()) {
            for (SynthonSpace.FragId frag : list) {
                map.put(frag.fragment_id, frag);
            }
        }
        return map;
    }
}
