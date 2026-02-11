package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Side-car structure that mirrors a SynthonSpace with downsampled sets.
 */
public class DownsampledSynthonSpace implements Serializable {

    private final String algorithmName;
    private final SynthonDownsamplingRequest request;
    private final Map<String, Map<Integer, List<SynthonSpace.FragId>>> downsampledSets;

    public DownsampledSynthonSpace(String algorithmName,
                                   SynthonDownsamplingRequest request,
                                   Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> representatives) {
        this.algorithmName = Objects.requireNonNull(algorithmName, "algorithmName");
        this.request = Objects.requireNonNull(request, "request");
        this.downsampledSets = buildNestedMap(Objects.requireNonNull(representatives, "representatives"));
    }

    private static Map<String, Map<Integer, List<SynthonSpace.FragId>>> buildNestedMap(Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> representatives) {
        Map<String, Map<Integer, List<SynthonSpace.FragId>>> nested = new HashMap<>();
        for (Map.Entry<SynthonSpace.FragType, List<SynthonSpace.FragId>> entry : representatives.entrySet()) {
            SynthonSpace.FragType fragType = entry.getKey();
            nested.computeIfAbsent(fragType.rxn_id, key -> new HashMap<>())
                    .put(fragType.frag, Collections.unmodifiableList(new ArrayList<>(entry.getValue())));
        }
        Map<String, Map<Integer, List<SynthonSpace.FragId>>> immutableOuter = new HashMap<>();
        for (Map.Entry<String, Map<Integer, List<SynthonSpace.FragId>>> entry : nested.entrySet()) {
            immutableOuter.put(entry.getKey(), Collections.unmodifiableMap(entry.getValue()));
        }
        return Collections.unmodifiableMap(immutableOuter);
    }

    public String getAlgorithmName() {
        return algorithmName;
    }

    public SynthonDownsamplingRequest getRequest() {
        return request;
    }

    public Map<String, Map<Integer, List<SynthonSpace.FragId>>> getDownsampledSets() {
        return downsampledSets;
    }

    public List<SynthonSpace.FragId> getSynthonSet(String rxnId, int fragIdx) {
        Map<Integer, List<SynthonSpace.FragId>> rxn = downsampledSets.get(rxnId);
        if (rxn == null) {
            return Collections.emptyList();
        }
        return rxn.getOrDefault(fragIdx, Collections.emptyList());
    }
}
