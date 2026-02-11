package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

public class SynthonDownsamplingResult implements Serializable {

    private final SynthonSpace sourceSpace;
    private final String algorithmName;
    private final SynthonDownsamplingRequest request;
    private final List<SynthonSetDownsamplingResult> setResults;

    public SynthonDownsamplingResult(SynthonSpace sourceSpace,
                                     String algorithmName,
                                     SynthonDownsamplingRequest request,
                                     List<SynthonSetDownsamplingResult> setResults) {
        this.sourceSpace = Objects.requireNonNull(sourceSpace, "sourceSpace");
        this.algorithmName = Objects.requireNonNull(algorithmName, "algorithmName");
        this.request = Objects.requireNonNull(request, "request");
        this.setResults = Collections.unmodifiableList(new ArrayList<>(setResults));
    }

    public SynthonSpace getSourceSpace() {
        return sourceSpace;
    }

    public String getAlgorithmName() {
        return algorithmName;
    }

    public SynthonDownsamplingRequest getRequest() {
        return request;
    }

    public List<SynthonSetDownsamplingResult> getSetResults() {
        return setResults;
    }

    public long getTotalOriginalSynthons() {
        return setResults.stream().mapToLong(SynthonSetDownsamplingResult::getOriginalSize).sum();
    }

    public long getTotalRetainedSynthons() {
        return setResults.stream().mapToLong(SynthonSetDownsamplingResult::getRetainedSize).sum();
    }

    public long getTotalDiscardedSynthons() {
        return Math.max(0L, getTotalOriginalSynthons() - getTotalRetainedSynthons());
    }

    public DownsampledSynthonSpace toDownsampledSpace() {
        Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> map = new HashMap<>();
        for (SynthonSetDownsamplingResult result : setResults) {
            map.put(result.getFragType(), result.getRepresentatives());
        }
        return new DownsampledSynthonSpace(algorithmName, request, map);
    }
}
