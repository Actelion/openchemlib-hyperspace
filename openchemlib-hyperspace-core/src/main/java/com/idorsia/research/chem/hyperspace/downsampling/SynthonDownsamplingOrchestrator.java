package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Traverses a SynthonSpace and applies a SynthonDownsampler set-by-set.
 */
public class SynthonDownsamplingOrchestrator {

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        List<SynthonSpace.FragType> allTypes = new ArrayList<>();
        for (String rxnId : space.getRxnIds()) {
            Map<Integer, SynthonSpace.FragType> ft = space.getFragTypes(rxnId);
            if (ft != null) {
                allTypes.addAll(ft.values());
            }
        }
        return downsample(space, allTypes, downsampler, request);
    }

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                Collection<SynthonSpace.FragType> fragTypes,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        List<SynthonSetDownsamplingResult> perSetResults = new ArrayList<>();
        for (SynthonSpace.FragType fragType : fragTypes) {
            List<SynthonSpace.FragId> synthons = new ArrayList<>(space.getSynthonSet(fragType.rxn_id, fragType.frag));
            if (synthons.isEmpty()) {
                continue;
            }
            SynthonSetDownsamplingResult result = downsampler.downsample(space, fragType, synthons, request);
            perSetResults.add(result);
        }
        return new SynthonDownsamplingResult(space, downsampler.getName(), request, perSetResults);
    }
}
