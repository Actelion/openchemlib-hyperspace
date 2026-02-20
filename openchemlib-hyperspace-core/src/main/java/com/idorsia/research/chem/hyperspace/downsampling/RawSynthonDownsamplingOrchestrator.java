package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSet;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Applies a {@link RawSynthonDownsampler} to every synthon set in a {@link RawSynthonSpace}.
 */
public class RawSynthonDownsamplingOrchestrator {

    public SynthonDownsamplingResult downsample(RawSynthonSpace rawSpace,
                                                RawSynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        List<SynthonSetDownsamplingResult> perSetResults = new ArrayList<>();
        for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : rawSpace.getReactions().entrySet()) {
            entry.getValue().getRawFragmentSets().forEach((fragIdx, synthons) -> {
                if (synthons == null || synthons.isEmpty()) {
                    return;
                }
                RawSynthonSet set = new RawSynthonSet(entry.getKey(), fragIdx, synthons);
                SynthonSetDownsamplingResult result = downsampler.downsample(set, request);
                perSetResults.add(result);
            });
        }
        return new SynthonDownsamplingResult(downsampler.getName(), request, perSetResults);
    }
}
