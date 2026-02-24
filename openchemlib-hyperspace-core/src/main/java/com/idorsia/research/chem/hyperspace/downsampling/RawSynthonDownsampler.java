package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSet;

/**
 * Downsampler operating directly on raw synthon sets.
 */
public interface RawSynthonDownsampler {

    String getName();

    SynthonSetDownsamplingResult downsample(RawSynthonSet synthonSet,
                                            SynthonDownsamplingRequest request);
}
