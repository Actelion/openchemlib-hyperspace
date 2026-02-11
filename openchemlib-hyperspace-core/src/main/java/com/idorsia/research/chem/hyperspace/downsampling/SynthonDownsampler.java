package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.List;

/**
 * Strategy interface that transforms a synthon set into a representative subset.
 */
public interface SynthonDownsampler {

    String getName();

    SynthonSetDownsamplingResult downsample(SynthonSpace sourceSpace,
                                            SynthonSpace.FragType fragType,
                                            List<SynthonSpace.FragId> synthons,
                                            SynthonDownsamplingRequest request);
}
