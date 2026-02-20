package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Strategy interface to select neighboring synthons for beam expansions.
 */
public interface NeighborSampler {

    List<SynthonSpace.FragId> sampleNeighbors(String reactionId,
                                              int fragIdx,
                                              SynthonSpace.FragId center,
                                              Map<Integer, List<SynthonSpace.FragId>> synthonSets,
                                              LocalOptimizationRequest request,
                                              Random rng);
}
