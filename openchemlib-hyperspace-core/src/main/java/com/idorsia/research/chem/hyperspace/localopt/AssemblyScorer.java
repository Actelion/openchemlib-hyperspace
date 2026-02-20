package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.List;

/**
 * Strategy interface that scores assembled synthon combinations.
 */
public interface AssemblyScorer {

    LocalOptimizationResult.BeamEntry score(String reactionId,
                                            List<SynthonSpace.FragId> fragments,
                                            int originatingRound,
                                            LocalOptimizationLogger logger);
}
