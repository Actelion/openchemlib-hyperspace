package com.idorsia.research.chem.hyperspace3d.batch;

import java.util.List;

public interface BatchAssemblyScorer {
    List<ScoredAssemblyCandidate> scoreBatch(List<AssemblyCandidate> candidates);
}
