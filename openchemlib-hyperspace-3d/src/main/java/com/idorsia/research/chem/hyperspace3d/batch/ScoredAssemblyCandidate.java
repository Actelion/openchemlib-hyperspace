package com.idorsia.research.chem.hyperspace3d.batch;

import java.util.Arrays;

public record ScoredAssemblyCandidate(AssemblyCandidate candidate, Status status,
                                      String rejectionReason, double objective,
                                      float[] scoreComponents) {
    public enum Status { SCORED, REJECTED }
    public ScoredAssemblyCandidate {
        scoreComponents = scoreComponents == null ? null
                : Arrays.copyOf(scoreComponents, scoreComponents.length);
    }
    public static ScoredAssemblyCandidate rejected(AssemblyCandidate candidate, String reason) {
        return new ScoredAssemblyCandidate(candidate, Status.REJECTED, reason,
                Double.NaN, null);
    }
}
