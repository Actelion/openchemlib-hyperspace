package com.idorsia.research.chem.hyperspace.localopt;

public enum LocalOptimizationLogLevel {
    NONE,
    IMPROVEMENTS,
    VERBOSE;

    public boolean isImprovementsEnabled() {
        return this == IMPROVEMENTS || this == VERBOSE;
    }

    public boolean isVerboseEnabled() {
        return this == VERBOSE;
    }
}
