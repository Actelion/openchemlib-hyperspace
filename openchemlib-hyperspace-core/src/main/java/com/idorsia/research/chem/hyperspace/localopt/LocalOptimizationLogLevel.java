package com.idorsia.research.chem.hyperspace.localopt;

public enum LocalOptimizationLogLevel {
    NONE,
    SUMMARY,
    IMPROVEMENTS,
    VERBOSE;

    public boolean isSummaryEnabled() {
        return this == SUMMARY || this == IMPROVEMENTS || this == VERBOSE;
    }

    public boolean isImprovementsEnabled() {
        return this == IMPROVEMENTS || this == VERBOSE;
    }

    public boolean isVerboseEnabled() {
        return this == VERBOSE;
    }
}
