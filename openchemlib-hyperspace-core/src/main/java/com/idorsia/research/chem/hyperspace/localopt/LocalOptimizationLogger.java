package com.idorsia.research.chem.hyperspace.localopt;

import java.util.List;

class LocalOptimizationLogger {

    private final LocalOptimizationLogLevel logLevel;
    private final String seedLabel;

    LocalOptimizationLogger(LocalOptimizationLogLevel logLevel, List<String> seedFragments) {
        this.logLevel = logLevel;
        this.seedLabel = String.join(";", seedFragments);
    }

    void logImprovement(int round, double score, String idcode) {
        if (logLevel.isImprovementsEnabled()) {
            System.out.println("[LocalOpt] seed=" + seedLabel + " round=" + round + " newBest=" +
                    String.format("%.4f", score) + " idcode=" + idcode);
        }
    }

    void logCandidate(List<String> fragments, double score, String idcode, int positionIndex) {
        if (logLevel.isVerboseEnabled()) {
            System.out.println("[LocalOpt] seed=" + seedLabel + " pos=" + positionIndex + " candidate=" +
                    String.join(";", fragments) + " score=" + String.format("%.4f", score) + " idcode=" + idcode);
        }
    }
}
