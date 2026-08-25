package com.idorsia.research.chem.hyperspace3d.index.build;

import java.nio.file.Path;

public record MoleculeFingerprintIndexBuildResult(
        Path manifestPath, long sourceRowCount, long recordCount,
        int shardCount, boolean resumed,
        MoleculeFingerprintIndexBuildMetrics metrics) {
}
