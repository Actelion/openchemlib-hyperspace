package com.idorsia.research.chem.hyperspace3d.index.build;

import java.nio.file.Path;

public record ProductFingerprintIndexBuildResult(Path manifestPath, long recordCount,
        int shardCount, boolean resumed, ProductFingerprintIndexBuildMetrics metrics) {}
