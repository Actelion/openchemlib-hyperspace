package com.idorsia.research.chem.hyperspace3d.screening;

import java.util.List;
import java.util.Map;

public record CompactSkelSpheresScreeningResult(List<CompactSkelSpheresHit> global,
        Map<String, List<CompactSkelSpheresHit>> byReaction, long recordsScanned,
        long scanNanos, long tupleResolutionNanos) {}
