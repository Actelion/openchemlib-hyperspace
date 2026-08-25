package com.idorsia.research.chem.hyperspace3d.screening;

import java.util.List;

public record MoleculeSkelSpheresScreeningResult(List<MoleculeSkelSpheresHit> hits,
        long recordsScanned, long scanNanos, long metadataResolutionNanos) {}
