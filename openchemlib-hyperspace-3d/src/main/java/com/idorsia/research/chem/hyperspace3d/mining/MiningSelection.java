package com.idorsia.research.chem.hyperspace3d.mining;

import java.util.LinkedHashMap;
import java.util.Map;

/** Stable selection-mask contract for query-centric PheSA mining. */
public final class MiningSelection {
    public static final int BUCKETS = 6;
    public static final int TOTAL_OFFSET = 0;
    public static final int SHAPE_OFFSET = 6;
    public static final int PHARMACOPHORE_OFFSET = 12;
    public static final long RANDOM = 1L << 18;
    public static final long FLEXIBLE = 1L << 19;
    public static final long HIGH_SP3 = 1L << 20;
    public static final long STEREOCHEMICAL = 1L << 21;
    public static final long SIZE_EXTREME = 1L << 22;
    public static final long MODEL_DISAGREEMENT = 1L << 23;

    private MiningSelection() {}

    public static long rankBucket(int channelOffset, int bucket) {
        if ((channelOffset != TOTAL_OFFSET && channelOffset != SHAPE_OFFSET
                && channelOffset != PHARMACOPHORE_OFFSET)
                || bucket < 0 || bucket >= BUCKETS) {
            throw new IllegalArgumentException("invalid rank selection bit");
        }
        return 1L << (channelOffset + bucket);
    }

    public static Map<String, Long> manifestNames() {
        Map<String, Long> values = new LinkedHashMap<>();
        String[] channels = {"TOTAL", "SHAPE", "PHARMACOPHORE"};
        int[] offsets = {TOTAL_OFFSET, SHAPE_OFFSET, PHARMACOPHORE_OFFSET};
        for (int channel = 0; channel < channels.length; channel++) {
            for (int bucket = 0; bucket < BUCKETS; bucket++) {
                values.put(channels[channel] + "_RANK_BUCKET_" + bucket,
                        rankBucket(offsets[channel], bucket));
            }
        }
        values.put("RANDOM", RANDOM);
        values.put("FLEXIBLE", FLEXIBLE);
        values.put("HIGH_SP3", HIGH_SP3);
        values.put("STEREOCHEMICAL", STEREOCHEMICAL);
        values.put("SIZE_EXTREME", SIZE_EXTREME);
        values.put("MODEL_DISAGREEMENT", MODEL_DISAGREEMENT);
        return values;
    }
}
