package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Side-car structure that mirrors a SynthonSpace with downsampled sets.
 */
public class DownsampledSynthonSpace implements Serializable {

    private final String algorithmName;
    private final SynthonDownsamplingRequest request;
    private final Map<String, Map<Integer, List<SynthonSpace.FragId>>> downsampledSets;

    public DownsampledSynthonSpace(String algorithmName,
                                   SynthonDownsamplingRequest request,
                                   Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> representatives) {
        this.algorithmName = Objects.requireNonNull(algorithmName, "algorithmName");
        this.request = Objects.requireNonNull(request, "request");
        this.downsampledSets = buildNestedMap(Objects.requireNonNull(representatives, "representatives"));
    }

    public static DownsampledSynthonSpace fromRawFragmentSets(RawSynthonSpace rawSpace) {
        Objects.requireNonNull(rawSpace, "rawSpace");
        Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> representatives = new HashMap<>();
        rawSpace.getReactions().forEach((rxnId, data) ->
                data.getFragmentSets().forEach((idx, synthons) ->
                        representatives.put(new SynthonSpace.FragType(rxnId, idx), synthons)));
        String algorithm = rawSpace.getMetadata().getOrDefault(
                RawSynthonSpace.MetadataKeys.DOWNSAMPLING_ALGORITHM,
                "RawSynthonSpace");
        return new DownsampledSynthonSpace(algorithm, requestFromMetadata(rawSpace), representatives);
    }

    private static SynthonDownsamplingRequest requestFromMetadata(RawSynthonSpace rawSpace) {
        Map<String, String> metadata = rawSpace.getMetadata();
        SynthonDownsamplingRequest.Builder builder = SynthonDownsamplingRequest.builder()
                .withMaxCenters(parseInt(metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_MAX_CENTERS), 0))
                .withSizeCapScale(parseDouble(metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SIZE_CAP_SCALE), 0.0))
                .withSizeCapOffset(parseDouble(metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SIZE_CAP_OFFSET), 0.0))
                .withMinSimilarity(parseDouble(metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_MIN_SIMILARITY), 0.0))
                .withRandomSeed(parseLong(metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_SEED), 1L))
                .enforceConnectorEquivalence(parseBoolean(
                        metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_ENFORCE_CONNECTOR_EQUIVALENCE), true))
                .includeClusterMembers(parseBoolean(
                        metadata.get(RawSynthonSpace.MetadataKeys.DOWNSAMPLING_INCLUDE_CLUSTER_MEMBERS), false));
        return builder.build();
    }

    private static int parseInt(String value, int defaultValue) {
        if (value == null || value.isBlank()) {
            return defaultValue;
        }
        return Integer.parseInt(value);
    }

    private static long parseLong(String value, long defaultValue) {
        if (value == null || value.isBlank()) {
            return defaultValue;
        }
        return Long.parseLong(value);
    }

    private static double parseDouble(String value, double defaultValue) {
        if (value == null || value.isBlank()) {
            return defaultValue;
        }
        return Double.parseDouble(value);
    }

    private static boolean parseBoolean(String value, boolean defaultValue) {
        if (value == null || value.isBlank()) {
            return defaultValue;
        }
        return Boolean.parseBoolean(value);
    }

    private static Map<String, Map<Integer, List<SynthonSpace.FragId>>> buildNestedMap(Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> representatives) {
        Map<String, Map<Integer, List<SynthonSpace.FragId>>> nested = new HashMap<>();
        for (Map.Entry<SynthonSpace.FragType, List<SynthonSpace.FragId>> entry : representatives.entrySet()) {
            SynthonSpace.FragType fragType = entry.getKey();
            nested.computeIfAbsent(fragType.rxn_id, key -> new HashMap<>())
                    .put(fragType.frag, Collections.unmodifiableList(new ArrayList<>(entry.getValue())));
        }
        Map<String, Map<Integer, List<SynthonSpace.FragId>>> immutableOuter = new HashMap<>();
        for (Map.Entry<String, Map<Integer, List<SynthonSpace.FragId>>> entry : nested.entrySet()) {
            immutableOuter.put(entry.getKey(), Collections.unmodifiableMap(entry.getValue()));
        }
        return Collections.unmodifiableMap(immutableOuter);
    }

    public String getAlgorithmName() {
        return algorithmName;
    }

    public SynthonDownsamplingRequest getRequest() {
        return request;
    }

    public Map<String, Map<Integer, List<SynthonSpace.FragId>>> getDownsampledSets() {
        return downsampledSets;
    }

    public List<SynthonSpace.FragId> getSynthonSet(String rxnId, int fragIdx) {
        Map<Integer, List<SynthonSpace.FragId>> rxn = downsampledSets.get(rxnId);
        if (rxn == null) {
            return Collections.emptyList();
        }
        return rxn.getOrDefault(fragIdx, Collections.emptyList());
    }
}
