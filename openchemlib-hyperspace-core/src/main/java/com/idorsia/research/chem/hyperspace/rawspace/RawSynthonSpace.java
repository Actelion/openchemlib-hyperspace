package com.idorsia.research.chem.hyperspace.rawspace;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingRequest;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Lightweight representation of synthon reactions and fragments enriched with optional metadata.
 */
public class RawSynthonSpace implements Serializable {

    public static final class MetadataKeys {
        public static final String DESCRIPTOR_SHORT_NAME = "descriptor.shortName";
        public static final String DESCRIPTOR_BITS = "descriptor.bits";
        public static final String SOURCE_FORMAT = "source.format";
        public static final String DESCRIPTOR_TAGS = "descriptor.tags";
        public static final String SPACE_ROLE = "space.role";
        public static final String SPACE_ROLE_FULL = "full";
        public static final String SPACE_ROLE_DOWNSAMPLED = "downsampled";
        public static final String DOWNSAMPLING_ALGORITHM = "downsampling.algorithm";
        public static final String DOWNSAMPLING_MAX_CENTERS = "downsampling.maxCenters";
        public static final String DOWNSAMPLING_MIN_SIMILARITY = "downsampling.minSimilarity";
        public static final String DOWNSAMPLING_SEED = "downsampling.seed";
        public static final String DOWNSAMPLING_ENFORCE_CONNECTOR_EQUIVALENCE = "downsampling.enforceConnectorEquivalence";
        public static final String DOWNSAMPLING_SIZE_CAP_SCALE = "downsampling.sizeCapScale";
        public static final String DOWNSAMPLING_SIZE_CAP_OFFSET = "downsampling.sizeCapOffset";
        public static final String DOWNSAMPLING_INCLUDE_CLUSTER_MEMBERS = "downsampling.includeClusterMembers";

        private MetadataKeys() {
        }
    }

    private final String name;
    private final String version;
    private final Map<String, String> metadata;
    private final Map<String, ReactionData> reactions;

    private RawSynthonSpace(Builder builder) {
        this.name = builder.name;
        this.version = builder.version;
        this.metadata = Collections.unmodifiableMap(new LinkedHashMap<>(builder.metadata));
        Map<String, ReactionData> map = new LinkedHashMap<>();
        builder.reactionBuilders.forEach((id, reactionBuilder) -> map.put(id, reactionBuilder.build()));
        this.reactions = Collections.unmodifiableMap(map);
    }

    public String getName() {
        return name;
    }

    public String getVersion() {
        return version;
    }

    public Map<String, String> getMetadata() {
        return metadata;
    }

    public Map<String, ReactionData> getReactions() {
        return reactions;
    }

    public Map<Integer, List<SynthonSpace.FragId>> getFragmentSets(String reactionId) {
        ReactionData data = reactions.get(reactionId);
        if (data == null) {
            return Collections.emptyMap();
        }
        return data.getFragmentSets();
    }

    public SynthonSpace.FragId findFragment(String reactionId, String fragmentId) {
        ReactionData data = reactions.get(reactionId);
        if (data == null) {
            return null;
        }
        return data.findFragment(fragmentId);
    }

    public static Builder builder(String name) {
        return new Builder(name);
    }

    public static final class Builder {
        private final String name;
        private String version = "1.0";
        private final Map<String, String> metadata = new LinkedHashMap<>();
        private final Map<String, ReactionBuilder> reactionBuilders = new LinkedHashMap<>();

        private Builder(String name) {
            this.name = Objects.requireNonNull(name, "name");
        }

        public Builder version(String version) {
            this.version = Objects.requireNonNull(version, "version");
            return this;
        }

        public Builder putMetadata(String key, String value) {
            if (key != null && value != null) {
                metadata.put(key, value);
            }
            return this;
        }

        public Builder addRawFragments(String reactionId, int fragIdx, List<RawSynthon> fragments) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addRawFragments(fragIdx, fragments);
            return this;
        }

        public Builder addFragments(String reactionId, int fragIdx, List<SynthonSpace.FragId> fragments) {
            return addRawFragments(reactionId, fragIdx, toRawSynthons(fragments));
        }

        public Builder addExampleScaffolds(String reactionId, List<String> idcodes) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addExampleScaffolds(idcodes);
            return this;
        }

        public Builder addPartialAssemblies(String reactionId, int missingIndex, List<String> assemblies) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addPartialAssemblies(missingIndex, assemblies);
            return this;
        }

        public Builder addRepresentativeCompounds(String reactionId, List<String> compounds) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addRepresentativeCompounds(compounds);
            return this;
        }

        public Builder addReactionDescriptor(String reactionId, String key, String value) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addDescriptor(key, value);
            return this;
        }

        public Builder addReactionMetadata(String reactionId, String key, String value) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addMetadata(key, value);
            return this;
        }

        public Builder addFragmentAttribute(String reactionId, String fragmentId, String key, String value) {
            reactionBuilders.computeIfAbsent(reactionId, ReactionBuilder::new)
                    .addFragmentAttribute(fragmentId, key, value);
            return this;
        }

        public Builder withFullSynthonSpace(SynthonSpace space) {
            for (String rxnId : space.getRxnIds()) {
                Map<Integer, SynthonSpace.FragType> fragTypes = space.getFragTypes(rxnId);
                if (fragTypes == null) {
                    continue;
                }
                for (Integer idx : fragTypes.keySet()) {
                    List<SynthonSpace.FragId> synthons = space.getSynthonSet(rxnId, idx);
                    if (synthons == null || synthons.isEmpty()) {
                        continue;
                    }
                    addRawFragments(rxnId, idx, toRawSynthons(synthons));
                }
            }
            return this;
        }

        public Builder withDownsamplingMetadata(String algorithm,
                                                SynthonDownsamplingRequest request) {
            putMetadata(MetadataKeys.SPACE_ROLE, MetadataKeys.SPACE_ROLE_DOWNSAMPLED);
            putMetadata(MetadataKeys.DOWNSAMPLING_ALGORITHM, algorithm);
            if (request != null) {
                putMetadata(MetadataKeys.DOWNSAMPLING_MAX_CENTERS, Integer.toString(request.getMaxCenters()));
                putMetadata(MetadataKeys.DOWNSAMPLING_MIN_SIMILARITY, Double.toString(request.getMinSimilarity()));
                putMetadata(MetadataKeys.DOWNSAMPLING_SEED, Long.toString(request.getRandomSeed()));
                putMetadata(MetadataKeys.DOWNSAMPLING_ENFORCE_CONNECTOR_EQUIVALENCE,
                        Boolean.toString(request.isEnforceConnectorEquivalence()));
                putMetadata(MetadataKeys.DOWNSAMPLING_SIZE_CAP_SCALE, Double.toString(request.getSizeCapScale()));
                putMetadata(MetadataKeys.DOWNSAMPLING_SIZE_CAP_OFFSET, Double.toString(request.getSizeCapOffset()));
                putMetadata(MetadataKeys.DOWNSAMPLING_INCLUDE_CLUSTER_MEMBERS,
                        Boolean.toString(request.isIncludeClusterMembers()));
            }
            return this;
        }

        public RawSynthonSpace build() {
            return new RawSynthonSpace(this);
        }

        private List<RawSynthon> toRawSynthons(List<SynthonSpace.FragId> fragments) {
            List<RawSynthon> raws = new ArrayList<>(fragments.size());
            for (SynthonSpace.FragId frag : fragments) {
                raws.add(RawSynthon.fromFragId(frag));
            }
            return raws;
        }
    }


    public static final class ReactionData implements Serializable {
        private final Map<Integer, List<RawSynthon>> rawFragmentSets;
        private final List<String> exampleScaffolds;
        private final Map<Integer, List<String>> partialAssemblies;
        private final List<String> representativeCompounds;
        private final Map<String, String> descriptors;
        private final Map<String, String> reactionMetadata;
        private final Map<String, RawSynthon> fragmentLookup;
        private final Map<String, Map<String, String>> fragmentAttributes;
        private transient Map<Integer, List<SynthonSpace.FragId>> fragmentSetsView;

        private ReactionData(Map<Integer, List<RawSynthon>> fragmentSets,
                             List<String> exampleScaffolds,
                             Map<Integer, List<String>> partialAssemblies,
                             List<String> representativeCompounds,
                             Map<String, String> descriptors,
                             Map<String, String> reactionMetadata,
                             Map<String, Map<String, String>> fragmentAttributes) {
            this.rawFragmentSets = fragmentSets;
            this.exampleScaffolds = exampleScaffolds;
            this.partialAssemblies = partialAssemblies;
            this.representativeCompounds = representativeCompounds;
            this.descriptors = descriptors;
            this.reactionMetadata = reactionMetadata;
            this.fragmentAttributes = fragmentAttributes;
            Map<String, RawSynthon> lookup = new HashMap<>();
            fragmentSets.values().forEach(list -> list.forEach(f -> lookup.put(f.getFragmentId(), f)));
            this.fragmentLookup = Collections.unmodifiableMap(lookup);
        }

        public Map<Integer, List<SynthonSpace.FragId>> getFragmentSets() {
            if (fragmentSetsView == null) {
                fragmentSetsView = buildFragIdView(rawFragmentSets);
            }
            return fragmentSetsView;
        }

        public Map<Integer, List<RawSynthon>> getRawFragmentSets() {
            return rawFragmentSets;
        }

        public List<String> getExampleScaffolds() {
            return exampleScaffolds;
        }

        public Map<Integer, List<String>> getPartialAssemblies() {
            return partialAssemblies;
        }

        public List<String> getRepresentativeCompounds() {
            return representativeCompounds;
        }

        public Map<String, String> getDescriptors() {
            return descriptors;
        }

        public Map<String, String> getReactionMetadata() {
            return reactionMetadata;
        }

        public Map<String, Map<String, String>> getFragmentAttributes() {
            return fragmentAttributes;
        }

        public SynthonSpace.FragId findFragment(String fragmentId) {
            RawSynthon raw = fragmentLookup.get(fragmentId);
            return raw == null ? null : raw.toFragId();
        }

        public RawSynthon findRawFragment(String fragmentId) {
            return fragmentLookup.get(fragmentId);
        }

        private Map<Integer, List<SynthonSpace.FragId>> buildFragIdView(Map<Integer, List<RawSynthon>> source) {
            Map<Integer, List<SynthonSpace.FragId>> result = new HashMap<>();
            source.forEach((idx, list) -> {
                List<SynthonSpace.FragId> converted = new ArrayList<>(list.size());
                list.forEach(raw -> converted.add(raw.toFragId()));
                result.put(idx, Collections.unmodifiableList(converted));
            });
            return Collections.unmodifiableMap(result);
        }
    }

    private static final class ReactionBuilder {
        private final String reactionId;
        private final Map<Integer, List<RawSynthon>> fragmentSets = new HashMap<>();
        private final List<String> exampleScaffolds = new ArrayList<>();
        private final Map<Integer, List<String>> partialAssemblies = new HashMap<>();
        private final List<String> representativeCompounds = new ArrayList<>();
        private final Map<String, String> descriptors = new HashMap<>();
        private final Map<String, String> reactionMetadata = new HashMap<>();
        private final Map<String, Map<String, String>> fragmentAttributes = new HashMap<>();

        private ReactionBuilder(String reactionId) {
            this.reactionId = reactionId;
        }

        private void addRawFragments(int fragIdx, List<RawSynthon> synthons) {
            fragmentSets.put(fragIdx, immutableCopyRaw(synthons));
        }

        private void addExampleScaffolds(List<String> idcodes) {
            if (idcodes != null) {
                exampleScaffolds.addAll(idcodes);
            }
        }

        private void addPartialAssemblies(int missingIdx, List<String> assemblies) {
            if (assemblies == null || assemblies.isEmpty()) {
                return;
            }
            partialAssemblies.computeIfAbsent(missingIdx, key -> new ArrayList<>()).addAll(assemblies);
        }

        private void addRepresentativeCompounds(List<String> compounds) {
            if (compounds != null) {
                representativeCompounds.addAll(compounds);
            }
        }

        private void addDescriptor(String key, String value) {
            if (key != null && value != null) {
                descriptors.put(key, value);
            }
        }

        private void addMetadata(String key, String value) {
            if (key != null && value != null) {
                reactionMetadata.put(key, value);
            }
        }

        private void addFragmentAttribute(String fragmentId, String key, String value) {
            if (fragmentId == null || key == null || value == null) {
                return;
            }
            fragmentAttributes.computeIfAbsent(fragmentId, id -> new LinkedHashMap<>()).put(key, value);
        }

        private ReactionData build() {
            Map<Integer, List<RawSynthon>> fragments = wrapRawMap(fragmentSets);
            Map<Integer, List<String>> partials = new HashMap<>();
            partialAssemblies.forEach((idx, list) -> partials.put(idx, Collections.unmodifiableList(new ArrayList<>(list))));
            Map<String, Map<String, String>> attrs = new HashMap<>();
            fragmentAttributes.forEach((fragId, attr) -> attrs.put(fragId,
                    Collections.unmodifiableMap(new LinkedHashMap<>(attr))));
            return new ReactionData(fragments,
                    Collections.unmodifiableList(new ArrayList<>(exampleScaffolds)),
                    Collections.unmodifiableMap(partials),
                    Collections.unmodifiableList(new ArrayList<>(representativeCompounds)),
                    Collections.unmodifiableMap(new HashMap<>(descriptors)),
                    Collections.unmodifiableMap(new HashMap<>(reactionMetadata)),
                    Collections.unmodifiableMap(attrs));
        }

        private Map<Integer, List<RawSynthon>> wrapRawMap(Map<Integer, List<RawSynthon>> source) {
            Map<Integer, List<RawSynthon>> map = new HashMap<>();
            source.forEach((idx, list) -> map.put(idx, immutableCopyRaw(list)));
            return Collections.unmodifiableMap(map);
        }

        private List<RawSynthon> immutableCopyRaw(List<RawSynthon> synthons) {
            return Collections.unmodifiableList(new ArrayList<>(synthons));
        }
    }

}
