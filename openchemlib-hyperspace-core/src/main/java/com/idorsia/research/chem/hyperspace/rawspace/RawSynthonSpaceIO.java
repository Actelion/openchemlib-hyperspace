package com.idorsia.research.chem.hyperspace.rawspace;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingRequest;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * JSON serializer/deserializer for {@link RawSynthonSpace}.
 */
public final class RawSynthonSpaceIO {

    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT);

    private RawSynthonSpaceIO() {
    }

    public static void write(RawSynthonSpace space, Path path) throws IOException {
        Path parent = path.toAbsolutePath().getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        JsonRawSpace json = toJson(space);
        try (OutputStream out = newOutputStream(path)) {
            MAPPER.writeValue(out, json);
        }
    }

    public static RawSynthonSpace read(Path path) throws IOException {
        JsonRawSpace json;
        try (InputStream in = newInputStream(path)) {
            json = MAPPER.readValue(in, JsonRawSpace.class);
        }
        return toDomain(json);
    }

    private static JsonRawSpace toJson(RawSynthonSpace space) {
        JsonRawSpace json = new JsonRawSpace();
        json.name = space.getName();
        json.version = space.getVersion();
        json.metadata = space.getMetadata();
        json.downsamplingAlgorithm = space.getDownsamplingAlgorithm().orElse(null);
        json.downsamplingRequest = space.getDownsamplingRequest()
                .map(RawSynthonSpaceIO::toJsonRequest)
                .orElse(null);
        json.reactions = new ArrayList<>();
        space.getReactions().forEach((rxnId, data) -> json.reactions.add(toJsonReaction(rxnId, data)));
        return json;
    }

    private static JsonReaction toJsonReaction(String reactionId, RawSynthonSpace.ReactionData data) {
        JsonReaction reaction = new JsonReaction();
        reaction.reactionId = reactionId;
        reaction.fragmentSets = encodeFragmentMap(data.getRawFragmentSets());
        reaction.downsampledFragmentSets = encodeFragmentMap(data.getRawDownsampledSets());
        reaction.exampleScaffolds = data.getExampleScaffolds();
        reaction.partialAssemblies = stringifyKeys(data.getPartialAssemblies());
        reaction.representativeCompounds = data.getRepresentativeCompounds();
        reaction.descriptors = data.getDescriptors();
        reaction.fragmentAttributes = data.getFragmentAttributes();
        return reaction;
    }

    private static Map<String, List<JsonFragment>> encodeFragmentMap(Map<Integer, List<RawSynthon>> map) {
        Map<String, List<JsonFragment>> encoded = new LinkedHashMap<>();
        map.forEach((idx, list) -> {
            List<JsonFragment> fragments = new ArrayList<>();
            for (RawSynthon frag : list) {
                fragments.add(toJsonFragment(idx, frag));
            }
            encoded.put(Integer.toString(idx), fragments);
        });
        return encoded;
    }

    private static Map<String, List<String>> stringifyKeys(Map<Integer, List<String>> map) {
        Map<String, List<String>> result = new LinkedHashMap<>();
        map.forEach((idx, list) -> result.put(Integer.toString(idx), list));
        return result;
    }

    private static JsonFragment toJsonFragment(int fragIdx, RawSynthon fragment) {
        JsonFragment json = new JsonFragment();
        json.reactionId = fragment.getReactionId();
        json.fragIndex = fragIdx;
        json.fragmentId = fragment.getFragmentId();
        json.idcode = fragment.getIdcode();
        json.connectors = encodeBitSet(fragment.getConnectors());
        return json;
    }

    private static RawSynthonSpace toDomain(JsonRawSpace json) {
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(json.name)
                .version(json.version);
        if (json.metadata != null) {
            json.metadata.forEach(builder::putMetadata);
        }
        if (json.reactions != null) {
            for (JsonReaction reaction : json.reactions) {
                if (reaction.fragmentSets != null) {
                    reaction.fragmentSets.forEach((key, list) -> builder.addRawFragments(reaction.reactionId,
                            Integer.parseInt(key),
                            decodeFragments(list)));
                }
                if (reaction.downsampledFragmentSets != null) {
                    reaction.downsampledFragmentSets.forEach((key, list) -> builder.addRawDownsampledFragments(reaction.reactionId,
                            Integer.parseInt(key),
                            decodeFragments(list)));
                }
                if (reaction.exampleScaffolds != null) {
                    builder.addExampleScaffolds(reaction.reactionId, reaction.exampleScaffolds);
                }
                if (reaction.partialAssemblies != null) {
                    reaction.partialAssemblies.forEach((key, assemblies) -> builder.addPartialAssemblies(reaction.reactionId,
                            Integer.parseInt(key),
                            assemblies));
                }
                if (reaction.representativeCompounds != null) {
                    builder.addRepresentativeCompounds(reaction.reactionId, reaction.representativeCompounds);
                }
                if (reaction.descriptors != null) {
                    reaction.descriptors.forEach((k, v) -> builder.addReactionDescriptor(reaction.reactionId, k, v));
                }
                if (reaction.fragmentAttributes != null) {
                    reaction.fragmentAttributes.forEach((fragmentId, attributes) -> {
                        if (attributes != null) {
                            attributes.forEach((key, value) ->
                                    builder.addFragmentAttribute(reaction.reactionId, fragmentId, key, value));
                        }
                    });
                }
            }
        }
        if (json.downsamplingAlgorithm != null || json.downsamplingRequest != null) {
            SynthonDownsamplingRequest request = json.downsamplingRequest == null
                    ? SynthonDownsamplingRequest.builder().build()
                    : fromJsonRequest(json.downsamplingRequest);
            builder.withDownsamplingMetadata(json.downsamplingAlgorithm, request);
        }
        return builder.build();
    }

    private static OutputStream newOutputStream(Path path) throws IOException {
        OutputStream out = new BufferedOutputStream(Files.newOutputStream(path));
        if (path.getFileName().toString().toLowerCase().endsWith(".gz")) {
            return new GZIPOutputStream(out);
        }
        return out;
    }

    private static InputStream newInputStream(Path path) throws IOException {
        InputStream in = new BufferedInputStream(Files.newInputStream(path));
        if (path.getFileName().toString().toLowerCase().endsWith(".gz")) {
            return new GZIPInputStream(in);
        }
        return in;
    }

    private static List<RawSynthon> decodeFragments(List<JsonFragment> fragments) {
        List<RawSynthon> list = new ArrayList<>();
        for (JsonFragment fragment : fragments) {
            list.add(new RawSynthon(fragment.reactionId,
                    fragment.fragIndex,
                    fragment.fragmentId,
                    fragment.idcode,
                    decodeBitSet(fragment.connectors)));
        }
        return list;
    }

    private static String encodeBitSet(java.util.BitSet bitSet) {
        if (bitSet == null || bitSet.isEmpty()) {
            return null;
        }
        return Base64.getEncoder().encodeToString(bitSet.toByteArray());
    }

    private static java.util.BitSet decodeBitSet(String encoded) {
        if (encoded == null || encoded.isEmpty()) {
            return new java.util.BitSet();
        }
        return java.util.BitSet.valueOf(Base64.getDecoder().decode(encoded));
    }

    private static final class JsonRawSpace {
        public String name;
        public String version;
        public Map<String, String> metadata;
        public String downsamplingAlgorithm;
        public JsonDownsamplingRequest downsamplingRequest;
        public List<JsonReaction> reactions;
    }

    private static final class JsonReaction {
        public String reactionId;
        public Map<String, List<JsonFragment>> fragmentSets;
        public Map<String, List<JsonFragment>> downsampledFragmentSets;
        public List<String> exampleScaffolds;
        public Map<String, List<String>> partialAssemblies;
        public List<String> representativeCompounds;
        public Map<String, String> descriptors;
        public Map<String, Map<String, String>> fragmentAttributes;
    }

    private static final class JsonFragment {
        public String reactionId;
        public int fragIndex;
        public String fragmentId;
        public String idcode;
        public String connectors;
    }

    private static final class JsonDownsamplingRequest {
        public int maxCenters;
        public double sizeCapScale;
        public double sizeCapOffset;
        public double minSimilarity;
        public long randomSeed;
        public boolean enforceConnectorEquivalence;
        public Map<String, Object> attributes;
    }

    private static JsonDownsamplingRequest toJsonRequest(SynthonDownsamplingRequest request) {
        JsonDownsamplingRequest json = new JsonDownsamplingRequest();
        json.maxCenters = request.getMaxCenters();
        json.sizeCapScale = request.getSizeCapScale();
        json.sizeCapOffset = request.getSizeCapOffset();
        json.minSimilarity = request.getMinSimilarity();
        json.randomSeed = request.getRandomSeed();
        json.enforceConnectorEquivalence = request.isEnforceConnectorEquivalence();
        json.attributes = request.getAttributes();
        return json;
    }

    private static SynthonDownsamplingRequest fromJsonRequest(JsonDownsamplingRequest json) {
        SynthonDownsamplingRequest.Builder builder = SynthonDownsamplingRequest.builder()
                .withMaxCenters(json.maxCenters)
                .withSizeCapScale(json.sizeCapScale)
                .withSizeCapOffset(json.sizeCapOffset)
                .withMinSimilarity(json.minSimilarity)
                .withRandomSeed(json.randomSeed)
                .enforceConnectorEquivalence(json.enforceConnectorEquivalence);
        if (json.attributes != null) {
            json.attributes.forEach(builder::putAttribute);
        }
        return builder.build();
    }
}
