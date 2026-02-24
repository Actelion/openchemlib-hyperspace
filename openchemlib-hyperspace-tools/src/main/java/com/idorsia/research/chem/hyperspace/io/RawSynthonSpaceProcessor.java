package com.idorsia.research.chem.hyperspace.io;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace.Builder;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace.ReactionData;
import com.idorsia.research.chem.hyperspace.downsampling.SynthonDownsamplingRequest;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Applies one or more processing tasks to a {@link RawSynthonSpace} and emits a new instance.
 */
public final class RawSynthonSpaceProcessor {

    private RawSynthonSpaceProcessor() {
    }

    public static RawSynthonSpace process(RawSynthonSpace source, List<Task> tasks) throws Exception {
        Objects.requireNonNull(source, "source");
        Objects.requireNonNull(tasks, "tasks");
        if (tasks.isEmpty()) {
            return source;
        }
        Builder builder = RawSynthonSpace.builder(source.getName()).version(source.getVersion());
        source.getMetadata().forEach(builder::putMetadata);
        source.getDownsamplingAlgorithm().ifPresent(algorithm -> {
            SynthonDownsamplingRequest request = source.getDownsamplingRequest().orElse(null);
            if (request != null) {
                builder.withDownsamplingMetadata(algorithm, request);
            }
        });
        copyReactions(source, builder);

        Context context = new Context(source, builder);
        for (Task task : tasks) {
            task.apply(context);
        }
        context.flushMetadata();
        return builder.build();
    }

    private static void copyReactions(RawSynthonSpace source, Builder builder) {
        source.getReactions().forEach((reactionId, data) -> {
            data.getRawFragmentSets().forEach((idx, fragments) -> builder.addRawFragments(reactionId, idx, fragments));
            data.getRawDownsampledSets().forEach((idx, fragments) -> builder.addRawDownsampledFragments(reactionId, idx, fragments));
            builder.addExampleScaffolds(reactionId, data.getExampleScaffolds());
            data.getPartialAssemblies().forEach((missingIdx, assemblies) ->
                    builder.addPartialAssemblies(reactionId, missingIdx, assemblies));
            builder.addRepresentativeCompounds(reactionId, data.getRepresentativeCompounds());
            data.getDescriptors().forEach((key, value) -> builder.addReactionDescriptor(reactionId, key, value));
            data.getFragmentAttributes().forEach((fragmentId, attributes) -> {
                if (attributes == null) {
                    return;
                }
                attributes.forEach((key, value) -> builder.addFragmentAttribute(reactionId, fragmentId, key, value));
            });
        });
    }

    public interface Task {
        void apply(Context context) throws Exception;
    }

    public static final class Context {
        private final RawSynthonSpace source;
        private final RawSynthonSpace.Builder builder;
        private final LinkedHashSet<String> descriptorTags;

        private Context(RawSynthonSpace source, RawSynthonSpace.Builder builder) {
            this.source = source;
            this.builder = builder;
            this.descriptorTags = new LinkedHashSet<>(parseDescriptorTags(
                    source.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS)));
        }

        public RawSynthonSpace getSource() {
            return source;
        }

        public RawSynthonSpace.Builder getBuilder() {
            return builder;
        }

        public void addDescriptorTag(String tag) {
            if (tag != null && !tag.isBlank()) {
                descriptorTags.add(tag.trim());
            }
        }

        public void putMetadata(String key, String value) {
            builder.putMetadata(key, value);
        }

        public boolean hasFragmentAttribute(String reactionId, String fragmentId, String key) {
            ReactionData data = source.getReactions().get(reactionId);
            if (data == null) {
                return false;
            }
            Map<String, Map<String, String>> attributes = data.getFragmentAttributes();
            Map<String, String> fragAttributes = attributes.get(fragmentId);
            if (fragAttributes == null) {
                return false;
            }
            return fragAttributes.containsKey(key);
        }

        public void addFragmentAttribute(String reactionId, String fragmentId, String key, String value) {
            builder.addFragmentAttribute(reactionId, fragmentId, key, value);
        }

        private void flushMetadata() {
            if (!descriptorTags.isEmpty()) {
                builder.putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS, String.join(",", descriptorTags));
            }
        }

        private static List<String> parseDescriptorTags(String raw) {
            List<String> tags = new ArrayList<>();
            if (raw == null || raw.isBlank()) {
                return tags;
            }
            for (String token : raw.split(",")) {
                if (token != null && !token.isBlank()) {
                    tags.add(token.trim());
                }
            }
            return tags;
        }
    }
}
