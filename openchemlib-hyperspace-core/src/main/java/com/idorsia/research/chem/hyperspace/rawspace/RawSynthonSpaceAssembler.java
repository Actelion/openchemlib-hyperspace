package com.idorsia.research.chem.hyperspace.rawspace;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Utility methods to rehydrate {@link SynthonSpace} instances from {@link RawSynthonSpace} dumps.
 */
public final class RawSynthonSpaceAssembler {

    private static final int DEFAULT_BITS = 1024;

    private RawSynthonSpaceAssembler() {
    }

    public static SynthonSpace buildSynthonSpace(RawSynthonSpace rawSpace) throws Exception {
        return buildSynthonSpace(rawSpace, BuildOptions.builder().build());
    }

    public static SynthonSpace buildSynthonSpace(RawSynthonSpace rawSpace, BuildOptions options) throws Exception {
        BuildOptions applied = options == null ? BuildOptions.builder().build() : options;
        if (applied.descriptorShortName == null || applied.descriptorShortName.isBlank()) {
            throw new IllegalArgumentException("Descriptor short name is required when materializing a SynthonSpace");
        }
        String descriptorName = applied.descriptorShortName;
        int bits = applied.descriptorBits > 0
                ? applied.descriptorBits
                : parseBits(rawSpace.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_BITS));

        SynthonSpace space = new SynthonSpace();
        space.setFP(SynthonSpace.resolveDescriptorHandlerFromName(descriptorName), bits);

        for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : rawSpace.getReactions().entrySet()) {
            Map<Integer, List<Object>> molecules = new HashMap<>();
            Map<String, String> idcodeToId = new HashMap<>();
            entry.getValue().getRawFragmentSets().forEach((fragIdx, fragments) -> {
                List<Object> set = new ArrayList<>(fragments.size());
                for (RawSynthon frag : fragments) {
                    set.add(frag.getIdcode());
                    if (frag.getFragmentId() != null) {
                        idcodeToId.put(frag.getIdcode(), frag.getFragmentId());
                    }
                }
                molecules.put(fragIdx, set);
            });

            if (applied.validateReactions) {
                SynthonReactionValidator.validate(molecules);
            }

            space.addReaction(entry.getKey(), molecules, idcodeToId, null);
        }

        space.initAfterJavaDeserialization();
        space.reinitHelperMaps();
        space.reinitBitTree();
        return space;
    }

    private static int parseBits(String rawBits) {
        if (rawBits == null || rawBits.isBlank()) {
            return DEFAULT_BITS;
        }
        try {
            return Integer.parseInt(rawBits);
        } catch (NumberFormatException ex) {
            return DEFAULT_BITS;
        }
    }

    public static final class BuildOptions {
        private final String descriptorShortName;
        private final int descriptorBits;
        private final boolean validateReactions;

        private BuildOptions(Builder builder) {
            this.descriptorShortName = builder.descriptorShortName;
            this.descriptorBits = builder.descriptorBits;
            this.validateReactions = builder.validateReactions;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static final class Builder {
            private String descriptorShortName;
            private int descriptorBits = -1;
            private boolean validateReactions = true;

            public Builder descriptorShortName(String name) {
                this.descriptorShortName = name;
                return this;
            }

            public Builder descriptorBits(int bits) {
                this.descriptorBits = bits;
                return this;
            }

            public Builder validateReactions(boolean validate) {
                this.validateReactions = validate;
                return this;
            }

            public BuildOptions build() {
                if (descriptorShortName == null || descriptorShortName.isBlank()) {
                    throw new IllegalArgumentException("descriptorShortName must be provided");
                }
                return new BuildOptions(this);
            }
        }
    }
}
