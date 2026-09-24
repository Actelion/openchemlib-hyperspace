package com.idorsia.research.chem.hyperspace.rawspace;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CancellationException;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

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

        int workers = Math.min(applied.threads, Math.max(1, rawSpace.getReactions().size()));
        System.out.println("Building " + rawSpace.getReactions().size() + " reactions with " + workers + " workers");
        java.util.concurrent.atomic.AtomicInteger workerIds = new java.util.concurrent.atomic.AtomicInteger();
        ExecutorService executor = Executors.newFixedThreadPool(workers,
                task -> new Thread(task, "rawspace-build-" + workerIds.incrementAndGet()));
        CompletionService<String> completed = new ExecutorCompletionService<>(executor);
        List<Future<String>> futures = new ArrayList<>();
        try {
            for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : rawSpace.getReactions().entrySet()) {
                futures.add(completed.submit(() -> {
                    try {
                        buildReaction(space, entry, applied);
                        return entry.getKey();
                    } catch (Exception ex) {
                        throw new IllegalStateException("Failed to build reaction " + entry.getKey(), ex);
                    }
                }));
            }
            for (int i = 0; i < futures.size(); i++) {
                String reaction = completed.take().get();
                System.out.println("Built reaction " + reaction + " (" + (i + 1) + "/" + futures.size() + ")");
            }
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();
            throw ex;
        } catch (ExecutionException ex) {
            throw new IllegalStateException("Reaction build failed", ex.getCause());
        } finally {
            for (Future<String> future : futures) {
                if (!future.isDone()) future.cancel(true);
            }
            executor.shutdownNow();
            // Join even on interruption so no worker can continue mutating an abandoned index.
            boolean interrupted = Thread.interrupted();
            while (!executor.isTerminated()) {
                try {
                    executor.awaitTermination(1, TimeUnit.SECONDS);
                } catch (InterruptedException ex) {
                    interrupted = true;
                }
            }
            if (interrupted) Thread.currentThread().interrupt();
        }

        space.initAfterJavaDeserialization();
        space.reinitHelperMaps();
        space.reinitBitTree();
        return space;
    }

    private static void buildReaction(SynthonSpace space,
                                      Map.Entry<String, RawSynthonSpace.ReactionData> entry,
                                      BuildOptions applied) throws Exception {
        Map<Integer, List<Object>> molecules = new HashMap<>();
        Map<String, String> idcodeToId = new HashMap<>();
        entry.getValue().getRawFragmentSets().forEach((fragIdx, fragments) -> {
            List<Object> set = new ArrayList<>(fragments.size());
            for (RawSynthon frag : fragments) {
                if (Thread.currentThread().isInterrupted()) throw new CancellationException();
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
        private final int threads;

        private BuildOptions(Builder builder) {
            this.descriptorShortName = builder.descriptorShortName;
            this.descriptorBits = builder.descriptorBits;
            this.validateReactions = builder.validateReactions;
            this.threads = builder.threads;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static final class Builder {
            private String descriptorShortName;
            private int descriptorBits = -1;
            private boolean validateReactions = true;
            private int threads = 1;

            public Builder threads(int threads) {
                if (threads < 1) throw new IllegalArgumentException("threads must be positive");
                this.threads = threads;
                return this;
            }

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
