package com.idorsia.research.chem.hyperspace.io;

import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongPFP512;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceProcessor.Context;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceProcessor.Task;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Base64;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Task that attaches serialized descriptors to each fragment as fragment attributes.
 */
public final class DescriptorAugmentationTask implements Task {

    private static final String ATTRIBUTE_PREFIX = "descriptor.";
    private static final int DEFAULT_BATCH_SIZE = 1_000;

    private final List<String> descriptorNames;
    private final int threads;
    private final boolean overwriteExisting;

    public DescriptorAugmentationTask(List<String> descriptorNames, int threads, boolean overwriteExisting) {
        this.descriptorNames = Objects.requireNonNull(descriptorNames, "descriptorNames");
        this.threads = Math.max(1, threads);
        this.overwriteExisting = overwriteExisting;
    }

    @Override
    public void apply(Context context) throws Exception {
        List<DescriptorJob> jobs = buildJobs(descriptorNames);
        if (jobs.isEmpty()) {
            return;
        }
        List<FragmentRef> fragments = collectFragments(context.getSource());
        if (fragments.isEmpty()) {
            return;
        }
        ExecutorService executor = threads > 1 ? Executors.newFixedThreadPool(threads) : null;
        long processedCount = 0;
        try {
            if (executor == null) {
                ChunkResult result = new DescriptorWorker(fragments, jobs, overwriteExisting, context).call();
                applyResult(result, context);
                processedCount += result.getProcessed();
                logProgress(processedCount, fragments.size());
            } else {
                CompletionService<ChunkResult> completionService = new ExecutorCompletionService<>(executor);
                List<Future<ChunkResult>> futures = new ArrayList<>();
                for (int start = 0; start < fragments.size(); start += DEFAULT_BATCH_SIZE) {
                    int end = Math.min(fragments.size(), start + DEFAULT_BATCH_SIZE);
                    List<FragmentRef> batch = fragments.subList(start, end);
                    futures.add(completionService.submit(
                            new DescriptorWorker(new ArrayList<>(batch), jobs, overwriteExisting, context)));
                }
                int completed = 0;
                while (completed < futures.size()) {
                    ChunkResult result = completionService.take().get();
                    applyResult(result, context);
                    completed++;
                    processedCount += result.getProcessed();
                    logProgress(processedCount, fragments.size());
                }
                executor.shutdown();
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
            }
        } finally {
            if (executor != null && !executor.isShutdown()) {
                executor.shutdownNow();
            }
        }
        for (DescriptorJob job : jobs) {
            context.addDescriptorTag(job.spec.canonicalName);
        }
    }

    private void applyResult(ChunkResult result, Context context) {
        for (FragmentAttribute attribute : result.getAttributes()) {
            context.addFragmentAttribute(attribute.reactionId, attribute.fragmentId, attribute.key, attribute.value);
        }
    }

    private static void logProgress(long processed, long total) {
        if (total <= 0) {
            System.out.println("[RawProcessor] processed " + processed + " fragments");
            return;
        }
        long percent = Math.min(100, (processed * 100) / total);
        System.out.println("[RawProcessor] " + percent + "% processed (" + processed + " fragments)");
    }

    private static List<DescriptorJob> buildJobs(List<String> names) {
        List<DescriptorJob> jobs = new ArrayList<>();
        LinkedHashSet<String> unique = new LinkedHashSet<>();
        for (String name : names) {
            if (name == null) {
                continue;
            }
            String normalized = name.trim();
            if (normalized.isEmpty() || !unique.add(normalized.toLowerCase(Locale.ROOT))) {
                continue;
            }
            DescriptorSpec spec = DescriptorSpec.forName(normalized);
            if (spec == null) {
                throw new IllegalArgumentException("Descriptor not supported for augmentation: " + name);
            }
            jobs.add(new DescriptorJob(spec));
        }
        return jobs;
    }

    private static List<FragmentRef> collectFragments(RawSynthonSpace source) {
        Map<String, FragmentRef> refs = new LinkedHashMap<>();
        source.getReactions().forEach((reactionId, data) -> data.getRawFragmentSets().forEach((idx, fragments) -> {
            for (RawSynthon frag : fragments) {
                String key = reactionId + "\u0000" + frag.getFragmentId();
                refs.putIfAbsent(key, new FragmentRef(reactionId, frag));
            }
        }));
        return new ArrayList<>(refs.values());
    }

    private static final class DescriptorWorker implements Callable<ChunkResult> {
        private final List<FragmentRef> fragments;
        private final List<DescriptorJob> jobs;
        private final boolean overwrite;
        private final Context context;
        private final IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();

        private DescriptorWorker(List<FragmentRef> fragments,
                                 List<DescriptorJob> jobs,
                                 boolean overwrite,
                                 Context context) {
            this.fragments = fragments;
            this.jobs = jobs;
            this.overwrite = overwrite;
            this.context = context;
        }

        @Override
        public ChunkResult call() throws Exception {
            List<FragmentAttribute> updates = new ArrayList<>();
            long processed = 0;
            for (FragmentRef ref : fragments) {
                StereoMolecule molecule = new StereoMolecule();
                parser.parse(molecule, ref.synthon.getIdcode());
                molecule.ensureHelperArrays(Molecule.cHelperCIP);
                for (DescriptorJob job : jobs) {
                    String key = ATTRIBUTE_PREFIX + job.spec.canonicalName;
                    if (!overwrite && context.hasFragmentAttribute(ref.reactionId, ref.synthon.getFragmentId(), key)) {
                        continue;
                    }
                    String value = job.compute(molecule);
                    if (value != null) {
                        updates.add(new FragmentAttribute(ref.reactionId, ref.synthon.getFragmentId(), key, value));
                    }
                }
                processed++;
            }
            return new ChunkResult(processed, updates);
        }
    }

    private static final class ChunkResult {
        private final long processed;
        private final List<FragmentAttribute> attributes;

        private ChunkResult(long processed, List<FragmentAttribute> attributes) {
            this.processed = processed;
            this.attributes = attributes;
        }

        long getProcessed() {
            return processed;
        }

        List<FragmentAttribute> getAttributes() {
            return attributes;
        }
    }

    private static final class FragmentAttribute {
        private final String reactionId;
        private final String fragmentId;
        private final String key;
        private final String value;

        private FragmentAttribute(String reactionId, String fragmentId, String key, String value) {
            this.reactionId = reactionId;
            this.fragmentId = fragmentId;
            this.key = key;
            this.value = value;
        }
    }

    private static final class FragmentRef {
        private final String reactionId;
        private final RawSynthon synthon;

        private FragmentRef(String reactionId, RawSynthon synthon) {
            this.reactionId = reactionId;
            this.synthon = synthon;
        }
    }

    private static final class DescriptorJob {
        private final DescriptorSpec spec;
        private final ThreadLocal<DescriptorHandler<?, StereoMolecule>> handlers;

        private DescriptorJob(DescriptorSpec spec) {
            this.spec = spec;
            this.handlers = ThreadLocal.withInitial(() -> spec.prototype.getThreadSafeCopy());
        }

        private String compute(StereoMolecule molecule) {
            DescriptorHandler handler = handlers.get();
            Object descriptor = handler.createDescriptor(molecule);
            if (handler.calculationFailed(descriptor)) {
                return null;
            }
            return spec.serializer.serialize(descriptor);
        }
    }

    private interface DescriptorSerializer {
        String serialize(Object descriptor);
    }

    private enum DescriptorSpec {
        FRAG_FP("FragFp", new DescriptorHandlerLongFFP512(), DescriptorSpec::encodeLongArray),
        PATH_FP("PathFp", new DescriptorHandlerLongPFP512(), DescriptorSpec::encodeLongArray),
        MODE_PFP("mode_pfp", new DescriptorHandlerLongFFP1024_plus("pfp"), DescriptorSpec::encodeLongArray),
        MODE_FFP("mode_ffp", new DescriptorHandlerLongFFP1024_plus("ffp"), DescriptorSpec::encodeLongArray),
        MODE_PPCORE("mode_ppcore", new DescriptorHandlerPPCore(), DescriptorSpec::encodeLongArray),
        SKELSPHERES("SkelSpheres", new DescriptorHandlerBinarySkelSpheres(), DescriptorSpec::encodeIntArray);

        private final String canonicalName;
        private final DescriptorHandler<?, StereoMolecule> prototype;
        private final DescriptorSerializer serializer;

        DescriptorSpec(String canonicalName,
                       DescriptorHandler<?, StereoMolecule> prototype,
                       DescriptorSerializer serializer) {
            this.canonicalName = canonicalName;
            this.prototype = prototype;
            this.serializer = serializer;
        }

        private static DescriptorSpec forName(String raw) {
            if (raw == null) {
                return null;
            }
            String token = raw.trim().toLowerCase(Locale.ROOT);
            switch (token) {
                case "fragfp":
                    return FRAG_FP;
                case "pathfp":
                    return PATH_FP;
                case "mode_pfp":
                case "pfp":
                    return MODE_PFP;
                case "mode_ffp":
                case "ffp":
                    return MODE_FFP;
                case "mode_ppcore":
                case "ppcore":
                case "ppc":
                    return MODE_PPCORE;
                case "skelspheres":
                case "skel":
                    return SKELSPHERES;
                default:
                    return null;
            }
        }

        private static String encodeLongArray(Object descriptor) {
            long[] data = (long[]) descriptor;
            if (data == null) {
                return null;
            }
            byte[] bytes = new byte[data.length * Long.BYTES];
            ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.BIG_ENDIAN);
            for (long value : data) {
                buffer.putLong(value);
            }
            return Base64.getEncoder().encodeToString(bytes);
        }

        private static String encodeIntArray(Object descriptor) {
            int[] data = (int[]) descriptor;
            if (data == null) {
                return null;
            }
            byte[] bytes = new byte[data.length * Integer.BYTES];
            ByteBuffer buffer = ByteBuffer.wrap(bytes).order(ByteOrder.BIG_ENDIAN);
            for (int value : data) {
                buffer.putInt(value);
            }
            return Base64.getEncoder().encodeToString(bytes);
        }
    }
}
