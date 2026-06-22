package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.*;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Lightweight CLI to trigger synthon space downsampling from the command line.
 */
public class SynthonSpaceDownsamplingCLI {

    public static void main(String[] args) {
        Options options = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("SynthonSpaceDownsamplingCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        String spaceIn = cmd.getOptionValue("spaceIn");
        String rawIn = cmd.getOptionValue("rawIn");
        if ((spaceIn == null && rawIn == null) || (spaceIn != null && rawIn != null)) {
            throw new IllegalArgumentException("Provide either --spaceIn or --rawIn");
        }
        int maxCenters = Integer.parseInt(cmd.getOptionValue("maxCenters", "0"));
        double sizeCapScale = Double.parseDouble(cmd.getOptionValue("sizeCapScale", "0"));
        double sizeCapOffset = Double.parseDouble(cmd.getOptionValue("sizeCapOffset", "0"));
        double minSimilarity = Double.parseDouble(cmd.getOptionValue("minSimilarity", "0.7"));
        long seed = Long.parseLong(cmd.getOptionValue("seed", "13"));
        int threads = Integer.parseInt(cmd.getOptionValue("threads", "1"));
        int setProgressIntervalSeconds = Integer.parseInt(cmd.getOptionValue("setProgressIntervalSeconds", "60"));
        boolean enforceConnectors = !cmd.hasOption("allowConnectorMixing");
        String rawOut = cmd.getOptionValue("rawOut");
        if (threads <= 0) {
            throw new IllegalArgumentException("--threads must be positive");
        }
        if (setProgressIntervalSeconds < 0) {
            throw new IllegalArgumentException("--setProgressIntervalSeconds must be >= 0");
        }

        SynthonSpace space = null;
        RawSynthonSpace rawSource = null;
        try {
            if (spaceIn != null) {
                space = HyperspaceIOUtils.loadSynthonSpace(spaceIn);
            } else {
                rawSource = HyperspaceIOUtils.loadRawSynthonSpace(rawIn);
            }
        } catch (IOException e) {
            throw new IllegalStateException("Failed to load input space", e);
        }

        SynthonDownsamplingRequest.Builder requestBuilder = SynthonDownsamplingRequest.builder()
                .withMaxCenters(maxCenters)
                .withSizeCapScale(sizeCapScale)
                .withSizeCapOffset(sizeCapOffset)
                .withMinSimilarity(minSimilarity)
                .withRandomSeed(seed)
                .enforceConnectorEquivalence(enforceConnectors)
                .withProgressReportIntervalSeconds(setProgressIntervalSeconds);
        if (setProgressIntervalSeconds > 0) {
            requestBuilder.withProgressReporter(SynthonSpaceDownsamplingCLI::printSetProgress);
        }
        SynthonDownsamplingRequest request = requestBuilder.build();

        SynthonDownsamplingResult result;
        String algorithmName;

        if (rawSource != null) {
            RawSynthonDownsampler downsampler = new SkelSpheresKCentersRawDownsampler();
            RawSynthonDownsamplingOrchestrator orchestrator = new RawSynthonDownsamplingOrchestrator();
            result = orchestrator.downsample(rawSource, downsampler, request, threads);
            algorithmName = downsampler.getName();
        } else {
            SynthonDownsampler downsampler = new SkelSpheresKCentersDownsampler();
            SynthonDownsamplingOrchestrator orchestrator = new SynthonDownsamplingOrchestrator();
            result = orchestrator.downsample(space, downsampler, request, threads);
            algorithmName = downsampler.getName();
        }
        if (rawOut != null) {
            RawSynthonSpace rawSpace = buildDownsampledRawSpace(rawSource,
                    spaceIn == null ? rawIn : new java.io.File(spaceIn).getName(),
                    result,
                    algorithmName,
                    request);
            try {
                HyperspaceIOUtils.saveRawSynthonSpace(rawSpace, rawOut);
                System.out.println("Raw synthon space written to: " + rawOut);
            } catch (IOException e) {
                throw new IllegalStateException("Failed to store raw synthon space", e);
            }
        }

        System.out.println("Downsampling finished with algorithm: " + algorithmName);
        System.out.println("Original synthons: " + result.getTotalOriginalSynthons());
        System.out.println("Retained synthons: " + result.getTotalRetainedSynthons());
        System.out.println("Discarded synthons: " + result.getTotalDiscardedSynthons());
        System.out.println("Downsampled rawspace stored in: " + rawOut);
    }

    private static synchronized void printSetProgress(SynthonDownsamplingRequest.ProgressEvent event) {
        int pct = event.getTotalSynthons() == 0
                ? 100
                : (int) Math.round(event.getProcessedSynthons() * 100.0 / event.getTotalSynthons());
        String centerLimit = event.getEffectiveMaxCenters() <= 0
                ? "unlimited"
                : Integer.toString(event.getEffectiveMaxCenters());
        System.out.println(String.format(Locale.ROOT,
                "Downsampling set rxn=%s frag=%d processed=%,d/%,d (%d%%) centers=%,d/%s skipped=%,d elapsed=%s%s",
                event.getReactionId(),
                event.getFragmentIndex(),
                event.getProcessedSynthons(),
                event.getTotalSynthons(),
                pct,
                event.getCenterCount(),
                centerLimit,
                event.getSkippedDescriptors(),
                formatDuration(event.getElapsedMillis()),
                event.isCompleted() ? " done" : ""));
    }

    private static String formatDuration(long elapsedMillis) {
        long totalSeconds = Math.max(0L, elapsedMillis / 1000L);
        long hours = totalSeconds / 3600L;
        long minutes = (totalSeconds % 3600L) / 60L;
        long seconds = totalSeconds % 60L;
        if (hours > 0) {
            return String.format(Locale.ROOT, "%dh%02dm%02ds", hours, minutes, seconds);
        }
        if (minutes > 0) {
            return String.format(Locale.ROOT, "%dm%02ds", minutes, seconds);
        }
        return seconds + "s";
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("spaceIn").hasArg()
                .desc("Path to serialized SynthonSpace (.synthonspace.gz)").build());
        options.addOption(Option.builder().longOpt("rawIn").hasArg()
                .desc("Path to RawSynthonSpace JSON (alternative to --spaceIn)").build());
        options.addOption(Option.builder().longOpt("maxCenters").hasArg()
                .desc("Maximum number of centers per synthon set (0 = unlimited)").build());
        options.addOption(Option.builder().longOpt("sizeCapScale").hasArg()
                .desc("Scale factor for size-based cap: scale * sqrt(n) + offset (0 = disabled)").build());
        options.addOption(Option.builder().longOpt("sizeCapOffset").hasArg()
                .desc("Offset for size-based cap: scale * sqrt(n) + offset (0 = disabled)").build());
        options.addOption(Option.builder().longOpt("minSimilarity").hasArg()
                .desc("Minimum SkelSpheres similarity required to reuse an existing center [0-1]").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Random seed for shuffling synthons").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Worker threads for downsampling (default 1)").build());
        options.addOption(Option.builder().longOpt("setProgressIntervalSeconds").hasArg()
                .desc("Seconds between progress reports for running synthon sets (default 60, 0 = disabled)").build());
        options.addOption(Option.builder().longOpt("allowConnectorMixing")
                .desc("Allow synthons with different connector patterns to join the same center").build());
        options.addOption(Option.builder().longOpt("rawOut").hasArg().required(true)
                .desc("Output reduced RawSynthonSpace file ( .rawspace or .rawspace.gz )").build());
        return options;
    }

    private static RawSynthonSpace buildDownsampledRawSpace(RawSynthonSpace source,
                                                            String fallbackName,
                                                            SynthonDownsamplingResult result,
                                                            String algorithm,
                                                            SynthonDownsamplingRequest request) {
        String name = source != null ? source.getName() : fallbackName;
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(name == null ? "downsampled_rawspace" : name);
        if (source != null) {
            builder.version(source.getVersion());
            source.getMetadata().forEach(builder::putMetadata);
        }
        builder.withDownsamplingMetadata(algorithm, request);

        for (SynthonSetDownsamplingResult setResult : result.getSetResults()) {
            SynthonSpace.FragType type = setResult.getFragType();
            String reactionId = type.rxn_id;
            builder.addRawFragments(reactionId, type.frag, toRawSynthons(setResult.getRepresentatives()));
            copyRetainedFragmentAttributes(source, builder, reactionId, setResult.getRepresentatives());
        }

        if (source != null) {
            source.getReactions().forEach((reactionId, data) -> {
                builder.addExampleScaffolds(reactionId, data.getExampleScaffolds());
                data.getPartialAssemblies().forEach((missingIdx, assemblies) ->
                        builder.addPartialAssemblies(reactionId, missingIdx, assemblies));
                builder.addRepresentativeCompounds(reactionId, data.getRepresentativeCompounds());
                data.getDescriptors().forEach((key, value) -> builder.addReactionDescriptor(reactionId, key, value));
                data.getReactionMetadata().forEach((key, value) -> builder.addReactionMetadata(reactionId, key, value));
            });
        }
        return builder.build();
    }

    private static void copyRetainedFragmentAttributes(RawSynthonSpace source,
                                                       RawSynthonSpace.Builder builder,
                                                       String reactionId,
                                                       List<SynthonSpace.FragId> representatives) {
        if (source == null) {
            return;
        }
        RawSynthonSpace.ReactionData data = source.getReactions().get(reactionId);
        if (data == null) {
            return;
        }
        Map<String, Map<String, String>> attributes = data.getFragmentAttributes();
        for (SynthonSpace.FragId representative : representatives) {
            Map<String, String> fragmentAttributes = attributes.get(representative.fragment_id);
            if (fragmentAttributes == null) {
                continue;
            }
            fragmentAttributes.forEach((key, value) ->
                    builder.addFragmentAttribute(reactionId, representative.fragment_id, key, value));
        }
    }

    private static List<RawSynthon> toRawSynthons(List<SynthonSpace.FragId> synthons) {
        List<RawSynthon> raws = new ArrayList<>(synthons.size());
        for (SynthonSpace.FragId frag : synthons) {
            raws.add(RawSynthon.fromFragId(frag));
        }
        return raws;
    }
}
