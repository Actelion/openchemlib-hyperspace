package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.downsampling.*;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceProcessor;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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
        double minSimilarity = Double.parseDouble(cmd.getOptionValue("minSimilarity", "0.7"));
        long seed = Long.parseLong(cmd.getOptionValue("seed", "13"));
        boolean enforceConnectors = !cmd.hasOption("allowConnectorMixing");
        String rawOut = cmd.getOptionValue("rawOut");

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

        SynthonDownsamplingRequest request = SynthonDownsamplingRequest.builder()
                .withMaxCenters(maxCenters)
                .withMinSimilarity(minSimilarity)
                .withRandomSeed(seed)
                .enforceConnectorEquivalence(enforceConnectors)
                .build();

        SynthonDownsamplingResult result;
        DownsampledSynthonSpace downsampledSpace;
        String algorithmName;

        if (rawSource != null) {
            RawSynthonDownsampler downsampler = new SkelSpheresKCentersRawDownsampler();
            RawSynthonDownsamplingOrchestrator orchestrator = new RawSynthonDownsamplingOrchestrator();
            result = orchestrator.downsample(rawSource, downsampler, request);
            downsampledSpace = result.toDownsampledSpace();
            algorithmName = downsampler.getName();
        } else {
            SynthonDownsampler downsampler = new SkelSpheresKCentersDownsampler();
            SynthonDownsamplingOrchestrator orchestrator = new SynthonDownsamplingOrchestrator();
            result = orchestrator.downsample(space, downsampler, request);
            downsampledSpace = result.toDownsampledSpace();
            algorithmName = downsampler.getName();
        }
        if (rawOut != null) {
            RawSynthonSpace rawSpace;
            if (rawSource != null) {
                try {
                    rawSpace = attachDownsampledSets(rawSource, downsampledSpace, algorithmName, request);
                } catch (Exception e) {
                    throw new IllegalStateException("Failed to merge downsampled sets into raw space", e);
                }
            } else {
                rawSpace = RawSynthonSpace.builder(new java.io.File(spaceIn).getName())
                        .withFullSynthonSpace(space)
                        .withDownsampledSpace(downsampledSpace, algorithmName, request)
                        .build();
            }
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
        System.out.println("Downsampled representatives stored in: " + rawOut);
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("spaceIn").hasArg()
                .desc("Path to serialized SynthonSpace (.synthonspace.gz)").build());
        options.addOption(Option.builder().longOpt("rawIn").hasArg()
                .desc("Path to RawSynthonSpace JSON (alternative to --spaceIn)").build());
        options.addOption(Option.builder().longOpt("maxCenters").hasArg()
                .desc("Maximum number of centers per synthon set (0 = unlimited)").build());
        options.addOption(Option.builder().longOpt("minSimilarity").hasArg()
                .desc("Minimum SkelSpheres similarity required to reuse an existing center [0-1]").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Random seed for shuffling synthons").build());
        options.addOption(Option.builder().longOpt("allowConnectorMixing")
                .desc("Allow synthons with different connector patterns to join the same center").build());
        options.addOption(Option.builder().longOpt("rawOut").hasArg().required(true)
                .desc("Output RawSynthonSpace file ( .rawspace or .rawspace.gz ) containing downsampled sets").build());
        return options;
    }

    private static RawSynthonSpace attachDownsampledSets(RawSynthonSpace source,
                                                         DownsampledSynthonSpace downsampledSpace,
                                                         String algorithm,
                                                         SynthonDownsamplingRequest request) throws Exception {
        return RawSynthonSpaceProcessor.process(source,
                List.of(new DownsampledSetsTask(downsampledSpace, algorithm, request)));
    }

    private static final class DownsampledSetsTask implements RawSynthonSpaceProcessor.Task {
        private final DownsampledSynthonSpace downsampledSpace;
        private final String algorithm;
        private final SynthonDownsamplingRequest request;

        private DownsampledSetsTask(DownsampledSynthonSpace downsampledSpace,
                                    String algorithm,
                                    SynthonDownsamplingRequest request) {
            this.downsampledSpace = downsampledSpace;
            this.algorithm = algorithm;
            this.request = request;
        }

        @Override
        public void apply(RawSynthonSpaceProcessor.Context context) {
            context.getBuilder().withDownsamplingMetadata(algorithm, request);
            downsampledSpace.getDownsampledSets().forEach((rxnId, fragMap) ->
                    fragMap.forEach((idx, synthons) ->
                            context.getBuilder().addRawDownsampledFragments(rxnId, idx, toRawSynthons(synthons))));
        }

        private List<RawSynthon> toRawSynthons(List<SynthonSpace.FragId> synthons) {
            List<RawSynthon> raws = new ArrayList<>(synthons.size());
            for (SynthonSpace.FragId frag : synthons) {
                raws.add(RawSynthon.fromFragId(frag));
            }
            return raws;
        }
    }
}
