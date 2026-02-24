package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceAssembler;
import org.apache.commons.cli.*;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.zip.GZIPOutputStream;

/**
 * CLI to convert RawSynthonSpace JSON dumps into fully initialized SynthonSpace binaries.
 */
public class RawSynthonSpaceBuildCLI {

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("RawSynthonSpaceBuildCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        String rawPath = cmd.getOptionValue("rawIn");
        String synthonOut = cmd.getOptionValue("synthonOut");
        String descriptor = cmd.getOptionValue("descriptor");
        String bitsArg = cmd.getOptionValue("bits");
        boolean skipValidation = cmd.hasOption("skipValidation");
        String similarityOut = cmd.getOptionValue("similarityOut");
        int similarityThreads = Integer.parseInt(cmd.getOptionValue("similarityThreads", "4"));

        RawSynthonSpace raw = HyperspaceIOUtils.loadRawSynthonSpace(rawPath);
        RawSynthonSpaceAssembler.BuildOptions.Builder builder = RawSynthonSpaceAssembler.BuildOptions.builder();
        if (descriptor == null || descriptor.isBlank()) {
            String available = raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS);
            String suggestion = available == null ? "unknown" : available;
            throw new IllegalArgumentException("Please supply --descriptor (available tags: " + suggestion + ")");
        }
        builder.descriptorShortName(descriptor.trim());
        if (bitsArg != null) {
            builder.descriptorBits(Integer.parseInt(bitsArg));
        }
        if (skipValidation) {
            builder.validateReactions(false);
        }

        SynthonSpace space = RawSynthonSpaceAssembler.buildSynthonSpace(raw, builder.build());
        HyperspaceIOUtils.saveSynthonSpace(space, synthonOut);
        System.out.println("Synthon space written to: " + synthonOut);

        if (similarityOut != null) {
            SynthonSimilaritySpace3 similaritySpace = new SynthonSimilaritySpace3(space);
            similaritySpace.initFastSimilaritySearchers(similarityThreads, (progress) -> {
            });
            writeSimilaritySpace(similaritySpace, similarityOut);
            System.out.println("Similarity space written to: " + similarityOut);
        }
    }

    private static void writeSimilaritySpace(SynthonSimilaritySpace3 similaritySpace, String output) {
        try (ObjectOutputStream out = new ObjectOutputStream(
                new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(output))))) {
            out.writeObject(similaritySpace);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to write similarity space to " + output, e);
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawIn").hasArg().required(true)
                .desc("Path to RawSynthonSpace JSON file").build());
        options.addOption(Option.builder().longOpt("synthonOut").hasArg().required(true)
                .desc("Output file (.data) for the reconstructed SynthonSpace").build());
        options.addOption(Option.builder().longOpt("descriptor").hasArg().required(true)
                .desc("Descriptor short name to use when materializing the SynthonSpace (e.g., FragFp, mode_pfp)").build());
        options.addOption(Option.builder().longOpt("bits").hasArg()
                .desc("Optional descriptor length in bits (defaults to metadata)").build());
        options.addOption(Option.builder().longOpt("skipValidation")
                .desc("Skip connector validation; trust RawSynthonSpace content.").build());
        options.addOption(Option.builder().longOpt("similarityOut").hasArg()
                .desc("Optional output for SynthonSimilaritySpace3 (.data)").build());
        options.addOption(Option.builder().longOpt("similarityThreads").hasArg()
                .desc("Worker threads for similarity initialization (default 4)").build());
        return options;
    }
}
