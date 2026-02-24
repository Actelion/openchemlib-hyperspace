package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.io.DescriptorAugmentationTask;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceProcessor;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.apache.commons.cli.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * CLI entry point to process RawSynthonSpace dumps (e.g., attach descriptors).
 */
public class RawSynthonSpaceProcessorCLI {

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("RawSynthonSpaceProcessorCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        String rawIn = cmd.getOptionValue("rawIn");
        String rawOut = cmd.getOptionValue("rawOut");
        boolean overwrite = cmd.hasOption("overwriteDescriptors");
        int threads = Integer.parseInt(cmd.getOptionValue("threads",
                Integer.toString(Math.max(1, Runtime.getRuntime().availableProcessors()))));

        List<RawSynthonSpaceProcessor.Task> tasks = new ArrayList<>();
        List<String> descriptorNames = parseList(cmd.getOptionValues("descriptor"));
        if (!descriptorNames.isEmpty()) {
            tasks.add(new DescriptorAugmentationTask(descriptorNames, threads, overwrite));
        }

        if (tasks.isEmpty()) {
            throw new IllegalArgumentException("No processing tasks selected. Specify --descriptor to attach descriptors.");
        }

        RawSynthonSpace raw = HyperspaceIOUtils.loadRawSynthonSpace(rawIn);
        RawSynthonSpace processed = RawSynthonSpaceProcessor.process(raw, tasks);
        HyperspaceIOUtils.saveRawSynthonSpace(processed, rawOut);
        System.out.println("Processed raw space written to: " + rawOut);
    }

    private static List<String> parseList(String[] rawValues) {
        List<String> tokens = new ArrayList<>();
        if (rawValues == null) {
            return tokens;
        }
        for (String value : rawValues) {
            if (value == null) {
                continue;
            }
            Arrays.stream(value.split(","))
                    .map(String::trim)
                    .filter(token -> !token.isEmpty())
                    .forEach(tokens::add);
        }
        return tokens;
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawIn").hasArg().required(true)
                .desc("Input RawSynthonSpace JSON file").build());
        options.addOption(Option.builder().longOpt("rawOut").hasArg().required(true)
                .desc("Output RawSynthonSpace JSON file after processing").build());
        options.addOption(Option.builder().longOpt("descriptor").hasArg()
                .desc("Descriptor short name(s) to attach (comma-separated or repeatable)").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Worker threads for descriptor computation (default = available processors)").build());
        options.addOption(Option.builder().longOpt("overwriteDescriptors")
                .desc("Overwrite descriptor attributes if they already exist on fragments").build());
        return options;
    }
}
