package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsExporter;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsOptions;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsReport;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.nio.file.Path;
import java.util.Locale;

public final class RawSynthonSpaceStatsCLI {
    private RawSynthonSpaceStatsCLI() {
    }

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd = parse(options, args);
        if (cmd.hasOption("help")) {
            new HelpFormatter().printHelp("RawSynthonSpaceStatsCLI", options, true);
            return;
        }

        String rawIn = cmd.getOptionValue("rawIn");
        String outDir = cmd.getOptionValue("outDir");
        if (isBlank(rawIn) || isBlank(outDir)) {
            new HelpFormatter().printHelp("RawSynthonSpaceStatsCLI", options, true);
            throw new IllegalArgumentException("Both --rawIn and --outDir are required");
        }

        int examplesPerReaction = parseInt(cmd.getOptionValue("examplesPerReaction"), 10);
        int productSamplesPerReaction = parseInt(cmd.getOptionValue("productSamplesPerReaction"), 100);
        long seed = parseLong(cmd.getOptionValue("seed"), 13L);
        int threads = parseInt(cmd.getOptionValue("threads"), Runtime.getRuntime().availableProcessors());
        int maxReactions = parseInt(cmd.getOptionValue("maxReactions"), 0);

        System.out.println("Loading rawspace: " + rawIn);
        RawSynthonSpace rawSpace = HyperspaceIOUtils.loadRawSynthonSpace(rawIn);
        System.out.printf(Locale.ROOT, "Loaded rawspace name=%s reactions=%d%n",
                rawSpace.getName(), rawSpace.getReactions().size());

        RawSynthonSpaceStatsOptions reportOptions = RawSynthonSpaceStatsOptions.builder()
                .examplesPerReaction(examplesPerReaction)
                .productSamplesPerReaction(productSamplesPerReaction)
                .seed(seed)
                .threads(threads)
                .maxReactions(maxReactions)
                .progressListener((completed, total, reactionId) -> {
                    if (completed == total || completed % Math.max(1, total / 20) == 0) {
                        System.out.printf(Locale.ROOT, "Stats progress: %d / %d (%d%%)%n",
                                completed, total, Math.round(100.0 * completed / Math.max(1, total)));
                    }
                })
                .build();

        RawSynthonSpaceStatsReport report = RawSynthonSpaceStatsExporter.analyze(rawSpace, reportOptions);
        Path outputDirectory = Path.of(outDir).toAbsolutePath().normalize();
        report.writeToDirectory(outputDirectory);
        System.out.printf(Locale.ROOT,
                "Wrote rawspace stats to %s (reactions=%d examples=%d)%n",
                outputDirectory,
                report.getReactionStats().size(),
                report.getExampleProducts().size());
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawIn").hasArg()
                .desc("Input RawSynthonSpace file (.rawspace or .rawspace.gz)").build());
        options.addOption(Option.builder().longOpt("outDir").hasArg()
                .desc("Output directory for the stats report bundle").build());
        options.addOption(Option.builder().longOpt("examplesPerReaction").hasArg()
                .desc("Random assembled examples per reaction (default 10)").build());
        options.addOption(Option.builder().longOpt("productSamplesPerReaction").hasArg()
                .desc("Random assembled products sampled per reaction for product-property stats (default 100)").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Random seed for deterministic sampling (default 13)").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Worker threads (default available processors)").build());
        options.addOption(Option.builder().longOpt("maxReactions").hasArg()
                .desc("Limit number of sorted reactions to process; 0 means all (default 0)").build());
        options.addOption(Option.builder().longOpt("help")
                .desc("Show this help message").build());
        return options;
    }

    private static CommandLine parse(Options options, String[] args) {
        try {
            return new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }
    }

    private static int parseInt(String value, int defaultValue) {
        if (isBlank(value)) {
            return defaultValue;
        }
        return Integer.parseInt(value.trim());
    }

    private static long parseLong(String value, long defaultValue) {
        if (isBlank(value)) {
            return defaultValue;
        }
        return Long.parseLong(value.trim());
    }

    private static boolean isBlank(String value) {
        return value == null || value.isBlank();
    }
}
