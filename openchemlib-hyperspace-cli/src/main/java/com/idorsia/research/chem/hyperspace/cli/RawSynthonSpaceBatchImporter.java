package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceImporter;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser3;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

/**
 * Batch driver that reads TSV manifests describing raw synthon spaces
 * (vendor format files) and materializes all of them into RawSynthonSpace
 * JSON dumps.
 */
public final class RawSynthonSpaceBatchImporter {

    private RawSynthonSpaceBatchImporter() {
    }

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd = parse(options, args);

        Path root = Path.of(cmd.getOptionValue("root"));
        Path enamineList = cmd.hasOption("enamine") ? root.resolve(cmd.getOptionValue("enamine")) : null;
        Path synpleList = cmd.hasOption("synple") ? root.resolve(cmd.getOptionValue("synple")) : null;
        Path outputDir = root.resolve(cmd.getOptionValue("out"));
        Files.createDirectories(outputDir);

        int importerThreads = Integer.parseInt(cmd.getOptionValue("threads", "4"));
        String mode = cmd.getOptionValue("mode");
        int maxSets = Integer.parseInt(cmd.getOptionValue("maxSets", "3"));

        RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
        List<ManifestEntry> entries = new ArrayList<>();
        if (enamineList != null) {
            entries.addAll(readManifest(enamineList, ManifestEntry.Format.ENAMINE));
        }
        if (synpleList != null) {
            entries.addAll(readManifest(synpleList, ManifestEntry.Format.SYNPLE));
        }
        if (entries.isEmpty()) {
            System.out.println("No entries found; nothing to do.");
            return;
        }

        for (ManifestEntry entry : entries) {
            System.out.printf("Importing %s from %s (%s)%n", entry.spaceName, entry.input, entry.format);
            Path inputPath = root.resolve(entry.input).normalize();
            if (!Files.exists(inputPath)) {
                System.err.printf("  Skipping %s: input file does not exist (%s)%n", entry.spaceName, inputPath);
                continue;
            }
            try {
                RawSynthonSpace rawSpace = importSingle(importer, entry, inputPath, mode, importerThreads, maxSets);
                Path outputFile = outputDir.resolve(entry.spaceName + ".json.gz");
                HyperspaceIOUtils.saveRawSynthonSpace(rawSpace, outputFile.toString());
                System.out.printf("  -> stored %s%n", outputFile);
            } catch (Exception ex) {
                System.err.printf("  Failed to import %s: %s%n", entry.spaceName, ex.getMessage());
                ex.printStackTrace(System.err);
            }
        }
    }

    private static RawSynthonSpace importSingle(RawSynthonSpaceImporter importer,
                                                ManifestEntry entry,
                                                Path input,
                                                String mode,
                                                int threads,
                                                int maxSets) throws Exception {
        switch (entry.format) {
            case ENAMINE -> {
                SynthonSpaceParser3.EnamineOptions.Builder builder = SynthonSpaceParser3.EnamineOptions.builder()
                        .input(input)
                        .spaceName(entry.spaceName)
                        .threads(threads)
                        .maxSynthonSets(maxSets)
                        .buildSynthonSpace(false);
                if (mode != null && !mode.isBlank()) {
                    builder.mode(mode);
                }
                entry.descriptorTags.forEach(builder::addDescriptorTag);
                return importer.importEnamine(builder.build()).getRawSpace();
            }
            case SYNPLE -> {
                SynthonSpaceParser3.EnamineOptions.Builder builder = SynthonSpaceParser3.EnamineOptions.builder()
                        .input(input)
                        .spaceName(entry.spaceName)
                        .threads(threads)
                        .maxSynthonSets(maxSets)
                        .buildSynthonSpace(false);
                if (mode != null && !mode.isBlank()) {
                    builder.mode(mode);
                }
                entry.descriptorTags.forEach(builder::addDescriptorTag);
                return importer.importEnamine(builder.build()).getRawSpace();
            }
            default -> throw new IllegalStateException("Unsupported format " + entry.format);
        }
    }

    private static List<ManifestEntry> readManifest(Path manifest, ManifestEntry.Format format) throws IOException {
        List<ManifestEntry> entries = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(manifest, StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                String[] tokens = line.split("\t");
                if (tokens.length < 2) {
                    continue;
                }
                String spaceName = tokens[0].trim();
                String input = tokens[1].trim();
                List<String> descriptorTags = new ArrayList<>();
                if (tokens.length >= 3 && !tokens[2].isBlank()) {
                    descriptorTags.add(tokens[2].trim());
                }
                entries.add(new ManifestEntry(spaceName, input, descriptorTags, format));
            }
        }
        return entries;
    }

    private static CommandLine parse(Options options, String[] args) throws ParseException {
        try {
            return new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("RawSynthonSpaceBatchImporter", options, true);
            throw e;
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("root").hasArg().required(true)
                .desc("Root directory that contains TSVs and vendor files").build());
        options.addOption(Option.builder().longOpt("enamine").hasArg()
                .desc("TSV manifest listing Enamine-format sources").build());
        options.addOption(Option.builder().longOpt("synple").hasArg()
                .desc("TSV manifest listing Synple CSV directories").build());
        options.addOption(Option.builder().longOpt("out").hasArg().required(true)
                .desc("Output directory for RawSynthonSpace JSON(.gz) files").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Worker threads per individual import (default 4)").build());
        options.addOption(Option.builder().longOpt("mode").hasArg()
                .desc("Descriptor mode passed to importers (omit to use importer defaults)").build());
        options.addOption(Option.builder().longOpt("maxSets").hasArg()
                .desc("Max synthon sets per reaction for Enamine imports (default 3)").build());
        return options;
    }

    private static final class ManifestEntry {
        private final String spaceName;
        private final String input;
        private final List<String> descriptorTags;
        private final Format format;

        private ManifestEntry(String spaceName, String input, List<String> descriptorTags, Format format) {
            this.spaceName = Objects.requireNonNull(spaceName, "spaceName");
            this.input = Objects.requireNonNull(input, "input");
            this.descriptorTags = descriptorTags == null ? List.of() : List.copyOf(descriptorTags);
            this.format = Objects.requireNonNull(format, "format");
        }

        private enum Format {
            ENAMINE,
            SYNPLE
        }
    }
}
