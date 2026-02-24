package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.io.RawSynthonSpaceImporter;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser3;
import org.apache.commons.cli.*;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.zip.GZIPOutputStream;

/**
 * CLI to parse vendor TSV files (Enamine format) into RawSynthonSpace JSON payloads.
 */
public class RawSynthonSpaceImportCLI {

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("RawSynthonSpaceImportCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        String format = cmd.getOptionValue("format", "enamine").toLowerCase(Locale.ROOT);
        String spaceName = cmd.getOptionValue("spaceName");
        String mode = cmd.getOptionValue("mode", "FragFp");
        String rawOut = cmd.getOptionValue("rawOut");
        String synthonOut = cmd.getOptionValue("synthonOut");
        String similarityOut = cmd.getOptionValue("similarityOut");
        int similarityThreads = Integer.parseInt(cmd.getOptionValue("similarityThreads", "4"));
        int threads = Integer.parseInt(cmd.getOptionValue("threads", "4"));
        LinkedHashSet<String> descriptorTags = parseDescriptorTags(cmd.getOptionValues("descriptorTag"));

        if (!descriptorTags.isEmpty() && rawOut == null) {
            throw new IllegalArgumentException("--descriptorTag can only be used when --rawOut is specified");
        }
        if (descriptorTags.size() > 1 && (synthonOut != null || similarityOut != null)) {
            throw new IllegalArgumentException("Multiple --descriptorTag values are only supported when producing raw output only");
        }

        boolean needSynthonSpace = synthonOut != null || similarityOut != null;
        RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
        RawSynthonSpaceImporter.Result result;
        if ("enamine".equals(format)) {
            String inputPath = cmd.getOptionValue("input");
            if (inputPath == null) {
                throw new IllegalArgumentException("--input is required for Enamine format");
            }
            int maxSets = Integer.parseInt(cmd.getOptionValue("maxSets", "3"));
            SynthonSpaceParser3.EnamineOptions.Builder builder = SynthonSpaceParser3.EnamineOptions.builder()
                    .input(Path.of(inputPath))
                    .spaceName(spaceName)
                    .mode(mode)
                    .threads(threads)
                    .maxSynthonSets(maxSets)
                    .buildSynthonSpace(needSynthonSpace);
            descriptorTags.forEach(builder::addDescriptorTag);
            result = importer.importEnamine(builder.build());
        } else if ("csv".equals(format)) {
            String directory = cmd.getOptionValue("directory");
            if (directory == null) {
                throw new IllegalArgumentException("--directory is required for CSV format");
            }
            SynthonSpaceParser3.CsvDirectoryOptions.Builder builder = SynthonSpaceParser3.CsvDirectoryOptions.builder()
                    .directory(Path.of(directory))
                    .spaceName(spaceName)
                    .mode(mode)
                    .smilesColumn(cmd.getOptionValue("smilesColumn", "SMILES"))
                    .idColumn(cmd.getOptionValue("idColumn", "bb1_parent_id"))
                    .priceColumn(cmd.getOptionValue("priceColumn", "Price"))
                    .priceAttributeKey(cmd.getOptionValue("priceAttribute", "price"))
                    .defaultSynthonSet(Integer.parseInt(cmd.getOptionValue("synthonDefaultSet", "0")))
                    .threads(threads)
                    .buildSynthonSpace(needSynthonSpace);
            if (cmd.hasOption("synthonSetColumn")) {
                builder.synthonSetColumn(cmd.getOptionValue("synthonSetColumn"));
            }
            descriptorTags.forEach(builder::addDescriptorTag);
            result = importer.importCsvDirectory(builder.build());
        } else {
            throw new IllegalArgumentException("Unsupported format: " + format);
        }
        RawSynthonSpace rawSpace = result.getRawSpace();
        HyperspaceIOUtils.saveRawSynthonSpace(rawSpace, rawOut);
        System.out.println("Raw synthon space written to: " + rawOut);

        if (needSynthonSpace) {
            SynthonSpace synthonSpace = result.getSynthonSpace();
            if (synthonSpace == null) {
                throw new IllegalStateException("SynthonSpace output requested but parser did not create one");
            }
            if (synthonOut != null) {
                HyperspaceIOUtils.saveSynthonSpace(synthonSpace, synthonOut);
                System.out.println("SynthonSpace written to: " + synthonOut);
            }
            if (similarityOut != null) {
                SynthonSimilaritySpace3 similarity = new SynthonSimilaritySpace3(synthonSpace);
                similarity.initFastSimilaritySearchers(similarityThreads, (progress) -> {
                });
                writeSimilaritySpace(similarity, similarityOut);
                System.out.println("Similarity space written to: " + similarityOut);
            }
        }
    }

    private static void writeSimilaritySpace(SynthonSimilaritySpace3 similarity, String output) {
        try (ObjectOutputStream out = new ObjectOutputStream(
                new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(output))))) {
            out.writeObject(similarity);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to store similarity space", e);
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("format").hasArg()
                .desc("Input format: enamine (default) or csv").build());
        options.addOption(Option.builder().longOpt("input").hasArg()
                .desc("Path to vendor TSV file (Enamine format)").build());
        options.addOption(Option.builder().longOpt("directory").hasArg()
                .desc("Directory containing per-reaction CSV files (CSV format)").build());
        options.addOption(Option.builder().longOpt("spaceName").hasArg().required(true)
                .desc("Logical name of the synthon space").build());
        options.addOption(Option.builder().longOpt("mode").hasArg()
                .desc("Descriptor mode (FragFp, PathFp, mode_pfp, mode_ffp, mode_ppcore)").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Worker threads for parsing (default 4)").build());
        options.addOption(Option.builder().longOpt("maxSets").hasArg()
                .desc("Max number of synthon sets per reaction to keep (default 3)").build());
        options.addOption(Option.builder().longOpt("descriptorTag").hasArg()
                .desc("Adds a descriptor short name tag to the raw output metadata (repeatable; raw-only)").build());
        options.addOption(Option.builder().longOpt("smilesColumn").hasArg()
                .desc("CSV column that stores SMILES (default SMILES)").build());
        options.addOption(Option.builder().longOpt("idColumn").hasArg()
                .desc("CSV column that stores fragment identifiers (default bb1_parent_id)").build());
        options.addOption(Option.builder().longOpt("priceColumn").hasArg()
                .desc("CSV column containing price metadata (default Price)").build());
        options.addOption(Option.builder().longOpt("priceAttribute").hasArg()
                .desc("Attribute key used to store price metadata (default price)").build());
        options.addOption(Option.builder().longOpt("synthonSetColumn").hasArg()
                .desc("CSV column containing synthon set indices (optional)").build());
        options.addOption(Option.builder().longOpt("synthonDefaultSet").hasArg()
                .desc("Default synthon set index if column missing/empty (default 0)").build());
        options.addOption(Option.builder().longOpt("rawOut").hasArg().required(true)
                .desc("Output JSON file for RawSynthonSpace").build());
        options.addOption(Option.builder().longOpt("synthonOut").hasArg()
                .desc("Optional .data file for immediate SynthonSpace serialization").build());
        options.addOption(Option.builder().longOpt("similarityOut").hasArg()
                .desc("Optional .data file for SynthonSimilaritySpace3").build());
        options.addOption(Option.builder().longOpt("similarityThreads").hasArg()
                .desc("Worker threads for similarity construction (default 4)").build());
        return options;
    }

    private static LinkedHashSet<String> parseDescriptorTags(String[] rawValues) {
        LinkedHashSet<String> tags = new LinkedHashSet<>();
        if (rawValues == null) {
            return tags;
        }
        Arrays.stream(rawValues)
                .filter(value -> value != null && !value.isBlank())
                .flatMap(value -> Arrays.stream(value.split(",")))
                .map(String::trim)
                .filter(token -> !token.isEmpty())
                .forEach(tags::add);
        return tags;
    }
}
