package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * CLI for merging multiple RawSynthonSpace dumps into one output or
 * into two split outputs (2-synthon-set reactions and 3-synthon-set reactions).
 */
public final class SynthonSpaceMergeCLI {

    private SynthonSpaceMergeCLI() {
    }

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd = parse(options, args);

        List<Path> inputFiles = parseInputFiles(cmd);
        if (inputFiles.isEmpty()) {
            throw new IllegalArgumentException("No inputs provided. Use --rawIn and/or --manifest and/or --inputDirectory");
        }

        boolean splitBySetCount = cmd.hasOption("splitBySetCount");
        String out = cmd.getOptionValue("out");
        String out2s = cmd.getOptionValue("out2s");
        String out3s = cmd.getOptionValue("out3s");
        String synthonListOut = cmd.getOptionValue("synthonListOut");

        if (splitBySetCount) {
            if (isBlank(out2s) || isBlank(out3s)) {
                throw new IllegalArgumentException("Split mode requires both --out2s and --out3s");
            }
        } else {
            if (isBlank(out)) {
                throw new IllegalArgumentException("Single-output mode requires --out");
            }
        }

        String baseName = cmd.getOptionValue("spaceName", "merged_rawspace");
        String version = cmd.getOptionValue("version", "1.0");
        String metadataPrefix = cmd.getOptionValue("sourceMetadataPrefix", "source");
        DuplicateStrategy duplicateStrategy = parseDuplicateStrategy(cmd.getOptionValue("duplicateStrategy", "rename"));

        List<LoadedSpace> loaded = new ArrayList<>();
        for (int i = 0; i < inputFiles.size(); i++) {
            Path file = inputFiles.get(i);
            RawSynthonSpace raw = HyperspaceIOUtils.loadRawSynthonSpace(file.toString());
            String sourceName = deriveSourceName(raw, file);
            loaded.add(new LoadedSpace(i + 1, file, raw, sourceName));
            System.out.printf(Locale.ROOT, "Loaded %-40s reactions=%d%n", file, raw.getReactions().size());
        }

        if (!isBlank(synthonListOut)) {
            Path outPath = Path.of(synthonListOut).toAbsolutePath().normalize();
            Map<String, SynthonStats> synthonStats = collectSynthonStats(loaded);
            writeSynthonList(outPath, synthonStats);
            System.out.printf(Locale.ROOT, "Wrote synthon list: %s (uniqueSynthons=%d)%n", outPath, synthonStats.size());
        }

        MergeResult merge = mergeSpaces(loaded, duplicateStrategy);
        System.out.printf(Locale.ROOT,
                "Merged reactions=%d (duplicates=%d, renamed=%d)%n",
                merge.reactions.size(),
                merge.duplicateConflicts,
                merge.renamedReactions);

        if (!splitBySetCount) {
            RawSynthonSpace merged = buildMergedSpace(baseName, version, loaded, merge.reactions, metadataPrefix, "single");
            HyperspaceIOUtils.saveRawSynthonSpace(merged, out);
            System.out.printf(Locale.ROOT, "Wrote merged raw space: %s (reactions=%d)%n", out, merged.getReactions().size());
            return;
        }

        List<MergedReaction> reactions2s = new ArrayList<>();
        List<MergedReaction> reactions3s = new ArrayList<>();
        Map<Integer, Integer> otherSetCounts = new LinkedHashMap<>();
        for (MergedReaction reaction : merge.reactions) {
            if (reaction.synthonSetCount == 2) {
                reactions2s.add(reaction);
            } else if (reaction.synthonSetCount == 3) {
                reactions3s.add(reaction);
            } else {
                otherSetCounts.merge(reaction.synthonSetCount, 1, Integer::sum);
            }
        }

        RawSynthonSpace merged2s = buildMergedSpace(baseName + "_2s", version, loaded, reactions2s, metadataPrefix, "split_2s");
        RawSynthonSpace merged3s = buildMergedSpace(baseName + "_3s", version, loaded, reactions3s, metadataPrefix, "split_3s");
        HyperspaceIOUtils.saveRawSynthonSpace(merged2s, out2s);
        HyperspaceIOUtils.saveRawSynthonSpace(merged3s, out3s);

        System.out.printf(Locale.ROOT, "Wrote 2-set merged raw space: %s (reactions=%d)%n", out2s, merged2s.getReactions().size());
        System.out.printf(Locale.ROOT, "Wrote 3-set merged raw space: %s (reactions=%d)%n", out3s, merged3s.getReactions().size());
        if (!otherSetCounts.isEmpty()) {
            System.out.printf(Locale.ROOT, "Skipped reactions with non-2/3 set counts: %s%n", otherSetCounts);
        }
    }

    private static MergeResult mergeSpaces(List<LoadedSpace> spaces, DuplicateStrategy duplicateStrategy) {
        Set<String> usedReactionIds = new LinkedHashSet<>();
        List<MergedReaction> merged = new ArrayList<>();
        int duplicateConflicts = 0;
        int renamedReactions = 0;

        for (LoadedSpace loaded : spaces) {
            List<String> reactionIds = new ArrayList<>(loaded.space.getReactions().keySet());
            Collections.sort(reactionIds);
            for (String originalReactionId : reactionIds) {
                RawSynthonSpace.ReactionData reactionData = loaded.space.getReactions().get(originalReactionId);
                if (reactionData == null) {
                    continue;
                }

                String mergedReactionId = originalReactionId;
                if (!usedReactionIds.add(mergedReactionId)) {
                    duplicateConflicts++;
                    if (duplicateStrategy == DuplicateStrategy.FAIL) {
                        throw new IllegalStateException("Duplicate reaction id detected: " + originalReactionId
                                + " (source=" + loaded.path + ")");
                    }
                    String suffix = sanitizeToken(loaded.sourceName);
                    String base = originalReactionId + "__" + suffix;
                    mergedReactionId = base;
                    int counter = 2;
                    while (!usedReactionIds.add(mergedReactionId)) {
                        mergedReactionId = base + "_" + counter;
                        counter++;
                    }
                    renamedReactions++;
                }

                int synthonSetCount = countNonEmptySets(reactionData.getRawFragmentSets());
                merged.add(new MergedReaction(
                        mergedReactionId,
                        originalReactionId,
                        synthonSetCount,
                        loaded.path,
                        loaded.sourceName,
                        reactionData
                ));
            }
        }

        return new MergeResult(merged, duplicateConflicts, renamedReactions);
    }

    private static RawSynthonSpace buildMergedSpace(String name,
                                                    String version,
                                                    List<LoadedSpace> sources,
                                                    List<MergedReaction> reactions,
                                                    String metadataPrefix,
                                                    String mergeMode) {
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(name).version(version);
        builder.putMetadata("merge.mode", mergeMode);
        builder.putMetadata("merge.sourceCount", Integer.toString(sources.size()));
        builder.putMetadata("merge.reactionCount", Integer.toString(reactions.size()));

        for (LoadedSpace source : sources) {
            String prefix = "merge.source." + source.index + ".";
            builder.putMetadata(prefix + "name", source.sourceName);
            builder.putMetadata(prefix + "path", source.path.toString());
        }

        for (MergedReaction reaction : reactions) {
            copyReaction(builder, reaction, metadataPrefix);
        }
        return builder.build();
    }

    private static void copyReaction(RawSynthonSpace.Builder builder,
                                     MergedReaction reaction,
                                     String metadataPrefix) {
        String reactionId = reaction.mergedReactionId;
        RawSynthonSpace.ReactionData source = reaction.sourceReactionData;

        source.getRawFragmentSets().forEach((setIndex, synthons) ->
                builder.addRawFragments(reactionId, setIndex, rewriteReactionIds(synthons, reactionId, setIndex)));
        source.getRawDownsampledSets().forEach((setIndex, synthons) ->
                builder.addRawDownsampledFragments(reactionId, setIndex, rewriteReactionIds(synthons, reactionId, setIndex)));
        source.getExampleScaffolds().forEach(scaffold ->
                builder.addExampleScaffolds(reactionId, List.of(scaffold)));
        source.getPartialAssemblies().forEach((missingIndex, assemblies) ->
                builder.addPartialAssemblies(reactionId, missingIndex, assemblies));
        if (!source.getRepresentativeCompounds().isEmpty()) {
            builder.addRepresentativeCompounds(reactionId, source.getRepresentativeCompounds());
        }
        source.getDescriptors().forEach((key, value) -> builder.addReactionDescriptor(reactionId, key, value));
        source.getReactionMetadata().forEach((key, value) -> builder.addReactionMetadata(reactionId, key, value));
        source.getFragmentAttributes().forEach((fragmentId, attributes) ->
                attributes.forEach((key, value) -> builder.addFragmentAttribute(reactionId, fragmentId, key, value)));

        builder.addReactionMetadata(reactionId, metadataPrefix + ".spaceName", reaction.sourceSpaceName);
        builder.addReactionMetadata(reactionId, metadataPrefix + ".spacePath", reaction.sourcePath.toString());
        builder.addReactionMetadata(reactionId, metadataPrefix + ".originalReactionId", reaction.originalReactionId);
        builder.addReactionMetadata(reactionId, metadataPrefix + ".synthonSetCount", Integer.toString(reaction.synthonSetCount));
        if (!reaction.mergedReactionId.equals(reaction.originalReactionId)) {
            builder.addReactionMetadata(reactionId, metadataPrefix + ".mergedReactionId", reaction.mergedReactionId);
        }
    }

    private static List<RawSynthon> rewriteReactionIds(List<RawSynthon> synthons,
                                                       String reactionId,
                                                       int setIndex) {
        List<RawSynthon> rewritten = new ArrayList<>(synthons.size());
        for (RawSynthon synthon : synthons) {
            rewritten.add(new RawSynthon(
                    reactionId,
                    setIndex,
                    synthon.getFragmentId(),
                    synthon.getIdcode(),
                    synthon.getConnectors()
            ));
        }
        return rewritten;
    }

    private static int countNonEmptySets(Map<Integer, List<RawSynthon>> sets) {
        int count = 0;
        for (List<RawSynthon> synthons : sets.values()) {
            if (synthons != null && !synthons.isEmpty()) {
                count++;
            }
        }
        return count;
    }

    private static String deriveSourceName(RawSynthonSpace raw, Path path) {
        if (raw.getName() != null && !raw.getName().isBlank()) {
            return raw.getName();
        }
        Path file = path.getFileName();
        return file == null ? "source" : file.toString();
    }

    private static List<Path> parseInputFiles(CommandLine cmd) throws Exception {
        LinkedHashSet<Path> deduplicated = new LinkedHashSet<>();
        String[] rawInValues = cmd.getOptionValues("rawIn");
        if (rawInValues != null) {
            for (String value : rawInValues) {
                if (value == null) {
                    continue;
                }
                for (String token : value.split(",")) {
                    String trimmed = token.trim();
                    if (!trimmed.isEmpty()) {
                        deduplicated.add(Path.of(trimmed).toAbsolutePath().normalize());
                    }
                }
            }
        }

        String inputDirectory = cmd.getOptionValue("inputDirectory");
        if (!isBlank(inputDirectory)) {
            Path directory = Path.of(inputDirectory).toAbsolutePath().normalize();
            if (!Files.isDirectory(directory)) {
                throw new IllegalArgumentException("Input directory does not exist or is not a directory: " + directory);
            }
            try (Stream<Path> stream = Files.list(directory)) {
                List<Path> discovered = stream
                        .filter(Files::isRegularFile)
                        .filter(path -> path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".rawspace.gz"))
                        .map(path -> path.toAbsolutePath().normalize())
                        .sorted(Comparator.comparing(Path::toString))
                        .collect(Collectors.toList());
                deduplicated.addAll(discovered);
            }
        }

        String manifest = cmd.getOptionValue("manifest");
        if (!isBlank(manifest)) {
            Path manifestPath = Path.of(manifest).toAbsolutePath().normalize();
            Path parent = manifestPath.getParent();
            for (String line : Files.readAllLines(manifestPath, StandardCharsets.UTF_8)) {
                String trimmed = line.trim();
                if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                    continue;
                }
                Path candidate = Path.of(trimmed);
                if (!candidate.isAbsolute() && parent != null) {
                    candidate = parent.resolve(candidate);
                }
                deduplicated.add(candidate.toAbsolutePath().normalize());
            }
        }

        return new ArrayList<>(deduplicated);
    }

    private static DuplicateStrategy parseDuplicateStrategy(String raw) {
        String value = raw == null ? "rename" : raw.trim().toLowerCase(Locale.ROOT);
        return switch (value) {
            case "rename" -> DuplicateStrategy.RENAME;
            case "fail" -> DuplicateStrategy.FAIL;
            default -> throw new IllegalArgumentException("Unsupported --duplicateStrategy: " + raw
                    + " (expected rename|fail)");
        };
    }

    private static String sanitizeToken(String token) {
        if (token == null || token.isBlank()) {
            return "source";
        }
        return token.replaceAll("[^a-zA-Z0-9._-]+", "_");
    }

    private static boolean isBlank(String value) {
        return value == null || value.isBlank();
    }

    private static CommandLine parse(Options options, String[] args) {
        try {
            return new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("SynthonSpaceMergeCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawIn").hasArg()
                .desc("Input RawSynthonSpace file(s), repeatable or comma-separated").build());
        options.addOption(Option.builder().longOpt("inputDirectory").hasArg()
                .desc("Input directory; all *.rawspace.gz files in this folder are included").build());
        options.addOption(Option.builder().longOpt("manifest").hasArg()
                .desc("Optional text file listing one input path per line").build());
        options.addOption(Option.builder().longOpt("splitBySetCount")
                .desc("Write split outputs: reactions with 2 sets and 3 sets").build());
        options.addOption(Option.builder().longOpt("out").hasArg()
                .desc("Output file for single merged space (required without --splitBySetCount)").build());
        options.addOption(Option.builder().longOpt("out2s").hasArg()
                .desc("Output file for merged 2-set reactions (required with --splitBySetCount)").build());
        options.addOption(Option.builder().longOpt("out3s").hasArg()
                .desc("Output file for merged 3-set reactions (required with --splitBySetCount)").build());
        options.addOption(Option.builder().longOpt("spaceName").hasArg()
                .desc("Base name for merged space(s) (default merged_rawspace)").build());
        options.addOption(Option.builder().longOpt("version").hasArg()
                .desc("Version for merged space(s) (default 1.0)").build());
        options.addOption(Option.builder().longOpt("duplicateStrategy").hasArg()
                .desc("How to handle duplicate reaction ids: rename (default) or fail").build());
        options.addOption(Option.builder().longOpt("sourceMetadataPrefix").hasArg()
                .desc("Prefix for provenance metadata keys at reaction level (default source)").build());
        options.addOption(Option.builder().longOpt("synthonListOut").hasArg()
                .desc("Optional TSV export of unique synthons across loaded spaces").build());
        return options;
    }

    private static Map<String, SynthonStats> collectSynthonStats(List<LoadedSpace> loadedSpaces) {
        Map<String, SynthonStatsAccumulator> accumulators = new TreeMap<>();
        for (LoadedSpace loaded : loadedSpaces) {
            loaded.space.getReactions().forEach((reactionId, reactionData) -> {
                if (reactionData == null) {
                    return;
                }
                reactionData.getRawFragmentSets().forEach((setIndex, synthons) -> {
                    if (synthons == null || synthons.isEmpty()) {
                        return;
                    }
                    String synthonSetKey = loaded.index + "|" + reactionId + "|" + setIndex;
                    for (RawSynthon synthon : synthons) {
                        if (synthon == null || isBlank(synthon.getIdcode())) {
                            continue;
                        }
                        accumulators.computeIfAbsent(synthon.getIdcode(),
                                        key -> new SynthonStatsAccumulator(computeNumAtoms(synthon.getIdcode())))
                                .addOccurrence(loaded.index, synthonSetKey);
                    }
                });
            });
        }

        Map<String, SynthonStats> result = new TreeMap<>();
        accumulators.forEach((idcode, accumulator) ->
                result.put(idcode, new SynthonStats(accumulator.numAtoms,
                        accumulator.spaceIndices.size(),
                        accumulator.synthonSetKeys.size())));
        return result;
    }

    private static void writeSynthonList(Path output, Map<String, SynthonStats> stats) throws Exception {
        Path parent = output.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("Structure[idcode]\tNumAtoms\tNumSpaces\tNumSynthonSets");
            writer.newLine();
            for (Map.Entry<String, SynthonStats> entry : stats.entrySet()) {
                SynthonStats value = entry.getValue();
                writer.write(entry.getKey());
                writer.write('\t');
                writer.write(Integer.toString(value.numAtoms()));
                writer.write('\t');
                writer.write(Integer.toString(value.numSpaces()));
                writer.write('\t');
                writer.write(Integer.toString(value.numSynthonSets()));
                writer.newLine();
            }
        }
    }

    private static int computeNumAtoms(String idcode) {
        try {
            StereoMolecule mol = new StereoMolecule();
            new IDCodeParser().parse(mol, idcode);
            return mol.getAtoms();
        } catch (Exception ignored) {
            return -1;
        }
    }

    private enum DuplicateStrategy {
        RENAME,
        FAIL
    }

    private record LoadedSpace(int index,
                               Path path,
                               RawSynthonSpace space,
                               String sourceName) {
    }

    private record MergedReaction(String mergedReactionId,
                                  String originalReactionId,
                                  int synthonSetCount,
                                  Path sourcePath,
                                  String sourceSpaceName,
                                  RawSynthonSpace.ReactionData sourceReactionData) {
    }

    private record MergeResult(List<MergedReaction> reactions,
                               int duplicateConflicts,
                               int renamedReactions) {
    }

    private record SynthonStats(int numAtoms,
                                int numSpaces,
                                int numSynthonSets) {
    }

    private static final class SynthonStatsAccumulator {
        private final int numAtoms;
        private final Set<Integer> spaceIndices = new LinkedHashSet<>();
        private final Set<String> synthonSetKeys = new LinkedHashSet<>();

        private SynthonStatsAccumulator(int numAtoms) {
            this.numAtoms = numAtoms;
        }

        private void addOccurrence(int spaceIndex, String synthonSetKey) {
            spaceIndices.add(spaceIndex);
            synthonSetKeys.add(synthonSetKey);
        }
    }
}
