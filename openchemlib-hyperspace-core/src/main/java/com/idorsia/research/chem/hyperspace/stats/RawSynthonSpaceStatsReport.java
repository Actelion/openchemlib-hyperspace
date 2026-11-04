package com.idorsia.research.chem.hyperspace.stats;

import java.io.BufferedWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public final class RawSynthonSpaceStatsReport {
    private final SpaceSummary spaceSummary;
    private final List<ReactionStats> reactionStats;
    private final List<SourceStats> sourceStats;
    private final List<SynthonSetStats> synthonSetStats;
    private final List<ExampleProduct> exampleProducts;

    RawSynthonSpaceStatsReport(SpaceSummary spaceSummary,
                               List<ReactionStats> reactionStats,
                               List<SourceStats> sourceStats,
                               List<SynthonSetStats> synthonSetStats,
                               List<ExampleProduct> exampleProducts) {
        this.spaceSummary = spaceSummary;
        this.reactionStats = List.copyOf(reactionStats);
        this.sourceStats = List.copyOf(sourceStats);
        this.synthonSetStats = List.copyOf(synthonSetStats);
        this.exampleProducts = List.copyOf(exampleProducts);
    }

    public SpaceSummary getSpaceSummary() {
        return spaceSummary;
    }

    public List<ReactionStats> getReactionStats() {
        return reactionStats;
    }

    public List<SourceStats> getSourceStats() {
        return sourceStats;
    }

    public List<SynthonSetStats> getSynthonSetStats() {
        return synthonSetStats;
    }

    public List<ExampleProduct> getExampleProducts() {
        return exampleProducts;
    }

    public void writeToDirectory(Path outputDirectory) throws IOException {
        Files.createDirectories(outputDirectory);
        writeSummaryMarkdown(outputDirectory.resolve("summary.md"));
        writeSpaceSummary(outputDirectory.resolve("space_summary.tsv"));
        writeReactionStats(outputDirectory.resolve("reaction_stats.tsv"));
        writeSourceStats(outputDirectory.resolve("source_stats.tsv"));
        writeSynthonSetStats(outputDirectory.resolve("synthon_set_stats.tsv"));
        writeExampleProducts(outputDirectory.resolve("example_products.tsv"));
    }

    private void writeSummaryMarkdown(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("# Raw Synthon Space Statistics");
            writer.newLine();
            writer.newLine();
            writer.write("Generated: " + Instant.now());
            writer.newLine();
            writer.newLine();
            writer.write("## Space");
            writer.newLine();
            writer.newLine();
            writer.write("- Name: " + display(spaceSummary.name()));
            writer.newLine();
            writer.write("- Version: " + display(spaceSummary.version()));
            writer.newLine();
            writer.write("- Reactions: " + spaceSummary.reactionCount());
            writer.newLine();
            writer.write("- 2-set reactions: " + spaceSummary.twoSetReactionCount());
            writer.newLine();
            writer.write("- 3-set reactions: " + spaceSummary.threeSetReactionCount());
            writer.newLine();
            writer.write("- Other reactions: " + spaceSummary.otherReactionCount());
            writer.newLine();
            writer.write("- Total synthons: " + spaceSummary.totalSynthons());
            writer.newLine();
            writer.write("- Unique synthon IDCodes: " + spaceSummary.uniqueSynthonIdcodes());
            writer.newLine();
            writer.write("- Total combinatorial products: " + spaceSummary.totalProductCount());
            writer.newLine();
            writer.write("- Product sampling: " + spaceSummary.productSamplesPerReaction() + " attempts per reaction");
            writer.newLine();
            writer.write("- Example products: up to " + spaceSummary.examplesPerReaction() + " per reaction");
            writer.newLine();
            writer.newLine();
            writer.write("## Files");
            writer.newLine();
            writer.newLine();
            writer.write("- `space_summary.tsv`: one-row space totals and rawspace metadata.");
            writer.newLine();
            writer.write("- `reaction_stats.tsv`: exact reaction sizes plus sampled product-property summaries.");
            writer.newLine();
            writer.write("- `source_stats.tsv`: reaction and synthon totals grouped by source metadata.");
            writer.newLine();
            writer.write("- `synthon_set_stats.tsv`: per-reaction/per-set synthon counts and atom statistics.");
            writer.newLine();
            writer.write("- `example_products.tsv`: deterministic random assembled products.");
            writer.newLine();
        }
    }

    private void writeSpaceSummary(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("spaceName\tversion\treactions\treactions2s\treactions3s\treactionsOther\tsynthonSets\ttotalSynthons\tuniqueSynthonIdcodes\ttotalProductCount\texamplesPerReaction\tproductSamplesPerReaction\tseed\tmetadata");
            writer.newLine();
            writer.write(tsv(spaceSummary.name()));
            writer.write('\t');
            writer.write(tsv(spaceSummary.version()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.reactionCount()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.twoSetReactionCount()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.threeSetReactionCount()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.otherReactionCount()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.synthonSetCount()));
            writer.write('\t');
            writer.write(Long.toString(spaceSummary.totalSynthons()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.uniqueSynthonIdcodes()));
            writer.write('\t');
            writer.write(spaceSummary.totalProductCount().toString());
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.examplesPerReaction()));
            writer.write('\t');
            writer.write(Integer.toString(spaceSummary.productSamplesPerReaction()));
            writer.write('\t');
            writer.write(Long.toString(spaceSummary.seed()));
            writer.write('\t');
            writer.write(tsv(joinMetadata(spaceSummary.metadata())));
            writer.newLine();
        }
    }

    private void writeReactionStats(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("reactionId\tsourceSpace\tsourcePath\toriginalReactionId\tsynthonSetCount\tsetSizes\tproductCount\ttotalSynthons\tuniqueSynthonIdcodes\tsynthonAtomsMin\tsynthonAtomsMean\tsynthonAtomsMax\tproductSampleCount\tproductSampleFailures\tproductAtomsMin\tproductAtomsMean\tproductAtomsMax\tproductRotatableMin\tproductRotatableMean\tproductRotatableMax");
            writer.newLine();
            for (ReactionStats stats : reactionStats) {
                writer.write(tsv(stats.reactionId()));
                writer.write('\t');
                writer.write(tsv(stats.sourceSpace()));
                writer.write('\t');
                writer.write(tsv(stats.sourcePath()));
                writer.write('\t');
                writer.write(tsv(stats.originalReactionId()));
                writer.write('\t');
                writer.write(Integer.toString(stats.synthonSetCount()));
                writer.write('\t');
                writer.write(tsv(stats.setSizes()));
                writer.write('\t');
                writer.write(stats.productCount().toString());
                writer.write('\t');
                writer.write(Long.toString(stats.totalSynthons()));
                writer.write('\t');
                writer.write(Integer.toString(stats.uniqueSynthonIdcodes()));
                writer.write('\t');
                writer.write(formatInt(stats.synthonAtomsMin()));
                writer.write('\t');
                writer.write(formatDouble(stats.synthonAtomsMean()));
                writer.write('\t');
                writer.write(formatInt(stats.synthonAtomsMax()));
                writer.write('\t');
                writer.write(Integer.toString(stats.productSampleCount()));
                writer.write('\t');
                writer.write(Integer.toString(stats.productSampleFailures()));
                writer.write('\t');
                writer.write(formatInt(stats.productAtomsMin()));
                writer.write('\t');
                writer.write(formatDouble(stats.productAtomsMean()));
                writer.write('\t');
                writer.write(formatInt(stats.productAtomsMax()));
                writer.write('\t');
                writer.write(formatInt(stats.productRotatableMin()));
                writer.write('\t');
                writer.write(formatDouble(stats.productRotatableMean()));
                writer.write('\t');
                writer.write(formatInt(stats.productRotatableMax()));
                writer.newLine();
            }
        }
    }

    private void writeSourceStats(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("sourceSpace\treactions\tsynthonSets\ttotalSynthons\tuniqueSynthonIdcodes\ttotalProductCount\treactions2s\treactions3s\treactionsOther");
            writer.newLine();
            for (SourceStats stats : sourceStats) {
                writer.write(tsv(stats.sourceSpace()));
                writer.write('\t');
                writer.write(Integer.toString(stats.reactionCount()));
                writer.write('\t');
                writer.write(Integer.toString(stats.synthonSetCount()));
                writer.write('\t');
                writer.write(Long.toString(stats.totalSynthons()));
                writer.write('\t');
                writer.write(Integer.toString(stats.uniqueSynthonIdcodes()));
                writer.write('\t');
                writer.write(stats.totalProductCount().toString());
                writer.write('\t');
                writer.write(Integer.toString(stats.twoSetReactionCount()));
                writer.write('\t');
                writer.write(Integer.toString(stats.threeSetReactionCount()));
                writer.write('\t');
                writer.write(Integer.toString(stats.otherReactionCount()));
                writer.newLine();
            }
        }
    }

    private void writeSynthonSetStats(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("reactionId\tsourceSpace\tsetIndex\tsynthonCount\tuniqueSynthonIdcodes\tsynthonAtomsMin\tsynthonAtomsMean\tsynthonAtomsMax");
            writer.newLine();
            for (SynthonSetStats stats : synthonSetStats) {
                writer.write(tsv(stats.reactionId()));
                writer.write('\t');
                writer.write(tsv(stats.sourceSpace()));
                writer.write('\t');
                writer.write(Integer.toString(stats.setIndex()));
                writer.write('\t');
                writer.write(Integer.toString(stats.synthonCount()));
                writer.write('\t');
                writer.write(Integer.toString(stats.uniqueSynthonIdcodes()));
                writer.write('\t');
                writer.write(formatInt(stats.synthonAtomsMin()));
                writer.write('\t');
                writer.write(formatDouble(stats.synthonAtomsMean()));
                writer.write('\t');
                writer.write(formatInt(stats.synthonAtomsMax()));
                writer.newLine();
            }
        }
    }

    private void writeExampleProducts(Path output) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
            writer.write("reactionId\tsourceSpace\texampleIndex\tfragIds\tStructure [idcode]\tatoms\trotatableBonds");
            writer.newLine();
            for (ExampleProduct example : exampleProducts) {
                writer.write(tsv(example.reactionId()));
                writer.write('\t');
                writer.write(tsv(example.sourceSpace()));
                writer.write('\t');
                writer.write(Integer.toString(example.exampleIndex()));
                writer.write('\t');
                writer.write(tsv(example.fragmentIds()));
                writer.write('\t');
                writer.write(tsv(example.idcode()));
                writer.write('\t');
                writer.write(Integer.toString(example.atoms()));
                writer.write('\t');
                writer.write(Integer.toString(example.rotatableBonds()));
                writer.newLine();
            }
        }
    }

    private static String joinMetadata(Map<String, String> metadata) {
        StringBuilder builder = new StringBuilder();
        metadata.forEach((key, value) -> {
            if (builder.length() > 0) {
                builder.append(';');
            }
            builder.append(key).append('=').append(value);
        });
        return builder.toString();
    }

    private static String tsv(String value) {
        if (value == null) {
            return "";
        }
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    private static String display(String value) {
        return value == null || value.isBlank() ? "(unknown)" : value;
    }

    private static String formatInt(int value) {
        return value < 0 ? "" : Integer.toString(value);
    }

    private static String formatDouble(double value) {
        if (Double.isNaN(value)) {
            return "";
        }
        return String.format(Locale.ROOT, "%.3f", value);
    }

    public record SpaceSummary(String name,
                               String version,
                               Map<String, String> metadata,
                               int reactionCount,
                               int twoSetReactionCount,
                               int threeSetReactionCount,
                               int otherReactionCount,
                               int synthonSetCount,
                               long totalSynthons,
                               int uniqueSynthonIdcodes,
                               BigInteger totalProductCount,
                               int examplesPerReaction,
                               int productSamplesPerReaction,
                               long seed) {
    }

    public record ReactionStats(String reactionId,
                                String sourceSpace,
                                String sourcePath,
                                String originalReactionId,
                                int synthonSetCount,
                                String setSizes,
                                BigInteger productCount,
                                long totalSynthons,
                                int uniqueSynthonIdcodes,
                                int synthonAtomsMin,
                                double synthonAtomsMean,
                                int synthonAtomsMax,
                                int productSampleCount,
                                int productSampleFailures,
                                int productAtomsMin,
                                double productAtomsMean,
                                int productAtomsMax,
                                int productRotatableMin,
                                double productRotatableMean,
                                int productRotatableMax) {
    }

    public record SourceStats(String sourceSpace,
                              int reactionCount,
                              int synthonSetCount,
                              long totalSynthons,
                              int uniqueSynthonIdcodes,
                              BigInteger totalProductCount,
                              int twoSetReactionCount,
                              int threeSetReactionCount,
                              int otherReactionCount) {
    }

    public record SynthonSetStats(String reactionId,
                                  String sourceSpace,
                                  int setIndex,
                                  int synthonCount,
                                  int uniqueSynthonIdcodes,
                                  int synthonAtomsMin,
                                  double synthonAtomsMean,
                                  int synthonAtomsMax) {
    }

    public record ExampleProduct(String reactionId,
                                 String sourceSpace,
                                 int exampleIndex,
                                 String fragmentIds,
                                 String idcode,
                                 int atoms,
                                 int rotatableBonds) {
    }
}
