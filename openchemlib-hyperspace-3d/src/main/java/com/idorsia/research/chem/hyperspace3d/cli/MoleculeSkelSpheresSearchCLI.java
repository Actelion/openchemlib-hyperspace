package com.idorsia.research.chem.hyperspace3d.cli;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeExactSkelSpheresReranker;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresHit;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresScreener;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Exhaustive learned SkelSpheres screen of a flat molecule fingerprint index. */
public final class MoleculeSkelSpheresSearchCLI {
    private MoleculeSkelSpheresSearchCLI() {}

    public static void main(String[] args) throws Exception {
        long started = System.nanoTime();
        Path configPath = config(args).toAbsolutePath().normalize();
        var mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);
        var config = mapper.readValue(configPath.toFile(),
                MoleculeSkelSpheresSearchConfig.class);
        config.validate();
        var paths = config.resolve(configPath);
        var index = MoleculeFingerprintIndexReader.loadManifest(paths.index());
        requirePath(index.modelBundle, paths.modelBundle(), "primary model bundle");
        requirePath(index.compactBundle, paths.compactBundle(), "compact model bundle");
        var source = DeepSpaceModelBundle.load(paths.modelBundle());
        var compact = CompactSkelSpheresModelBundle.load(paths.compactBundle(),
                source.manifest());
        StereoMolecule query = parse(config.query.structure, config.query.format);

        long queryStarted = System.nanoTime();
        float[] queryVector;
        var runtime = new DeepSpaceOnnxEnvironment(config.device(),
                config.runtime.cudaDeviceId);
        try (var encoder = new DeepSpaceV1Encoder(runtime, source);
             var projector = new CompactSkelSpheresProjector(runtime, compact)) {
            var batch = new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer())
                    .build(List.of(query));
            queryVector = projector.project(encoder.encode(batch)[0]);
        }
        long queryNanos = System.nanoTime() - queryStarted;

        var screening = new MoleculeSkelSpheresScreener(
                new CompactSkelSpheresScorer(compact.manifest())).screen(
                        paths.index(), index, queryVector, config.output.learnedTopK,
                        config.runtime.scanBatchSize);
        Map<MoleculeVectorReference, Integer> learnedRanks = new HashMap<>();
        for (int i = 0; i < screening.hits().size(); i++) {
            learnedRanks.put(screening.hits().get(i).reference(), i + 1);
        }
        long exactStarted = System.nanoTime();
        List<MoleculeSkelSpheresHit> exactHits =
                new MoleculeExactSkelSpheresReranker().rerank(query, screening.hits());
        long exactNanos = System.nanoTime() - exactStarted;
        double spearman = spearman(screening.hits(), exactHits);

        createParents(paths);
        writeTsv(paths.hitsTsv(), exactHits, learnedRanks);
        writeSdf(paths.hitsSdf(), exactHits, learnedRanks);
        long totalNanos = System.nanoTime() - started;
        writeSummary(paths.summaryMarkdown(), config, index.recordCount, queryNanos,
                screening.scanNanos(), screening.metadataResolutionNanos(), exactNanos,
                totalNanos, spearman, exactHits, learnedRanks);

        Map<String, Object> run = new LinkedHashMap<>();
        run.put("artifactType", "hyperspace-molecule-skelspheres16-search-run");
        run.put("formatVersion", 1);
        run.put("queryIdentifier", config.query.identifier);
        run.put("queryStructure", config.query.structure);
        run.put("queryFormat", config.query.format.toLowerCase(Locale.ROOT));
        run.put("queryIdcode", new Canonizer(query).getIDCode());
        run.put("index", paths.index().toString());
        run.put("indexManifestSha256", ProductFingerprintIndexBuilder.sha256(
                paths.index().resolve("manifest.json")));
        run.put("modelBundle", paths.modelBundle().toString());
        run.put("compactBundle", paths.compactBundle().toString());
        run.put("recordsScanned", screening.recordsScanned());
        run.put("learnedTopK", config.output.learnedTopK);
        run.put("exactScoresComputed", exactHits.size());
        run.put("shortlistSpearman", spearman);
        run.put("queryEncodingNanos", queryNanos);
        run.put("vectorScanNanos", screening.scanNanos());
        run.put("metadataResolutionNanos", screening.metadataResolutionNanos());
        run.put("exactRerankNanos", exactNanos);
        run.put("totalNanos", totalNanos);
        run.put("comparisonsPerSecond", screening.recordsScanned()
                / (screening.scanNanos() / 1.0e9));
        run.put("outputs", Map.of("hitsTsv", paths.hitsTsv().toString(),
                "hitsSdf", paths.hitsSdf().toString(),
                "summaryMarkdown", paths.summaryMarkdown().toString()));
        mapper.writeValue(paths.runManifest().toFile(), run);
        System.out.printf(Locale.ROOT,
                "Scanned %,d molecules at %,.0f comparisons/s; shortlist Spearman %.4f%n",
                screening.recordsScanned(), screening.recordsScanned()
                        / (screening.scanNanos() / 1.0e9), spearman);
        System.out.printf("Results: %s%n", paths.hitsTsv());
    }

    private static StereoMolecule parse(String structure, String format) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        if ("smiles".equalsIgnoreCase(format)) new SmilesParser().parse(molecule, structure);
        else new IDCodeParser().parse(molecule, structure);
        return molecule;
    }

    private static void requirePath(String indexed, Path selected, String label) {
        if (!Path.of(indexed).toAbsolutePath().normalize().equals(selected)) {
            throw new IllegalArgumentException(label + " does not match index manifest");
        }
    }

    private static void createParents(MoleculeSkelSpheresSearchConfig.ResolvedPaths paths)
            throws Exception {
        for (Path path : List.of(paths.hitsTsv(), paths.hitsSdf(),
                paths.summaryMarkdown(), paths.runManifest())) {
            if (path.getParent() != null) Files.createDirectories(path.getParent());
        }
    }

    private static void writeTsv(Path path, List<MoleculeSkelSpheresHit> hits,
            Map<MoleculeVectorReference, Integer> learnedRanks) throws Exception {
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("exact_rank\tlearned_rank\tmolecule_id\tsmiles\tsource_row\theavy_atoms\t"
                    + "dot_product\tpredicted_skelspheres\texact_skelspheres\n");
            for (int i = 0; i < hits.size(); i++) {
                var hit = hits.get(i);
                out.write((i + 1) + "\t" + learnedRanks.get(hit.reference()) + "\t"
                        + hit.molecule().moleculeId() + "\t" + hit.molecule().smiles() + "\t"
                        + hit.molecule().sourceRow() + "\t" + hit.molecule().heavyAtomCount()
                        + "\t" + format(hit.dotProduct()) + "\t"
                        + format(hit.predictedSimilarity()) + "\t"
                        + format(hit.exactSimilarity()) + "\n");
            }
        }
    }

    private static void writeSdf(Path path, List<MoleculeSkelSpheresHit> hits,
            Map<MoleculeVectorReference, Integer> learnedRanks) throws Exception {
        var parser = new SmilesParser();
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            for (int i = 0; i < hits.size(); i++) {
                var hit = hits.get(i);
                StereoMolecule molecule = new StereoMolecule();
                parser.parse(molecule, hit.molecule().smiles());
                molecule.setName(hit.molecule().moleculeId());
                new CoordinateInventor().invent(molecule);
                new MolfileCreator(molecule).writeMolfile(out);
                property(out, "ENAMINE_ID", hit.molecule().moleculeId());
                property(out, "SMILES", hit.molecule().smiles());
                property(out, "EXACT_RANK", Integer.toString(i + 1));
                property(out, "LEARNED_RANK",
                        Integer.toString(learnedRanks.get(hit.reference())));
                property(out, "PREDICTED_SKELSPHERES", format(hit.predictedSimilarity()));
                property(out, "EXACT_SKELSPHERES", format(hit.exactSimilarity()));
                out.write("$$$$\n");
            }
        }
    }

    private static void property(BufferedWriter out, String name, String value)
            throws Exception {
        out.write("> <" + name + ">\n" + value + "\n\n");
    }

    private static void writeSummary(Path path, MoleculeSkelSpheresSearchConfig config,
            long records, long queryNanos, long scanNanos, long resolutionNanos,
            long exactNanos, long totalNanos, double spearman,
            List<MoleculeSkelSpheresHit> hits,
            Map<MoleculeVectorReference, Integer> learnedRanks) throws Exception {
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("# Enamine 13.2M learned SkelSpheres screen\n\n");
            out.write("- Query: `" + config.query.identifier + "`\n");
            out.write("- Structure: `" + config.query.structure + "`\n");
            out.write(String.format(Locale.ROOT,
                    "- Scanned: %,d molecules in %.3f s (%,.0f comparisons/s)\n",
                    records, scanNanos / 1.0e9, records / (scanNanos / 1.0e9)));
            out.write(String.format(Locale.ROOT,
                    "- Query encoding: %.3f s; metadata: %.3f s; exact top-%d: %.3f s; total: %.3f s\n",
                    queryNanos / 1.0e9, resolutionNanos / 1.0e9, hits.size(),
                    exactNanos / 1.0e9, totalNanos / 1.0e9));
            out.write(String.format(Locale.ROOT,
                    "- Predicted/exact Spearman within learned shortlist: %.4f\n\n", spearman));
            out.write("Exact scores cover only the learned shortlist and are not a global exact-search recall estimate.\n\n");
            out.write("| Exact rank | Learned rank | Enamine ID | Predicted | Exact | SMILES |\n");
            out.write("|---:|---:|---|---:|---:|---|\n");
            int shown = Math.min(config.output.reportTopK, hits.size());
            for (int i = 0; i < shown; i++) {
                var hit = hits.get(i);
                out.write("| " + (i + 1) + " | " + learnedRanks.get(hit.reference())
                        + " | " + hit.molecule().moleculeId() + " | "
                        + format(hit.predictedSimilarity()) + " | "
                        + format(hit.exactSimilarity()) + " | `"
                        + hit.molecule().smiles().replace("|", "\\|") + "` |\n");
            }
        }
    }

    private static double spearman(List<MoleculeSkelSpheresHit> learned,
            List<MoleculeSkelSpheresHit> exact) {
        Map<MoleculeVectorReference, Double> exactScores = new HashMap<>();
        for (var hit : exact) exactScores.put(hit.reference(), hit.exactSimilarity());
        double[] predicted = new double[learned.size()];
        double[] observed = new double[learned.size()];
        for (int i = 0; i < learned.size(); i++) {
            predicted[i] = learned.get(i).predictedSimilarity();
            observed[i] = exactScores.get(learned.get(i).reference());
        }
        return pearson(ranks(predicted), ranks(observed));
    }

    private static double[] ranks(double[] values) {
        Integer[] order = new Integer[values.length];
        for (int i = 0; i < values.length; i++) order[i] = i;
        Arrays.sort(order, Comparator.comparingDouble(i -> values[i]));
        double[] ranks = new double[values.length];
        int start = 0;
        while (start < order.length) {
            int end = start + 1;
            while (end < order.length
                    && Double.compare(values[order[start]], values[order[end]]) == 0) end++;
            double rank = (start + end - 1) / 2.0 + 1.0;
            for (int i = start; i < end; i++) ranks[order[i]] = rank;
            start = end;
        }
        return ranks;
    }

    private static double pearson(double[] left, double[] right) {
        double leftMean = Arrays.stream(left).average().orElse(Double.NaN);
        double rightMean = Arrays.stream(right).average().orElse(Double.NaN);
        double numerator = 0, leftSquared = 0, rightSquared = 0;
        for (int i = 0; i < left.length; i++) {
            double a = left[i] - leftMean;
            double b = right[i] - rightMean;
            numerator += a * b;
            leftSquared += a * a;
            rightSquared += b * b;
        }
        return numerator / Math.sqrt(leftSquared * rightSquared);
    }

    private static String format(double value) {
        return String.format(Locale.ROOT, "%.8f", value);
    }

    private static Path config(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        throw new IllegalArgumentException(
                "Usage: MoleculeSkelSpheresSearchCLI --config <search-molecules.json>");
    }
}
