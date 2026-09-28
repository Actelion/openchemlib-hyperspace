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
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSources;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Comparator;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeExactSkelSpheresReranker;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculePheSAHit;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculePheSAScreener;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresHit;
import com.idorsia.research.chem.hyperspace3d.screening.MoleculeSkelSpheresScreener;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Unified flat-library learned 2D, learned 3D, and cascade search command. */
public final class MoleculeSimilaritySearchCLI {
    private MoleculeSimilaritySearchCLI() {}

    public static void main(String[] args) throws Exception {
        long started = System.nanoTime();
        Path configPath = config(args).toAbsolutePath().normalize();
        var mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);
        var config = mapper.readValue(configPath.toFile(), MoleculeSimilaritySearchConfig.class);
        config.validate();
        var paths = config.resolve(configPath);
        var model = DeepSpaceModelBundle.load(paths.modelBundle());
        CompactSkelSpheresModelBundle compact = config.usesCompact()
                ? CompactSkelSpheresModelBundle.load(paths.compactBundle(), model.manifest()) : null;
        StereoMolecule query = parse(config.query.structure, config.query.format);
        var runtime = new DeepSpaceOnnxEnvironment(config.device(), config.runtime.cudaDeviceId);

        long queryStarted = System.nanoTime();
        float[] query128;
        try (var encoder = new DeepSpaceV1Encoder(runtime, model)) {
            var batch = new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer())
                    .build(List.of(query));
            query128 = encoder.encode(batch)[0];
        }
        float[] query16 = null;
        if (compact != null) {
            try (var projector = new CompactSkelSpheresProjector(runtime, compact)) {
                query16 = projector.project(query128);
            }
        }
        long queryNanos = System.nanoTime() - queryStarted;

        List<RenderedHit> hits;
        Map<String, Object> timings = new LinkedHashMap<>();
        long recordsScanned;
        try (var source = MoleculeFingerprintDataSources.open(paths.index())) {
            source.validateCompatibility(model, compact, config.usesCompact());
            switch (config.mode()) {
                case SKELSPHERES_2D -> {
                    var result = new MoleculeSkelSpheresScreener(
                            new CompactSkelSpheresScorer(compact.manifest())).screen(
                                    source, query16, config.screening.resultTopK,
                                    config.runtime.scanBatchSize, config.runtime.prefetchDepth);
                    var exact = new MoleculeExactSkelSpheresReranker().rerank(query, result.hits());
                    hits = render2d(exact);
                    recordsScanned = result.recordsScanned();
                    timings.put("compactScanNanos", result.scanNanos());
                    timings.put("metadataResolutionNanos", result.metadataResolutionNanos());
                }
                case PHESA_3D -> {
                    try (var comparator = new DeepSpaceV1Comparator(runtime, model)) {
                        var result = new MoleculePheSAScreener(comparator, model.manifest(),
                                config.objective()).screen(source, query128,
                                        config.screening.resultTopK,
                                        config.runtime.comparatorBatchSize, config.runtime.prefetchDepth);
                        hits = render3d(result.hits(), model, Map.of());
                        recordsScanned = result.recordsScanned();
                        add3dTimings(timings, result);
                    }
                }
                case SKELSPHERES_THEN_PHESA -> {
                    var compactResult = new MoleculeSkelSpheresScreener(
                            new CompactSkelSpheresScorer(compact.manifest())).scanUnresolved(
                                    source, query16, config.screening.cascadeShortlistTopK,
                                    config.runtime.scanBatchSize, config.runtime.prefetchDepth);
                    try (var comparator = new DeepSpaceV1Comparator(runtime, model)) {
                        var result = new MoleculePheSAScreener(comparator, model.manifest(),
                                config.objective()).rerank(source, query128, compactResult.hits(),
                                        config.screening.resultTopK,
                                        config.runtime.comparatorBatchSize);
                        Map<com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference,
                                Double> exactScores = exactSkelSpheres(query, result.hits());
                        hits = render3d(result.hits(), model, exactScores);
                        recordsScanned = compactResult.recordsScanned();
                        timings.put("compactScanNanos", compactResult.scanNanos());
                        add3dTimings(timings, result);
                    }
                }
                default -> throw new IllegalStateException("unhandled search mode");
            }
            timings.put("queryEncodingNanos", queryNanos);
            timings.put("totalNanos", System.nanoTime() - started);
            createParents(paths);
            writeTsv(paths.hitsTsv(), hits);
            writeSdf(paths.hitsSdf(), hits);
            writeSummary(paths.summaryMarkdown(), config, source.artifactType(),
                    source.recordCount(), recordsScanned, hits, timings);
            writeRunManifest(mapper, paths, config, source.artifactType(), source.recordCount(),
                    recordsScanned, model, compact, query, hits.size(), timings);
        }
        System.out.printf(Locale.ROOT, "%s screened %,d molecules; wrote %d hits to %s%n",
                config.mode(), recordsScanned, hits.size(), paths.hitsTsv());
    }

    private static void add3dTimings(Map<String, Object> timings,
            com.idorsia.research.chem.hyperspace3d.screening.MoleculePheSAScreeningResult result) {
        timings.put("baseVectorPreparationNanos", result.vectorPreparationNanos());
        timings.put("comparatorNanos", result.comparatorNanos());
        timings.put("metadataResolutionNanos", result.metadataResolutionNanos());
    }

    private static Map<com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference, Double>
            exactSkelSpheres(StereoMolecule query, List<MoleculePheSAHit> hits) throws Exception {
        List<MoleculeSkelSpheresHit> compact = hits.stream().map(hit ->
                new MoleculeSkelSpheresHit(hit.reference(), hit.molecule(),
                        hit.skelSpheresDotProduct(), hit.predictedSkelSpheres(), null)).toList();
        var exact = new MoleculeExactSkelSpheresReranker().rerank(query, compact);
        Map<com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference, Double> result =
                new HashMap<>();
        for (var hit : exact) result.put(hit.reference(), hit.exactSimilarity());
        return result;
    }

    private static List<RenderedHit> render2d(List<MoleculeSkelSpheresHit> hits) {
        List<RenderedHit> result = new ArrayList<>(hits.size());
        for (int index = 0; index < hits.size(); index++) {
            var hit = hits.get(index);
            result.add(new RenderedHit(index + 1, null, hit.reference(), hit.molecule(),
                    hit.predictedSimilarity(), hit.predictedSimilarity(), hit.exactSimilarity(),
                    null, null, null));
        }
        return result;
    }

    private static List<RenderedHit> render3d(List<MoleculePheSAHit> hits,
            DeepSpaceModelBundle model,
            Map<com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference, Double> exact) {
        int total = model.manifest().targetIndex("phesa_total");
        int shape = model.manifest().targetIndex("phesa_shape");
        int pharmacophore = model.manifest().targetIndex("phesa_pharmacophore");
        List<RenderedHit> result = new ArrayList<>(hits.size());
        for (int index = 0; index < hits.size(); index++) {
            var hit = hits.get(index);
            result.add(new RenderedHit(index + 1, hit.compactRank(), hit.reference(), hit.molecule(),
                    hit.objectiveScore(), hit.predictedSkelSpheres(), exact.get(hit.reference()),
                    (double) hit.modelScores()[total], (double) hit.modelScores()[shape],
                    (double) hit.modelScores()[pharmacophore]));
        }
        return result;
    }

    private static void writeTsv(Path path, List<RenderedHit> hits) throws Exception {
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("rank\tcompact_rank\tmolecule_id\tsmiles\tcanonical_smiles\tsource_row\t"
                    + "heavy_atoms\tshard\tlocal_row\tobjective_score\tpredicted_skelspheres\t"
                    + "exact_skelspheres\tpredicted_phesa_total\tpredicted_phesa_shape\t"
                    + "predicted_phesa_pharmacophore\n");
            for (var hit : hits) {
                out.write(hit.rank + "\t" + value(hit.compactRank) + "\t"
                        + tsv(hit.molecule.moleculeId()) + "\t" + tsv(hit.molecule.smiles()) + "\t"
                        + tsv(hit.molecule.canonicalSmiles()) + "\t" + hit.molecule.sourceRow() + "\t"
                        + hit.molecule.heavyAtomCount() + "\t" + hit.reference.shardIndex() + "\t"
                        + hit.reference.localRow() + "\t" + value(hit.objective) + "\t"
                        + value(hit.predictedSkelSpheres) + "\t" + value(hit.exactSkelSpheres) + "\t"
                        + value(hit.phesaTotal) + "\t" + value(hit.phesaShape) + "\t"
                        + value(hit.phesaPharmacophore) + "\n");
            }
        }
    }

    private static void writeSdf(Path path, List<RenderedHit> hits) throws Exception {
        var parser = new SmilesParser();
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            for (var hit : hits) {
                StereoMolecule molecule = new StereoMolecule();
                String smiles = hit.molecule.canonicalSmiles() == null
                        || hit.molecule.canonicalSmiles().isBlank()
                        ? hit.molecule.smiles() : hit.molecule.canonicalSmiles();
                parser.parse(molecule, smiles);
                molecule.setName(hit.molecule.moleculeId());
                new CoordinateInventor().invent(molecule);
                new MolfileCreator(molecule).writeMolfile(out);
                property(out, "RANK", Integer.toString(hit.rank));
                property(out, "MOLECULE_ID", hit.molecule.moleculeId());
                property(out, "SMILES", hit.molecule.smiles());
                property(out, "OBJECTIVE_SCORE", value(hit.objective));
                property(out, "PREDICTED_SKELSPHERES", value(hit.predictedSkelSpheres));
                property(out, "EXACT_SKELSPHERES", value(hit.exactSkelSpheres));
                property(out, "PREDICTED_PHESA_TOTAL", value(hit.phesaTotal));
                property(out, "PREDICTED_PHESA_SHAPE", value(hit.phesaShape));
                property(out, "PREDICTED_PHESA_PHARMACOPHORE", value(hit.phesaPharmacophore));
                out.write("$$$$\n");
            }
        }
    }

    private static void writeSummary(Path path, MoleculeSimilaritySearchConfig config,
            String artifact, long datasetRecords, long recordsScanned, List<RenderedHit> hits,
            Map<String, Object> timings) throws Exception {
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("# Hyperspace3D molecule similarity screen\n\n");
            out.write("- Mode: `" + config.mode() + "`\n");
            out.write("- Query: `" + config.query.identifier + "`\n");
            out.write("- Dataset artifact: `" + artifact + "`\n");
            out.write(String.format(Locale.ROOT, "- Dataset records: %,d; scanned: %,d\n",
                    datasetRecords, recordsScanned));
            out.write("- Objective: `" + objectiveLabel(config) + "`\n");
            out.write("- Timings (ns): `" + timings + "`\n\n");
            if (config.mode() == MoleculeSimilaritySearchConfig.Mode.SKELSPHERES_2D) {
                out.write("The learned 16D metric selected the shortlist, which was reranked "
                        + "with exact OCL SkelSpheres. Exact ranks apply within that shortlist.\n\n");
                out.write("| Rank | Molecule | Learned SkelSpheres | Exact SkelSpheres | SMILES |\n");
                out.write("|---:|---|---:|---:|---|\n");
                for (int index = 0; index < Math.min(config.output.reportTopK, hits.size()); index++) {
                    var hit = hits.get(index);
                    out.write("| " + hit.rank + " | " + hit.molecule.moleculeId() + " | "
                            + value(hit.predictedSkelSpheres) + " | "
                            + value(hit.exactSkelSpheres) + " | `"
                            + hit.molecule.smiles().replace("|", "\\|") + "` |\n");
                }
            } else {
                out.write("Learned PheSA scores are graph-only predictions; no conformers or "
                        + "exact PheSA were generated.\n\n");
                out.write("| Rank | Molecule | Objective | PheSA total | Shape | Pharmacophore | SMILES |\n");
                out.write("|---:|---|---:|---:|---:|---:|---|\n");
                for (int index = 0; index < Math.min(config.output.reportTopK, hits.size()); index++) {
                    var hit = hits.get(index);
                    out.write("| " + hit.rank + " | " + hit.molecule.moleculeId() + " | "
                            + value(hit.objective) + " | " + value(hit.phesaTotal) + " | "
                            + value(hit.phesaShape) + " | " + value(hit.phesaPharmacophore)
                            + " | `" + hit.molecule.smiles().replace("|", "\\|") + "` |\n");
                }
            }
        }
    }

    private static void writeRunManifest(ObjectMapper mapper,
            MoleculeSimilaritySearchConfig.ResolvedPaths paths,
            MoleculeSimilaritySearchConfig config, String artifact, long datasetRecords,
            long recordsScanned, DeepSpaceModelBundle model,
            CompactSkelSpheresModelBundle compact, StereoMolecule query,
            int resultCount, Map<String, Object> timings) throws Exception {
        Map<String, Object> run = new LinkedHashMap<>();
        run.put("artifactType", "hyperspace-molecule-similarity-search-run");
        run.put("formatVersion", 1);
        run.put("mode", config.mode().name().toLowerCase(Locale.ROOT));
        run.put("queryIdentifier", config.query.identifier);
        run.put("queryIdcode", new Canonizer(query).getIDCode());
        run.put("querySha256", sha256(new Canonizer(query).getIDCode().getBytes(StandardCharsets.UTF_8)));
        run.put("index", paths.index().toString());
        run.put("indexArtifactType", artifact);
        run.put("indexManifestSha256", sha256(paths.index().resolve("manifest.json")));
        run.put("datasetRecords", datasetRecords);
        run.put("recordsScanned", recordsScanned);
        run.put("modelBundleHash", model.bundleHash());
        if (compact != null) run.put("compactBundleHash", compact.bundleHash());
        run.put("objective", objectiveProvenance(config));
        run.put("resultCount", resultCount);
        run.put("timings", timings);
        run.put("outputs", Map.of("hitsTsv", paths.hitsTsv().toString(),
                "hitsSdf", paths.hitsSdf().toString(),
                "summaryMarkdown", paths.summaryMarkdown().toString()));
        mapper.writeValue(paths.runManifest().toFile(), run);
    }

    private static StereoMolecule parse(String structure, String format) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        if ("smiles".equalsIgnoreCase(format)) new SmilesParser().parse(molecule, structure);
        else new IDCodeParser().parse(molecule, structure);
        return molecule;
    }
    private static void createParents(MoleculeSimilaritySearchConfig.ResolvedPaths paths)
            throws Exception {
        for (Path path : List.of(paths.hitsTsv(), paths.hitsSdf(),
                paths.summaryMarkdown(), paths.runManifest())) {
            if (path.getParent() != null) Files.createDirectories(path.getParent());
        }
    }
    private static void property(BufferedWriter out, String name, String value) throws Exception {
        if (value == null || value.isEmpty()) return;
        out.write("> <" + name + ">\n" + value + "\n\n");
    }
    private static String tsv(String value) {
        return value == null ? "" : value.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
    }
    private static String value(Number value) {
        return value == null ? "" : String.format(Locale.ROOT, "%.8f", value.doubleValue());
    }
    private static String objectiveLabel(MoleculeSimilaritySearchConfig config) {
        if (config.mode() == MoleculeSimilaritySearchConfig.Mode.SKELSPHERES_2D) {
            return "direct:skelspheres_similarity";
        }
        return config.screening.objective.type + ":" + config.screening.objective.target;
    }
    private static Object objectiveProvenance(MoleculeSimilaritySearchConfig config) {
        if (config.mode() == MoleculeSimilaritySearchConfig.Mode.SKELSPHERES_2D) {
            return Map.of("type", "direct", "target", "skelspheres_similarity");
        }
        return config.screening.objective;
    }
    private static String sha256(Path path) throws Exception {
        try (InputStream input = Files.newInputStream(path)) {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] buffer = new byte[1024 * 1024];
            for (int count; (count = input.read(buffer)) >= 0;) digest.update(buffer, 0, count);
            return HexFormat.of().formatHex(digest.digest());
        }
    }
    private static String sha256(byte[] value) throws Exception {
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value));
    }
    private static Path config(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        throw new IllegalArgumentException("Usage: MoleculeSimilaritySearchCLI --config <search.json>");
    }

    private record RenderedHit(int rank, Integer compactRank,
            com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference reference,
            com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintMetadata molecule,
            Double objective, Double predictedSkelSpheres, Double exactSkelSpheres,
            Double phesaTotal, Double phesaShape, Double phesaPharmacophore) {}
}
