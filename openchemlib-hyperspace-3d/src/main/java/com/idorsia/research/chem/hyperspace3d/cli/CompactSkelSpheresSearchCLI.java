package com.idorsia.research.chem.hyperspace3d.cli;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatchBuilder;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresScorer;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import com.idorsia.research.chem.hyperspace3d.screening.CompactSkelSpheresHit;
import com.idorsia.research.chem.hyperspace3d.screening.CompactSkelSpheresScreener;
import com.idorsia.research.chem.hyperspace3d.screening.ExactSkelSpheresReranker;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Standalone learned 2D search. Hybrid 2D-to-3D handoff is intentionally a later workflow. */
public final class CompactSkelSpheresSearchCLI {
    private CompactSkelSpheresSearchCLI() {}
    public static void main(String[] args) throws Exception {
        long started = System.nanoTime(); Path configPath = config(args).toAbsolutePath().normalize();
        var mapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);
        var config = mapper.readValue(configPath.toFile(), CompactSkelSpheresSearchConfig.class);
        config.validate(); var paths = config.resolve(configPath);
        var source = DeepSpaceModelBundle.load(paths.modelBundle());
        var compact = CompactSkelSpheresModelBundle.load(paths.compactBundle(), source.manifest());
        var index = ProductVectorIndexReader.loadManifest(paths.index(), config.runtime.verifyIndexChecksums);
        String sourceHash = ProductFingerprintIndexBuilder.sha256(paths.modelBundle().resolve("manifest.json"));
        String compactHash = ProductFingerprintIndexBuilder.sha256(paths.compactBundle().resolve("manifest.json"));
        if (!sourceHash.equals(index.sourceModelBundleHash)
                || !compactHash.equals(index.projectionModelBundleHash))
            throw new IllegalArgumentException("search model bundles do not match the index provenance");
        StereoMolecule query = parse(config.query);
        var runtime = new DeepSpaceOnnxEnvironment(config.device(), config.runtime.cudaDeviceId);
        float[] queryVector;
        long queryStarted = System.nanoTime();
        try (var encoder = new DeepSpaceV1Encoder(runtime, source);
             var projector = new CompactSkelSpheresProjector(runtime, compact)) {
            var batch = new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(query));
            queryVector = projector.project(encoder.encode(batch)[0]);
        }
        long queryNanos = System.nanoTime() - queryStarted;
        int globalK = config.exactRerank.enabled ? config.exactRerank.globalShortlist : config.output.globalTopK;
        int reactionK = config.exactRerank.enabled ? config.exactRerank.perReactionShortlist : config.output.perReactionTopK;
        var screening = new CompactSkelSpheresScreener(new CompactSkelSpheresScorer(compact.manifest()))
                .screen(paths.index(), index, queryVector, globalK, reactionK, config.runtime.scanBatchSize);
        List<CompactSkelSpheresHit> global = screening.global();
        Map<String,List<CompactSkelSpheresHit>> reactions = screening.byReaction();
        long rerankNanos = 0;
        if (config.exactRerank.enabled) {
            long rerankStarted = System.nanoTime(); Map<String,CompactSkelSpheresHit> union = new LinkedHashMap<>();
            for (var hit : global) union.put(key(hit), hit);
            for (var hits : reactions.values()) for (var hit : hits) union.putIfAbsent(key(hit), hit);
            var reranked = new ExactSkelSpheresReranker().rerank(query, new ArrayList<>(union.values()),
                    RawSynthonSpaceIO.read(paths.rawFull()));
            global = reranked.stream().limit(config.output.globalTopK).toList();
            Map<String,List<CompactSkelSpheresHit>> exactByReaction = new LinkedHashMap<>();
            for (var hit : reranked) exactByReaction.computeIfAbsent(hit.tuple().reactionId(), ignored -> new ArrayList<>());
            for (var hit : reranked) {
                var list = exactByReaction.get(hit.tuple().reactionId());
                if (list.size() < config.output.perReactionTopK) list.add(hit);
            }
            reactions = exactByReaction; rerankNanos = System.nanoTime() - rerankStarted;
        } else {
            global = global.stream().limit(config.output.globalTopK).toList();
        }
        writeHits(paths.hitsTsv(), global, reactions);
        Map<String,Object> run = new LinkedHashMap<>();
        run.put("artifactType", "hyperspace-skelspheres16-search-run");
        run.put("formatVersion", 1); run.put("queryIdentifier", config.query.identifier);
        run.put("queryIdcode", new Canonizer(query).getIDCode());
        run.put("indexManifestSha256", ProductFingerprintIndexBuilder.sha256(paths.index().resolve("manifest.json")));
        run.put("sourceModelBundleHash", sourceHash); run.put("projectionModelBundleHash", compactHash);
        run.put("recordsScanned", screening.recordsScanned()); run.put("exactRerankEnabled", config.exactRerank.enabled);
        run.put("queryEncodingNanos", queryNanos); run.put("vectorScanNanos", screening.scanNanos());
        run.put("tupleResolutionNanos", screening.tupleResolutionNanos()); run.put("exactRerankNanos", rerankNanos);
        run.put("totalNanos", System.nanoTime() - started); mapper.writeValue(paths.runManifest().toFile(), run);
        System.out.printf("Scanned %,d products; wrote %s%n", screening.recordsScanned(), paths.hitsTsv());
    }
    private static StereoMolecule parse(CompactSkelSpheresSearchConfig.Query query) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        if ("smiles".equalsIgnoreCase(query.format)) new SmilesParser().parse(molecule, query.structure);
        else new IDCodeParser().parse(molecule, query.structure);
        return molecule;
    }
    private static void writeHits(Path path, List<CompactSkelSpheresHit> global,
            Map<String,List<CompactSkelSpheresHit>> reactions) throws Exception {
        if (path.getParent() != null) Files.createDirectories(path.getParent());
        try (BufferedWriter out = Files.newBufferedWriter(path)) {
            out.write("source\treaction_id\tsynthon_ids\tsynthon_ordinals\tdot_product\tpredicted_similarity\texact_skelspheres\n");
            for (var hit : global) write(out, "global", hit);
            for (var entry : reactions.entrySet()) for (var hit : entry.getValue()) write(out, "reaction", hit);
        }
    }
    private static void write(BufferedWriter out, String source, CompactSkelSpheresHit hit) throws Exception {
        out.write(source + "\t" + hit.tuple().reactionId() + "\t"
                + String.join(",", hit.tuple().synthonIds()) + "\t"
                + hit.tuple().synthonOrdinals().stream().map(String::valueOf).collect(java.util.stream.Collectors.joining(","))
                + "\t" + hit.dotProduct() + "\t" + hit.predictedSimilarity() + "\t"
                + (hit.exactSimilarity() == null ? "" : hit.exactSimilarity()) + "\n");
    }
    private static String key(CompactSkelSpheresHit hit) {
        return hit.tuple().reactionId() + "\u0000" + String.join("\u0000", hit.tuple().synthonIds());
    }
    private static Path config(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        throw new IllegalArgumentException("Usage: CompactSkelSpheresSearchCLI --config <search-2d.json>");
    }
}
