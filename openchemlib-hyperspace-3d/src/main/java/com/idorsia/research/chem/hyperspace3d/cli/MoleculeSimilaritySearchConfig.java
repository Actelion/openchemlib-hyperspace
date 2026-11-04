package com.idorsia.research.chem.hyperspace3d.cli;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.screening.ScreeningObjective;
import java.nio.file.Path;
import java.util.Locale;
import java.util.Set;

/** Configuration for flat-molecule learned 2D, learned 3D, and cascade screens. */
public final class MoleculeSimilaritySearchConfig {
    private static final Set<String> PHESA_TARGETS = Set.of(
            "phesa_total", "phesa_shape", "phesa_pharmacophore");
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Query query = new Query();
    public Screening screening = new Screening();
    public Output output = new Output();
    public Runtime runtime = new Runtime();

    public enum Mode { SKELSPHERES_2D, PHESA_3D, SKELSPHERES_THEN_PHESA, POSE_PHESA_3D }

    public void validate() {
        require(formatVersion == 1, "formatVersion must be 1");
        require(text(inputs.index) && text(inputs.modelBundle), "index and modelBundle are required");
        if (usesCompact()) require(text(inputs.compactBundle), "compactBundle is required for this mode");
        require(text(query.identifier), "query identifier is required");
        if (usesPose()) {
            require(text(query.file) && text(inputs.poseBundle), "pose query file and poseBundle are required");
            require(query.record == null || query.record >= 0, "query record must be nonnegative (zero based)");
            require(Double.isFinite(screening.posePpWeight) && screening.posePpWeight >= 0
                    && screening.posePpWeight <= 1, "posePpWeight must be in [0,1]");
        } else {
            require(text(query.structure), "query structure is required");
            require("smiles".equalsIgnoreCase(query.format) || "idcode".equalsIgnoreCase(query.format),
                    "query format must be smiles or idcode");
        }
        mode();
        objective();
        require(screening.resultTopK > 0, "resultTopK must be positive");
        require(screening.cascadeShortlistTopK >= screening.resultTopK,
                "cascadeShortlistTopK must be at least resultTopK");
        require(output.reportTopK > 0 && output.reportTopK <= screening.resultTopK,
                "reportTopK must be within 1..resultTopK");
        require(text(output.hitsTsv) && text(output.hitsSdf)
                        && text(output.summaryMarkdown) && text(output.runManifest),
                "all output paths are required");
        require(runtime.scanBatchSize > 0 && runtime.comparatorBatchSize > 0
                        && runtime.prefetchDepth > 0 && runtime.intraOpThreads >= 0,
                "batch sizes must be positive");
        device();
    }

    public Mode mode() {
        try { return Mode.valueOf(screening.mode.toUpperCase(Locale.ROOT)); }
        catch (RuntimeException error) {
            throw new IllegalArgumentException("mode must be skelspheres_2d, phesa_3d, skelspheres_then_phesa, or pose_phesa_3d", error);
        }
    }

    public boolean usesPose() { return mode() == Mode.POSE_PHESA_3D; }
    public boolean usesCompact() { return mode() == Mode.SKELSPHERES_2D || mode() == Mode.SKELSPHERES_THEN_PHESA; }
    public boolean usesComparator() { return mode() != Mode.SKELSPHERES_2D; }

    public ScreeningObjective objective() {
        if (!usesComparator()) return ScreeningObjective.direct("skelspheres_similarity");
        if ("direct".equalsIgnoreCase(screening.objective.type)) {
            require(PHESA_TARGETS.contains(screening.objective.target),
                    "direct 3D target must be phesa_total, phesa_shape, or phesa_pharmacophore");
            return ScreeningObjective.direct(screening.objective.target);
        }
        if ("composite".equalsIgnoreCase(screening.objective.type)) {
            require(Double.isFinite(screening.objective.shapeWeight)
                            && Double.isFinite(screening.objective.pharmacophoreWeight),
                    "composite weights must be finite");
            return ScreeningObjective.phesaComposite(screening.objective.shapeWeight,
                    screening.objective.pharmacophoreWeight);
        }
        throw new IllegalArgumentException("objective type must be direct or composite");
    }

    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().normalize().getParent();
        return new ResolvedPaths(resolve(base, inputs.index), resolve(base, inputs.modelBundle),
                text(inputs.compactBundle) ? resolve(base, inputs.compactBundle) : null,
                resolve(base, output.hitsTsv), resolve(base, output.hitsSdf),
                resolve(base, output.summaryMarkdown), resolve(base, output.runManifest));
    }

    public DeepSpaceOnnxEnvironment.Device device() {
        try { return DeepSpaceOnnxEnvironment.Device.valueOf(runtime.device.toUpperCase(Locale.ROOT)); }
        catch (RuntimeException error) { throw new IllegalArgumentException("device must be CPU or CUDA", error); }
    }

    public Path poseBundlePath(Path config) { return resolve(config.toAbsolutePath().getParent(), inputs.poseBundle); }
    public Path queryFilePath(Path config) { return resolve(config.toAbsolutePath().getParent(), query.file); }

    private static Path resolve(Path base, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : base.resolve(path)).normalize().toAbsolutePath();
    }
    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    public record ResolvedPaths(Path index, Path modelBundle, Path compactBundle,
            Path hitsTsv, Path hitsSdf, Path summaryMarkdown, Path runManifest) {}
    public static final class Inputs {
        public String index;
        public String modelBundle;
        public String compactBundle;
        public String poseBundle;
    }
    public static final class Query {
        public String structure;
        public String format = "smiles";
        public String identifier;
        public String file;
        public Integer record;
    }
    public static final class Screening {
        public String mode = "phesa_3d";
        public double posePpWeight = 0.5;
        public Objective objective = new Objective();
        public int resultTopK = 1000;
        public int cascadeShortlistTopK = 100000;
    }
    public static final class Objective {
        public String type = "direct";
        public String target = "phesa_total";
        public double shapeWeight = 1.0;
        public double pharmacophoreWeight = 1.0;
    }
    public static final class Output {
        public String hitsTsv;
        public String hitsSdf;
        public String summaryMarkdown;
        public String runManifest;
        public int reportTopK = 25;
    }
    public static final class Runtime {
        public int intraOpThreads;
        public String device = "CPU";
        public int cudaDeviceId;
        public int scanBatchSize = 65536;
        public int comparatorBatchSize = 32768;
        public int prefetchDepth = 2;
    }
}
