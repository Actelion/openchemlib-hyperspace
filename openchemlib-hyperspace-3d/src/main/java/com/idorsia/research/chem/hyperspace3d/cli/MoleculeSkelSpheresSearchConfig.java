package com.idorsia.research.chem.hyperspace3d.cli;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.Locale;

/** Configuration for exhaustive learned SkelSpheres search of a flat molecule index. */
public final class MoleculeSkelSpheresSearchConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Query query = new Query();
    public Output output = new Output();
    public Runtime runtime = new Runtime();

    public void validate() {
        require(formatVersion == 1, "formatVersion must be 1");
        require(text(inputs.index) && text(inputs.modelBundle)
                && text(inputs.compactBundle), "invalid search inputs");
        require(text(query.structure) && ("smiles".equalsIgnoreCase(query.format)
                || "idcode".equalsIgnoreCase(query.format)),
                "query format must be smiles or idcode");
        require(text(query.identifier), "query identifier is required");
        require(text(output.hitsTsv) && text(output.hitsSdf)
                && text(output.summaryMarkdown) && text(output.runManifest),
                "all output paths are required");
        require(output.learnedTopK > 0, "learnedTopK must be positive");
        require(output.reportTopK > 0 && output.reportTopK <= output.learnedTopK,
                "reportTopK must be within 1..learnedTopK");
        require(runtime.scanBatchSize > 0, "scanBatchSize must be positive");
        device();
    }

    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().normalize().getParent();
        return new ResolvedPaths(resolve(base, inputs.index),
                resolve(base, inputs.modelBundle), resolve(base, inputs.compactBundle),
                resolve(base, output.hitsTsv), resolve(base, output.hitsSdf),
                resolve(base, output.summaryMarkdown), resolve(base, output.runManifest));
    }

    public DeepSpaceOnnxEnvironment.Device device() {
        try {
            return DeepSpaceOnnxEnvironment.Device.valueOf(
                    runtime.device.toUpperCase(Locale.ROOT));
        } catch (RuntimeException error) {
            throw new IllegalArgumentException("device must be CPU or CUDA", error);
        }
    }

    private static Path resolve(Path base, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : base.resolve(path)).normalize().toAbsolutePath();
    }
    private static boolean text(String value) {
        return value != null && !value.isBlank();
    }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    public record ResolvedPaths(Path index, Path modelBundle, Path compactBundle,
            Path hitsTsv, Path hitsSdf, Path summaryMarkdown, Path runManifest) {}
    public static final class Inputs {
        public String index;
        public String modelBundle;
        public String compactBundle;
    }
    public static final class Query {
        public String structure;
        public String format = "smiles";
        public String identifier;
    }
    public static final class Output {
        public String hitsTsv;
        public String hitsSdf;
        public String summaryMarkdown;
        public String runManifest;
        public int learnedTopK = 1000;
        public int reportTopK = 25;
    }
    public static final class Runtime {
        public String device = "CPU";
        public int cudaDeviceId;
        public int scanBatchSize = 65536;
    }
}
