package com.idorsia.research.chem.hyperspace3d.index.build;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

/** JSON configuration for streaming a flat molecule archive into a dual-vector index. */
public final class MoleculeFingerprintIndexBuildConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Output output = new Output();
    public Runtime runtime = new Runtime();

    public void validate() {
        if (runtime.encoderBatchSize == 0) {
            runtime.encoderBatchSize =
                    device() == DeepSpaceOnnxEnvironment.Device.CUDA ? 512 : 64;
        }
        require(formatVersion == 1, "formatVersion must be 1");
        require(text(inputs.library), "inputs.library is required");
        require(text(inputs.smilesColumn), "inputs.smilesColumn is required");
        require(text(inputs.idColumn), "inputs.idColumn is required");
        require(!inputs.smilesColumn.equals(inputs.idColumn),
                "SMILES and ID columns must differ");
        require(text(inputs.modelBundle), "inputs.modelBundle is required");
        require(text(inputs.compactBundle), "inputs.compactBundle is required");
        require(text(output.directory), "output.directory is required");
        require(output.sourceRowsPerShard > 0
                        && output.sourceRowsPerShard <= 8_000_000,
                "output.sourceRowsPerShard must be within 1..8,000,000");
        require(runtime.cpuWorkers > 0, "runtime.cpuWorkers must be positive");
        require(runtime.encoderBatchSize > 0,
                "runtime.encoderBatchSize must be positive");
        require(runtime.queueCapacity > 0, "runtime.queueCapacity must be positive");
        require(runtime.progressIntervalSeconds >= 0,
                "runtime.progressIntervalSeconds must be non-negative");
        device();
    }

    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().normalize().getParent();
        return new ResolvedPaths(resolve(base, inputs.library),
                resolve(base, inputs.modelBundle),
                resolve(base, inputs.compactBundle),
                resolve(base, output.directory));
    }

    @JsonIgnore public DeepSpaceOnnxEnvironment.Device device() {
        try {
            return DeepSpaceOnnxEnvironment.Device.valueOf(
                    runtime.device.toUpperCase(Locale.ROOT));
        } catch (RuntimeException error) {
            throw new IllegalArgumentException("runtime.device must be CPU or CUDA", error);
        }
    }

    @JsonIgnore public Map<String, Object> runtimeDescription() {
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("device", runtime.device.toUpperCase(Locale.ROOT));
        result.put("cudaDeviceId", runtime.cudaDeviceId);
        result.put("cpuWorkers", runtime.cpuWorkers);
        result.put("encoderBatchSize", runtime.encoderBatchSize);
        result.put("queueCapacity", runtime.queueCapacity);
        return result;
    }

    private static Path resolve(Path base, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : base.resolve(path))
                .normalize().toAbsolutePath();
    }
    private static boolean text(String value) {
        return value != null && !value.isBlank();
    }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    public record ResolvedPaths(Path library, Path modelBundle, Path compactBundle,
                                Path outputDirectory) {}

    public static final class Inputs {
        public String library;
        public String smilesColumn = "smiles";
        public String idColumn = "id";
        public String modelBundle;
        public String compactBundle;
    }

    public static final class Output {
        public String directory;
        public long sourceRowsPerShard = 1_000_000;
        public boolean resume = true;
    }

    public static final class Runtime {
        public String device = "CUDA";
        public int cudaDeviceId;
        public int cpuWorkers = Math.max(1,
                java.lang.Runtime.getRuntime().availableProcessors() - 2);
        public int encoderBatchSize;
        public int queueCapacity = 4;
        public int progressIntervalSeconds = 30;
    }
}
