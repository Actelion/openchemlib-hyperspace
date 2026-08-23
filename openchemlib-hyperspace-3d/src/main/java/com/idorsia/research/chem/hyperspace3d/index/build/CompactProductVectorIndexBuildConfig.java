package com.idorsia.research.chem.hyperspace3d.index.build;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.Locale;

/** Configuration for direct graph-to-16D index construction. */
public final class CompactProductVectorIndexBuildConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Output output = new Output();
    public ProductFingerprintIndexBuildConfig.Sampling sampling = new ProductFingerprintIndexBuildConfig.Sampling();
    public ProductFingerprintIndexBuildConfig.Filters filters = new ProductFingerprintIndexBuildConfig.Filters();
    public ProductFingerprintIndexBuildConfig.Runtime runtime = new ProductFingerprintIndexBuildConfig.Runtime();

    public void validate() {
        if (runtime.encoderBatchSize == 0) runtime.encoderBatchSize = device() == DeepSpaceOnnxEnvironment.Device.CUDA ? 512 : 64;
        require(formatVersion == 1, "formatVersion must be 1");
        require(text(inputs.rawFull) && text(inputs.rawDownsampled)
                && text(inputs.modelBundle) && text(inputs.compactBundle), "all input paths are required");
        require(text(output.directory) && output.targetRecordCount > 0 && output.recordsPerShard > 0,
                "valid output directory and record counts are required");
        require(sampling.minimumAcceptedPerReaction >= 0 && sampling.coverageAttemptsPerReaction > 0
                && sampling.maxAttemptsMultiplier >= 1, "invalid sampling limits");
        sampling.reactionWeighting.validate();
        require(filters.minHeavyAtoms >= 1 && filters.maxHeavyAtoms <= 32
                && filters.maxHeavyAtoms >= filters.minHeavyAtoms && filters.maxRotatableBonds >= 0,
                "invalid molecular filters");
        require(runtime.cpuWorkers > 0 && runtime.encoderBatchSize > 0
                && runtime.queueCapacity > 0, "invalid runtime settings");
    }

    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().getParent();
        return new ResolvedPaths(resolve(base, inputs.rawFull), resolve(base, inputs.rawDownsampled),
                resolve(base, inputs.modelBundle), resolve(base, inputs.compactBundle),
                resolve(base, output.directory));
    }
    @JsonIgnore public DeepSpaceOnnxEnvironment.Device device() {
        try { return DeepSpaceOnnxEnvironment.Device.valueOf(runtime.device.toUpperCase(Locale.ROOT)); }
        catch (RuntimeException error) { throw new IllegalArgumentException("runtime.device must be CPU or CUDA", error); }
    }
    private static Path resolve(Path base, String value) {
        Path path = Path.of(value); return (path.isAbsolute() ? path : base.resolve(path)).normalize().toAbsolutePath();
    }
    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean condition, String message) { if (!condition) throw new IllegalArgumentException(message); }
    public record ResolvedPaths(Path rawFull, Path rawDownsampled, Path modelBundle,
                                Path compactBundle, Path outputDirectory) {}
    public static final class Inputs {
        public String rawFull;
        public String rawDownsampled;
        public String modelBundle;
        public String compactBundle;
    }
    public static final class Output {
        public String directory;
        public long targetRecordCount;
        public long recordsPerShard = 1_000_000;
    }
}
