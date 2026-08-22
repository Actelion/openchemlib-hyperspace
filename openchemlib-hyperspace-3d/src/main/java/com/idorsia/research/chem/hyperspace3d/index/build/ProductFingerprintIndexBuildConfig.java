package com.idorsia.research.chem.hyperspace3d.index.build;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.idorsia.research.chem.hyperspace.screening.ReactionScheduler;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/** Versioned JSON configuration for reproducible product-index builds. */
public final class ProductFingerprintIndexBuildConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Output output = new Output();
    public Sampling sampling = new Sampling();
    public Filters filters = new Filters();
    public Runtime runtime = new Runtime();

    public void validate() {
        if (runtime.encoderBatchSize == 0) runtime.encoderBatchSize = device() == DeepSpaceOnnxEnvironment.Device.CUDA ? 512 : 64;
        require(formatVersion == 1, "formatVersion must be 1");
        require(text(inputs.rawFull), "inputs.rawFull is required");
        require(text(inputs.rawDownsampled), "inputs.rawDownsampled is required");
        require(text(inputs.modelBundle), "inputs.modelBundle is required");
        require(text(output.directory), "output.directory is required");
        require(output.targetRecordCount > 0, "output.targetRecordCount must be positive");
        require(output.recordsPerShard > 0, "output.recordsPerShard must be positive");
        require(sampling.minimumAcceptedPerReaction >= 0,
                "sampling.minimumAcceptedPerReaction must be non-negative");
        require(sampling.coverageAttemptsPerReaction > 0,
                "sampling.coverageAttemptsPerReaction must be positive");
        require(sampling.maxAttemptsMultiplier >= 1,
                "sampling.maxAttemptsMultiplier must be at least 1");
        sampling.reactionWeighting.validate();
        require(filters.minHeavyAtoms >= 1 && filters.maxHeavyAtoms <= 32
                        && filters.maxHeavyAtoms >= filters.minHeavyAtoms,
                "filters heavy-atom interval must be within 1..32");
        require(filters.maxRotatableBonds >= 0,
                "filters.maxRotatableBonds must be non-negative");
        require(runtime.cpuWorkers > 0, "runtime.cpuWorkers must be positive");
        require(runtime.encoderBatchSize > 0, "runtime.encoderBatchSize must be positive");
        require(runtime.queueCapacity > 0, "runtime.queueCapacity must be positive");
        require(runtime.progressIntervalSeconds >= 0,
                "runtime.progressIntervalSeconds must be non-negative");
        device();
    }

    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().getParent();
        return new ResolvedPaths(resolve(base, inputs.rawFull), resolve(base, inputs.rawDownsampled),
                resolve(base, inputs.modelBundle), resolve(base, output.directory));
    }

    @JsonIgnore
    public DeepSpaceOnnxEnvironment.Device device() {
        try {
            return DeepSpaceOnnxEnvironment.Device.valueOf(runtime.device.toUpperCase(Locale.ROOT));
        } catch (RuntimeException error) {
            throw new IllegalArgumentException("runtime.device must be CPU or CUDA", error);
        }
    }

    private static Path resolve(Path base, String value) {
        Path path = Path.of(value);
        return (path.isAbsolute() ? path : base.resolve(path)).normalize().toAbsolutePath();
    }

    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean condition, String message) {
        if (!condition) throw new IllegalArgumentException(message);
    }

    public record ResolvedPaths(Path rawFull, Path rawDownsampled, Path modelBundle,
                                Path outputDirectory) {}

    public static final class Inputs {
        public String rawFull;
        public String rawDownsampled;
        public String modelBundle;
    }

    public static final class Output {
        public String directory;
        public long targetRecordCount;
        public long recordsPerShard = 1_000_000;
        public boolean resume = true;
    }

    public static final class Sampling {
        public long seed = 13;
        public int minimumAcceptedPerReaction = 1;
        public int coverageAttemptsPerReaction = 100;
        public long maxAttemptsMultiplier = 20;
        public ReactionWeighting reactionWeighting = new ReactionWeighting();
    }

    public static final class ReactionWeighting {
        public String mode = ReactionScheduler.Mode.EXPONENT.name();
        public double minWeight = 0.01;
        public double exponent = 1.0;
        public List<Bucket> buckets = new ArrayList<>();

        public void validate() { toWeighting(); }

        @JsonIgnore
        public ReactionScheduler.Weighting toWeighting() {
            ReactionScheduler.Mode parsed;
            try {
                parsed = ReactionScheduler.Mode.valueOf(mode.toUpperCase(Locale.ROOT));
            } catch (RuntimeException error) {
                throw new IllegalArgumentException("sampling.reactionWeighting.mode is invalid", error);
            }
            if (parsed == ReactionScheduler.Mode.EXPONENT) {
                return ReactionScheduler.Weighting.exponent(minWeight, exponent);
            }
            List<ReactionScheduler.Bucket> converted = new ArrayList<>();
            for (Bucket bucket : buckets) {
                if (bucket == null) throw new IllegalArgumentException("reaction weight bucket is null");
                converted.add(new ReactionScheduler.Bucket(bucket.maxProductExclusive, bucket.weight));
            }
            return ReactionScheduler.Weighting.bucketedProduct(converted);
        }
    }

    public static final class Bucket {
        public Double maxProductExclusive;
        public double weight;
    }

    public static final class Filters {
        public int minHeavyAtoms = 6;
        public int maxHeavyAtoms = 32;
        /** Zero disables this additional filter. */
        public int maxRotatableBonds = 15;
    }

    public static final class Runtime {
        public String device = "CPU";
        public int cudaDeviceId = 0;
        public int cpuWorkers = Math.max(1, java.lang.Runtime.getRuntime().availableProcessors() - 1);
        /** Zero in JSON is replaced by the device-specific default by normalizeDefaults(). */
        public int encoderBatchSize = 0;
        public int queueCapacity = 4;
        public int progressIntervalSeconds = 30;
    }
}
