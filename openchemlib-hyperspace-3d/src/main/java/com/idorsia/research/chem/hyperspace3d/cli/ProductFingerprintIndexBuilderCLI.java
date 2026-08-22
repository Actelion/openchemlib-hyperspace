package com.idorsia.research.chem.hyperspace3d.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.nio.file.Path;

/** Command-line entry point for persistent sampled-product index construction. */
public final class ProductFingerprintIndexBuilderCLI {
    private ProductFingerprintIndexBuilderCLI() {}

    public static void main(String[] args) throws Exception {
        Path configPath = parseConfig(args).toAbsolutePath().normalize();
        ProductFingerprintIndexBuildConfig config = new ObjectMapper().readValue(
                configPath.toFile(), ProductFingerprintIndexBuildConfig.class);
        config.validate();
        ProductFingerprintIndexBuildConfig.ResolvedPaths paths = config.resolve(configPath);

        System.out.println("Loading full rawspace: " + paths.rawFull());
        RawSynthonSpace full = RawSynthonSpaceIO.read(paths.rawFull());
        System.out.println("Loading downsampled rawspace: " + paths.rawDownsampled());
        RawSynthonSpace downsampled = RawSynthonSpaceIO.read(paths.rawDownsampled());
        DeepSpaceModelBundle bundle = DeepSpaceModelBundle.load(paths.modelBundle());
        String rawHash = ProductFingerprintIndexBuilder.sha256(paths.rawFull());
        String downsampledHash = ProductFingerprintIndexBuilder.sha256(paths.rawDownsampled());
        String configHash = ProductFingerprintIndexBuilder.configurationHash(config);
        var provenance = new ProductFingerprintIndexBuilder.Provenance(
                ProductFingerprintIndexBuilder.identity(full), rawHash,
                ProductFingerprintIndexBuilder.identity(downsampled), downsampledHash,
                bundle.bundleHash(), configHash);

        System.out.printf("Building %,d records with %s, batch %,d, %,d CPU workers%n",
                config.output.targetRecordCount, config.device(), config.runtime.encoderBatchSize,
                config.runtime.cpuWorkers);
        DeepSpaceOnnxEnvironment runtime = new DeepSpaceOnnxEnvironment(
                config.device(), config.runtime.cudaDeviceId);
        try (DeepSpaceV1Encoder encoder = new DeepSpaceV1Encoder(runtime, bundle)) {
            var result = new ProductFingerprintIndexBuilder(full, downsampled, encoder,
                    config, provenance).build(paths.outputDirectory());
            System.out.printf("Index complete: %,d records in %,d shards%nManifest: %s%n",
                    result.recordCount(), result.shardCount(), result.manifestPath());
        }
    }

    private static Path parseConfig(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        if (args.length == 1 && ("--help".equals(args[0]) || "-h".equals(args[0]))) {
            System.out.println("Usage: ProductFingerprintIndexBuilderCLI --config <build-index.json>");
            System.exit(0);
        }
        throw new IllegalArgumentException(
                "Usage: ProductFingerprintIndexBuilderCLI --config <build-index.json>");
    }
}
