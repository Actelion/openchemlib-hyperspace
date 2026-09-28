package com.idorsia.research.chem.hyperspace3d.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import com.idorsia.research.chem.hyperspace3d.index.build.CompactProductVectorIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.CompactProductVectorIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresBatchEncoder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.nio.file.Path;

/** Builds the columnar FP16 index directly from assembled product graphs. */
public final class CompactSkelSpheresIndexBuilderCLI {
    private CompactSkelSpheresIndexBuilderCLI() {}
    public static void main(String[] args) throws Exception {
        Path configPath = config(args).toAbsolutePath().normalize();
        var config = new ObjectMapper().readValue(configPath.toFile(),
                CompactProductVectorIndexBuildConfig.class);
        config.validate(); var paths = config.resolve(configPath);
        var full = RawSynthonSpaceIO.read(paths.rawFull());
        var downsampled = RawSynthonSpaceIO.read(paths.rawDownsampled());
        var source = DeepSpaceModelBundle.load(paths.modelBundle());
        var compact = CompactSkelSpheresModelBundle.load(paths.compactBundle(), source.manifest());
        String configurationHash = CompactProductVectorIndexBuilder.configurationHash(config);
        var provenance = new CompactProductVectorIndexBuilder.Provenance(
                ProductFingerprintIndexBuilder.identity(full),
                ProductFingerprintIndexBuilder.sha256(paths.rawFull()),
                ProductFingerprintIndexBuilder.identity(downsampled),
                ProductFingerprintIndexBuilder.sha256(paths.rawDownsampled()),
                ProductFingerprintIndexBuilder.sha256(paths.modelBundle().resolve("manifest.json")),
                ProductFingerprintIndexBuilder.sha256(paths.compactBundle().resolve("manifest.json")),
                compact.manifest().modelVersion, configurationHash);
        var runtime = new DeepSpaceOnnxEnvironment(config.device(), config.runtime.cudaDeviceId);
        try (var base = new DeepSpaceV1Encoder(runtime, source);
             var projection = new CompactSkelSpheresProjector(runtime, compact)) {
            var encoder = new CompactSkelSpheresBatchEncoder(base, projection);
            var result = new CompactProductVectorIndexBuilder(full, downsampled, encoder,
                    config, provenance).build(paths.outputDirectory());
            System.out.printf("Compact index complete: %,d records in %,d shards%nManifest: %s%n",
                    result.recordCount(), result.shardCount(), result.manifestPath());
        }
    }
    private static Path config(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        throw new IllegalArgumentException("Usage: CompactSkelSpheresIndexBuilderCLI --config <build-2d-index.json>");
    }
}
