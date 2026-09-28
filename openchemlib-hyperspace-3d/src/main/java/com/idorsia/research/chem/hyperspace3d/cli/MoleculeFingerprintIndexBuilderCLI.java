package com.idorsia.research.chem.hyperspace3d.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceDualFingerprintBatchEncoder;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.nio.file.Path;

/** Builds a resumable 128D/16D FP16 index from a flat molecule table. */
public final class MoleculeFingerprintIndexBuilderCLI {
    private MoleculeFingerprintIndexBuilderCLI() {}

    public static void main(String[] args) throws Exception {
        Path configPath = config(args).toAbsolutePath().normalize();
        MoleculeFingerprintIndexBuildConfig config = new ObjectMapper().readValue(
                configPath.toFile(), MoleculeFingerprintIndexBuildConfig.class);
        config.validate();
        var paths = config.resolve(configPath);
        DeepSpaceModelBundle source = DeepSpaceModelBundle.load(paths.modelBundle());
        CompactSkelSpheresModelBundle compact = CompactSkelSpheresModelBundle.load(
                paths.compactBundle(), source.manifest());
        DeepSpaceOnnxEnvironment runtime =
                new DeepSpaceOnnxEnvironment(config.device(), config.runtime.cudaDeviceId);
        try (DeepSpaceV1Encoder base = new DeepSpaceV1Encoder(runtime, source);
             CompactSkelSpheresProjector projection =
                     new CompactSkelSpheresProjector(runtime, compact)) {
            var encoder = new DeepSpaceDualFingerprintBatchEncoder(base, projection);
            var result = new MoleculeFingerprintIndexBuilder(
                    encoder, config, paths).build();
            System.out.printf(
                    "Molecule index complete: %,d/%,d accepted rows in %,d shards%nManifest: %s%n",
                    result.recordCount(), result.sourceRowCount(), result.shardCount(),
                    result.manifestPath());
        }
    }

    private static Path config(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) return Path.of(args[0]);
        if (args.length == 2 && "--config".equals(args[0])) return Path.of(args[1]);
        throw new IllegalArgumentException(
                "Usage: MoleculeFingerprintIndexBuilderCLI --config <build-molecule-index.json>");
    }
}
