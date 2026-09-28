package com.idorsia.research.chem.hyperspace3d;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.CompactSkelSpheresProjector;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceDualFingerprintBatchEncoder;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class MoleculeFingerprintIndexOnnxBuildTest {
    @TempDir Path temporary;

    @Test void realDualOnnxEncoderBuildsReadableFp16Vectors() throws Exception {
        String primaryPath = System.getProperty("hyperspace3d.modelBundle");
        String compactPath = System.getProperty("hyperspace3d.compactBundle");
        Assumptions.assumeTrue(primaryPath != null && compactPath != null,
                "set primary and compact bundle properties to run flat-index ONNX integration");
        DeepSpaceOnnxEnvironment.Device device =
                DeepSpaceOnnxEnvironment.Device.valueOf(
                        System.getProperty("hyperspace3d.testDevice", "CPU"));

        Path input = temporary.resolve("molecules.tsv");
        Files.writeString(input, "smiles\tid\nCCCCCCC\tmol-a\nc1ccccc1\tmol-b\n");
        Path output = temporary.resolve("index");
        MoleculeFingerprintIndexBuildConfig config =
                new MoleculeFingerprintIndexBuildConfig();
        config.inputs.library = input.toString();
        config.inputs.modelBundle = primaryPath;
        config.inputs.compactBundle = compactPath;
        config.output.directory = output.toString();
        config.output.sourceRowsPerShard = 1;
        config.runtime.device = device.name();
        config.runtime.cpuWorkers = 2;
        config.runtime.encoderBatchSize = 2;
        config.runtime.queueCapacity = 2;
        config.runtime.progressIntervalSeconds = 0;
        config.validate();
        Path syntheticConfig = temporary.resolve("build.json");
        var paths = config.resolve(syntheticConfig);

        DeepSpaceModelBundle primary = DeepSpaceModelBundle.load(Path.of(primaryPath));
        CompactSkelSpheresModelBundle compact = CompactSkelSpheresModelBundle.load(
                Path.of(compactPath), primary.manifest());
        DeepSpaceOnnxEnvironment runtime = new DeepSpaceOnnxEnvironment(device);
        try (DeepSpaceV1Encoder base = new DeepSpaceV1Encoder(runtime, primary);
             CompactSkelSpheresProjector projection =
                     new CompactSkelSpheresProjector(runtime, compact)) {
            var result = new MoleculeFingerprintIndexBuilder(
                    new DeepSpaceDualFingerprintBatchEncoder(base, projection),
                    config, paths).build();
            assertEquals(2, result.recordCount());
        }

        var manifest = MoleculeFingerprintIndexReader.loadManifest(output);
        for (int shard = 0; shard < manifest.shards.size(); shard++) {
            try (MoleculeFingerprintIndexReader reader =
                         new MoleculeFingerprintIndexReader(output, manifest, shard)) {
                var record = reader.read(0);
                assertEquals(128, record.base128().length);
                assertEquals(16, record.compact16().length);
                assertEquals(1.0, norm(record.base128()), 0.002);
                assertEquals(1.0, norm(record.compact16()), 0.002);
            }
        }
    }

    private static double norm(float[] vector) {
        double squared = 0;
        for (float value : vector) squared += value * value;
        return Math.sqrt(squared);
    }
}
