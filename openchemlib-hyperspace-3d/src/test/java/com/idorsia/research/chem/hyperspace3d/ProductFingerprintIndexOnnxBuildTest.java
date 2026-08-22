package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.index.ProductFingerprintIndexReader;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceV1Encoder;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class ProductFingerprintIndexOnnxBuildTest {
    @TempDir Path temporary;

    @Test void realOnnxEncoderBuildsReadableNormalizedFingerprint() throws Exception {
        String configuredBundle = System.getProperty("hyperspace3d.modelBundle");
        Assumptions.assumeTrue(configuredBundle != null && !configuredBundle.isBlank());
        DeepSpaceModelBundle bundle = DeepSpaceModelBundle.load(Path.of(configuredBundle));
        RawSynthonSpace space = RawSynthonSpace.builder("onnx-test")
                .addRawFragments("r", 0, List.of(synthon(0, "a", "[U]CCC")))
                .addRawFragments("r", 1, List.of(synthon(1, "b", "[U]CCO"))).build();
        ProductFingerprintIndexBuildConfig config = new ProductFingerprintIndexBuildConfig();
        config.inputs.rawFull = "full";
        config.inputs.rawDownsampled = "reduced";
        config.inputs.modelBundle = configuredBundle;
        config.output.directory = "index";
        config.output.targetRecordCount = 1;
        config.output.recordsPerShard = 1;
        config.runtime.cpuWorkers = 1;
        config.runtime.encoderBatchSize = 1;
        config.runtime.progressIntervalSeconds = 0;
        config.filters.maxRotatableBonds = 0;
        config.validate();
        var provenance = new ProductFingerprintIndexBuilder.Provenance(
                "onnx-test:1.0", "raw", "onnx-test:1.0", "reduced", bundle.bundleHash(),
                ProductFingerprintIndexBuilder.configurationHash(config));
        DeepSpaceOnnxEnvironment runtime = new DeepSpaceOnnxEnvironment(
                DeepSpaceOnnxEnvironment.Device.CPU);
        try (DeepSpaceV1Encoder encoder = new DeepSpaceV1Encoder(runtime, bundle)) {
            var result = new ProductFingerprintIndexBuilder(space, space, encoder,
                    config, provenance).build(temporary.resolve("index"));
            try (ProductFingerprintIndexReader reader = new ProductFingerprintIndexReader(
                    result.manifestPath().getParent().resolve("shard-00000.h3di"))) {
                float[] embedding = reader.readBatch(2).get(0).embedding();
                double norm = 0d;
                for (float value : embedding) norm += value * value;
                assertEquals(1.0, Math.sqrt(norm), 1e-4);
            }
        }
    }

    private static RawSynthon synthon(int position, String id, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule("r", position, id, molecule);
    }
}
