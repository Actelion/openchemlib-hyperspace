package com.idorsia.research.chem.hyperspace3d.cli;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.Locale;

/** Configuration for query-centric exact-PheSA mining. */
public final class PheSAQueryPairMiningConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs(); public Output output = new Output();
    public Runtime runtime = new Runtime(); public Mining mining = new Mining();

    public void validate() {
        require(formatVersion == 1 && text(inputs.candidateUniverse) && text(inputs.modelBundle)
                && text(inputs.queries) && text(output.dataset), "mining paths are incomplete");
        require(runtime.comparatorBatchSize > 0 && runtime.exactLabelThreads > 0
                && mining.maxConformers == 8 && mining.phesaPpWeight == 0.5
                && mining.round >= 0, "unsupported initial mining contract");
        device();
    }

    public Resolved resolve(Path config) {
        Path base = config.toAbsolutePath().normalize().getParent();
        Path dataset = path(base, output.dataset);
        Path work = text(output.workDirectory) ? path(base, output.workDirectory)
                : dataset.resolveSibling(dataset.getFileName() + ".work");
        return new Resolved(path(base, inputs.candidateUniverse), path(base, inputs.modelBundle),
                path(base, inputs.queries), dataset, work);
    }

    public DeepSpaceOnnxEnvironment.Device device() {
        return DeepSpaceOnnxEnvironment.Device.valueOf(runtime.device.toUpperCase(Locale.ROOT));
    }
    private static Path path(Path base, String value) {
        Path path = Path.of(value); return (path.isAbsolute() ? path : base.resolve(path)).normalize();
    }
    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean value, String message) {
        if (!value) throw new IllegalArgumentException(message);
    }

    public record Resolved(Path universe, Path modelBundle, Path queries, Path dataset, Path work) {}
    public static final class Inputs {
        public String candidateUniverse; public String modelBundle; public String queries;
    }
    public static final class Output { public String dataset; public String workDirectory; }
    public static final class Runtime {
        public String device = "CPU"; public int cudaDeviceId;
        public int comparatorBatchSize = 32768; public int exactLabelThreads = 8;
    }
    public static final class Mining {
        public long seed = 17; public int round; public int maxConformers = 8;
        public double phesaPpWeight = 0.5;
    }
}
