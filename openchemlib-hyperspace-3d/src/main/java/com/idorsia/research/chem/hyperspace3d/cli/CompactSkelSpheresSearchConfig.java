package com.idorsia.research.chem.hyperspace3d.cli;

import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceOnnxEnvironment;
import java.nio.file.Path;
import java.util.Locale;

public final class CompactSkelSpheresSearchConfig {
    public int formatVersion = 1;
    public Inputs inputs = new Inputs();
    public Query query = new Query();
    public Output output = new Output();
    public Runtime runtime = new Runtime();
    public ExactRerank exactRerank = new ExactRerank();
    public void validate() {
        require(formatVersion == 1 && text(inputs.index) && text(inputs.modelBundle)
                && text(inputs.compactBundle), "invalid search inputs");
        require(text(query.structure) && ("smiles".equalsIgnoreCase(query.format)
                || "idcode".equalsIgnoreCase(query.format)), "query format must be smiles or idcode");
        require(text(output.hitsTsv) && text(output.runManifest)
                && output.globalTopK > 0 && output.perReactionTopK >= 0, "invalid output settings");
        require(runtime.scanBatchSize > 0, "scanBatchSize must be positive");
        if (exactRerank.enabled) {
            require(text(inputs.rawFull), "rawFull is required for exact reranking");
            require(exactRerank.globalShortlist >= output.globalTopK
                    && exactRerank.perReactionShortlist >= output.perReactionTopK,
                    "exact rerank shortlists must cover requested outputs");
        }
        device();
    }
    public ResolvedPaths resolve(Path configPath) {
        Path base = configPath.toAbsolutePath().getParent();
        return new ResolvedPaths(resolve(base, inputs.index), resolve(base, inputs.modelBundle),
                resolve(base, inputs.compactBundle), text(inputs.rawFull) ? resolve(base, inputs.rawFull) : null,
                resolve(base, output.hitsTsv), resolve(base, output.runManifest));
    }
    public DeepSpaceOnnxEnvironment.Device device() {
        try { return DeepSpaceOnnxEnvironment.Device.valueOf(runtime.device.toUpperCase(Locale.ROOT)); }
        catch (RuntimeException error) { throw new IllegalArgumentException("device must be CPU or CUDA", error); }
    }
    private static Path resolve(Path base, String value) { Path p = Path.of(value); return (p.isAbsolute() ? p : base.resolve(p)).normalize().toAbsolutePath(); }
    private static boolean text(String value) { return value != null && !value.isBlank(); }
    private static void require(boolean value, String message) { if (!value) throw new IllegalArgumentException(message); }
    public record ResolvedPaths(Path index, Path modelBundle, Path compactBundle, Path rawFull,
                                Path hitsTsv, Path runManifest) {}
    public static final class Inputs { public String index; public String modelBundle; public String compactBundle; public String rawFull; }
    public static final class Query { public String structure; public String format = "smiles"; public String identifier; }
    public static final class Output { public String hitsTsv; public String runManifest; public int globalTopK = 1000; public int perReactionTopK = 50; }
    public static final class Runtime { public String device = "CPU"; public int cudaDeviceId; public int scanBatchSize = 65536; public boolean verifyIndexChecksums = true; }
    public static final class ExactRerank { public boolean enabled; public int globalShortlist = 10000; public int perReactionShortlist = 100; }
}
