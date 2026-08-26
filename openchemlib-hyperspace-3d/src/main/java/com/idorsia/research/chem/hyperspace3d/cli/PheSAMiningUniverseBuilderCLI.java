package com.idorsia.research.chem.hyperspace3d.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.index.MoleculeFingerprintDataSources;
import com.idorsia.research.chem.hyperspace3d.mining.MiningCandidateUniverse;
import com.idorsia.research.chem.hyperspace3d.model.DeepSpaceModelBundle;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.List;
import java.util.Set;

/** Builds the reusable deterministic candidate subset consumed by the query miner. */
public final class PheSAMiningUniverseBuilderCLI {
    private PheSAMiningUniverseBuilderCLI() {}
    public static void main(String[] args) throws Exception {
        if (args.length != 2 || !"--config".equals(args[0]))
            throw new IllegalArgumentException("usage: --config FILE");
        Path configPath = Path.of(args[1]).toAbsolutePath().normalize();
        Config config = new ObjectMapper().readValue(configPath.toFile(), Config.class);
        config.validate(); Path base = configPath.getParent();
        Path index = resolve(base, config.index); Path bundlePath = resolve(base, config.modelBundle);
        Path output = resolve(base, config.output); Set<String> exclusions = new HashSet<>();
        for (String file : config.exclusionFiles) {
            for (String line : Files.readAllLines(resolve(base, file), StandardCharsets.UTF_8))
                if (!line.isBlank() && !line.startsWith("#")) exclusions.add(line.trim());
        }
        var bundle = DeepSpaceModelBundle.load(bundlePath);
        String modelHash = sha256(bundle.bundleHash().getBytes(StandardCharsets.UTF_8));
        String sourceHash = sha256(Files.readAllBytes(index.resolve("manifest.json")));
        try (var source = MoleculeFingerprintDataSources.open(index)) {
            source.validateCompatibility(bundle, null, false);
            MiningCandidateUniverse.build(output, source, config.candidateCount,
                    config.exclusionOversample, config.seed, exclusions, sourceHash, modelHash);
        }
        System.out.printf("wrote %,d deterministic mining candidates to %s%n",
                config.candidateCount, output);
    }
    private static Path resolve(Path base, String value) { Path path = Path.of(value); return (path.isAbsolute() ? path : base.resolve(path)).normalize(); }
    private static String sha256(byte[] value) throws Exception { return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value)); }
    public static final class Config {
        public int formatVersion = 1; public String index; public String modelBundle; public String output;
        public int candidateCount = 1_000_000; public int exclusionOversample = 20_000;
        public long seed = 17; public List<String> exclusionFiles = List.of();
        public void validate() {
            if (formatVersion != 1 || index == null || modelBundle == null || output == null
                    || candidateCount < 1 || exclusionOversample < 0 || exclusionFiles == null)
                throw new IllegalArgumentException("invalid mining-universe config");
        }
    }
}
