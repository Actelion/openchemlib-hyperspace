package com.idorsia.research.chem.hyperspace3d.index;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/** Auto-detects supported flat-molecule fingerprint storage formats. */
public final class MoleculeFingerprintDataSources {
    private MoleculeFingerprintDataSources() {}

    public static MoleculeFingerprintDataSource open(Path directory) throws IOException {
        Path normalized = directory.toAbsolutePath().normalize();
        Path manifest = normalized.resolve("manifest.json");
        if (!Files.isRegularFile(manifest)) throw new IOException("missing index manifest: " + manifest);
        JsonNode root = new ObjectMapper().readTree(manifest.toFile());
        String artifact = root.path("artifactType").asText(root.path("artifact_type").asText());
        return switch (artifact) {
            case "hyperspace-molecule-fingerprint-index" ->
                    new JavaMoleculeFingerprintDataSource(normalized);
            case DeepSpace7SupplierCacheDataSource.ARTIFACT_TYPE ->
                    new DeepSpace7SupplierCacheDataSource(normalized);
            default -> throw new IOException("unsupported molecule fingerprint artifact: " + artifact);
        };
    }
}
