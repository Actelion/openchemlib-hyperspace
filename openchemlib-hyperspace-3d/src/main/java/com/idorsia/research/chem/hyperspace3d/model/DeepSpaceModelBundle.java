package com.idorsia.research.chem.hyperspace3d.model;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HexFormat;

public record DeepSpaceModelBundle(Path directory, DeepSpaceModelManifest manifest,
                                   String bundleHash) {
    private static final ObjectMapper MAPPER = new ObjectMapper();

    public static DeepSpaceModelBundle load(Path directory) {
        try {
            Path manifestPath = required(directory, "manifest.json");
            Path schemaPath = required(directory, "feature-schema.json");
            Path encoderPath = required(directory, "encoder.onnx");
            Path comparatorPath = required(directory, "comparator.onnx");
            required(directory, "checksums.sha256");
            DeepSpaceModelManifest manifest =
                    MAPPER.readValue(manifestPath.toFile(), DeepSpaceModelManifest.class);
            manifest.validate();
            verify(encoderPath, manifest.encoderSha256);
            verify(comparatorPath, manifest.comparatorSha256);
            verify(schemaPath, manifest.featureSchemaSha256);
            return new DeepSpaceModelBundle(directory.toAbsolutePath(), manifest,
                    sha256(manifestPath) + ":" + manifest.encoderSha256 + ":" + manifest.comparatorSha256);
        } catch (IOException | IllegalArgumentException e) {
            throw new DeepSpaceInferenceException("invalid model bundle at " + directory, e);
        }
    }

    public Path encoderPath() { return directory.resolve("encoder.onnx"); }
    public Path comparatorPath() { return directory.resolve("comparator.onnx"); }

    private static Path required(Path directory, String name) throws IOException {
        Path path = directory.resolve(name);
        if (!Files.isRegularFile(path)) throw new IOException("missing " + name);
        return path;
    }

    private static void verify(Path path, String expected) throws IOException {
        String actual = sha256(path);
        if (!actual.equals(expected)) throw new IOException("SHA-256 mismatch for " + path.getFileName());
    }

    private static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int n; (n = input.read(buffer)) >= 0;) digest.update(buffer, 0, n);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }
}
