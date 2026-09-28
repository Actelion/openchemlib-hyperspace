package com.idorsia.research.chem.hyperspace3d.model;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HexFormat;

public record CompactSkelSpheresModelBundle(
        Path directory,
        CompactSkelSpheresManifest manifest,
        Path projectionPath,
        String bundleHash) {

    public static CompactSkelSpheresModelBundle load(
            Path directory, DeepSpaceModelManifest sourceManifest) throws IOException {
        Path manifestPath = required(directory, "manifest.json");
        Path projectionPath = required(directory, "projection.onnx");
        required(directory, "checksums.sha256");
        CompactSkelSpheresManifest manifest = new ObjectMapper().readValue(
                manifestPath.toFile(), CompactSkelSpheresManifest.class);
        manifest.validateSource(sourceManifest);
        verify(projectionPath, manifest.projectionSha256);
        String bundleHash = sha256(manifestPath) + ":" + manifest.projectionSha256;
        return new CompactSkelSpheresModelBundle(
                directory.toAbsolutePath().normalize(), manifest, projectionPath, bundleHash);
    }

    private static Path required(Path directory, String name) throws IOException {
        Path path = directory.resolve(name);
        if (!Files.isRegularFile(path)) throw new IOException("missing compact bundle file: " + name);
        return path;
    }

    private static void verify(Path path, String expected) throws IOException {
        if (!sha256(path).equals(expected)) {
            throw new IOException("SHA-256 mismatch for " + path.getFileName());
        }
    }

    static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream input = Files.newInputStream(path)) {
                byte[] buffer = new byte[1024 * 1024];
                for (int read; (read = input.read(buffer)) >= 0;) {
                    digest.update(buffer, 0, read);
                }
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException impossible) {
            throw new IllegalStateException(impossible);
        }
    }
}
