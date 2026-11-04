package com.idorsia.research.chem.hyperspace.mcp;

import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

final class JsonFiles {
    static final ObjectMapper JSON = new ObjectMapper();

    static Map<String, Object> read(Path path) throws IOException {
        return JSON.readValue(
                path.toFile(),
                new com.fasterxml.jackson.core.type.TypeReference<Map<String, Object>>() {});
    }

    static void write(Path path, Object value) throws IOException {
        Files.createDirectories(path.toAbsolutePath().getParent());
        Path tmp = Files.createTempFile(path.toAbsolutePath().getParent(), ".json-", ".tmp");
        try {
            JSON.writerWithDefaultPrettyPrinter().writeValue(tmp.toFile(), value);
            try {
                Files.move(
                        tmp,
                        path,
                        StandardCopyOption.ATOMIC_MOVE,
                        StandardCopyOption.REPLACE_EXISTING);
            } catch (AtomicMoveNotSupportedException e) {
                Files.move(tmp, path, StandardCopyOption.REPLACE_EXISTING);
            }
        } finally {
            Files.deleteIfExists(tmp);
        }
    }

    static String string(Map<String, Object> values, String key) {
        Object v = values.get(key);
        if (!(v instanceof String) || ((String) v).isBlank())
            throw new IllegalArgumentException(key + " must be a nonempty string");
        return (String) v;
    }

    static int integer(Map<String, Object> values, String key, int fallback, int min, int max) {
        Object v = values.getOrDefault(key, fallback);
        if (!(v instanceof Number) || ((Number) v).doubleValue() != ((Number) v).intValue())
            throw new IllegalArgumentException(key + " must be an integer");
        int n = ((Number) v).intValue();
        if (n < min || n > max)
            throw new IllegalArgumentException(key + " must be between " + min + " and " + max);
        return n;
    }
}
