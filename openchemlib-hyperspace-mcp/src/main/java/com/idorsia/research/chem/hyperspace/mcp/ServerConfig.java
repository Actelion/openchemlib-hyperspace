package com.idorsia.research.chem.hyperspace.mcp;

import java.nio.file.*;
import java.util.*;

record ServerConfig(
        Path workspace,
        Path cliJar,
        String javaExecutable,
        int threads,
        String heap,
        String guiHeap,
        Path fingerprintJar,
        Path fingerprintLibDirectory,
        Path modelBundle,
        Path compactBundle) {
    ServerConfig(Path workspace, Path cliJar, String javaExecutable, int threads, String heap, String guiHeap) {
        this(workspace, cliJar, javaExecutable, threads, heap, guiHeap, null, null, null, null);
    }
    static ServerConfig load(Path path) throws Exception {
        Map<String, Object> m = JsonFiles.read(path);
        Path base = path.toAbsolutePath().getParent();
        return new ServerConfig(
                base.resolve(JsonFiles.string(m, "workspace")).normalize(),
                base.resolve(JsonFiles.string(m, "cliJar")).normalize(),
                (String)
                        m.getOrDefault(
                                "javaExecutable",
                                Path.of(System.getProperty("java.home"), "bin", "java").toString()),
                JsonFiles.integer(
                        m,
                        "threads",
                        Math.min(4, Runtime.getRuntime().availableProcessors()),
                        1,
                        4096),
                validHeap((String) m.getOrDefault("heap", "8G")),
                validHeap((String) m.getOrDefault("guiHeap", "8G")),
                optionalPath(m, base, "fingerprintJar"),
                optionalPath(m, base, "fingerprintLibDirectory"),
                optionalPath(m, base, "modelBundle"),
                optionalPath(m, base, "compactBundle"));
    }

    private static Path optionalPath(Map<String, Object> m, Path base, String key) {
        return m.containsKey(key) ? base.resolve(JsonFiles.string(m, key)).toAbsolutePath().normalize() : null;
    }

    static String validHeap(String value) {
        if (value == null || !value.matches("[1-9][0-9]*[mMgG]"))
            throw new IllegalArgumentException(
                    "Heap must be a positive integer followed by M or G");
        return value;
    }

    Map<String, Object> asMap() {
        Map<String, Object> result = new LinkedHashMap<>(Map.of(
                "workspace",
                workspace.toString(),
                "cliJar",
                cliJar.toString(),
                "javaExecutable",
                javaExecutable,
                "threads",
                threads,
                "heap",
                heap,
                "guiHeap",
                guiHeap));
        if (fingerprintJar != null) result.put("fingerprintJar", fingerprintJar.toString());
        if (fingerprintLibDirectory != null) result.put("fingerprintLibDirectory", fingerprintLibDirectory.toString());
        if (modelBundle != null) result.put("modelBundle", modelBundle.toString());
        if (compactBundle != null) result.put("compactBundle", compactBundle.toString());
        return result;
    }

    static String setsid() {
        for (String p : List.of("/usr/bin/setsid", "/bin/setsid"))
            if (Files.isExecutable(Path.of(p))) return p;
        throw new IllegalStateException("Linux setsid is required for independent jobs");
    }
}
