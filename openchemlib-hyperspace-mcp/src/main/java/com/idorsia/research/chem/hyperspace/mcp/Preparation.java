package com.idorsia.research.chem.hyperspace.mcp;

import java.nio.file.*;
import java.util.*;

final class Preparation {
    private static final Set<String> FORMATS = Set.of("enamine", "xtalpi", "csv", "rawspace");
    static final Set<String> IMPORT_OPTIONS =
            Set.of(
                    "zipEntry",
                    "reactionZipEntry",
                    "sourceFormat",
                    "maxSets",
                    "smilesColumn",
                    "idColumn",
                    "reactionColumn",
                    "priceColumn",
                    "priceAttribute",
                    "synthonSetColumn",
                    "synthonDefaultSet",
                    "metadata");

    static Map<String, Object> normalize(Map<String, Object> input, ServerConfig config) {
        Set<String> allowed =
                Set.of(
                        "input",
                        "format",
                        "spaceName",
                        "searchModes",
                        "importOptions",
                        "threads",
                        "heap",
                        "guiHeap",
                        "guiVersion",
                        "launchGui");
        for (String k : input.keySet())
            if (!allowed.contains(k)) throw new IllegalArgumentException("Unknown parameter: " + k);
        Map<String, Object> request = new LinkedHashMap<>(input);
        Path source = Path.of(JsonFiles.string(input, "input")).toAbsolutePath().normalize();
        String format = JsonFiles.string(input, "format");
        if (!FORMATS.contains(format))
            throw new IllegalArgumentException("format must be enamine, xtalpi, csv, or rawspace");
        if (!Files.isReadable(source)
                || (format.equals("csv")
                        ? !Files.isDirectory(source)
                        : !Files.isRegularFile(source)))
            throw new IllegalArgumentException(
                    "Input does not exist or has the wrong file/directory type: " + source);
        JsonFiles.string(input, "spaceName");
        Object modes = input.getOrDefault("searchModes", List.of("substructure", "similarity"));
        if (!(modes instanceof List<?>)
                || ((List<?>) modes).isEmpty()
                || !Set.of("substructure", "similarity").containsAll((List<?>) modes))
            throw new IllegalArgumentException(
                    "searchModes must contain substructure and/or similarity");
        request.put("searchModes", new ArrayList<>(new LinkedHashSet<>((List<?>) modes)));
        Object options = input.getOrDefault("importOptions", Map.of());
        if (!(options instanceof Map<?, ?>))
            throw new IllegalArgumentException("importOptions must be an object");
        Map<?, ?> opts = (Map<?, ?>) options;
        Set<String> valid =
                switch (format) {
                    case "enamine" ->
                            Set.of(
                                    "zipEntry",
                                    "reactionZipEntry",
                                    "sourceFormat",
                                    "maxSets",
                                    "metadata");
                    case "xtalpi" ->
                            Set.of(
                                    "sourceFormat",
                                    "maxSets",
                                    "smilesColumn",
                                    "idColumn",
                                    "reactionColumn",
                                    "synthonSetColumn",
                                    "metadata");
                    case "csv" ->
                            Set.of(
                                    "smilesColumn",
                                    "idColumn",
                                    "priceColumn",
                                    "priceAttribute",
                                    "synthonSetColumn",
                                    "synthonDefaultSet",
                                    "metadata");
                    default -> Set.of();
                };
        for (var e : opts.entrySet()) {
            if (!valid.contains(e.getKey()))
                throw new IllegalArgumentException(
                        "Unsupported option for " + format + ": " + e.getKey());
            if (e.getKey().equals("metadata")) {
                if (!(e.getValue() instanceof Map<?, ?>))
                    throw new IllegalArgumentException("metadata must be an object of strings");
                for (var a : ((Map<?, ?>) e.getValue()).entrySet()) {
                    if (!(a.getKey() instanceof String)
                            || !(a.getValue() instanceof String)
                            || a.getKey().toString().contains("="))
                        throw new IllegalArgumentException(
                                "metadata keys and values must be strings; keys cannot contain"
                                    + " '='");
                }
            } else if (e.getKey().equals("maxSets") || e.getKey().equals("synthonDefaultSet")) {
                Object n = e.getValue();
                if (!(n instanceof Number)
                        || ((Number) n).doubleValue() != ((Number) n).intValue()
                        || ((Number) n).intValue() < 0)
                    throw new IllegalArgumentException(
                            e.getKey() + " must be a nonnegative integer");
            } else if (!(e.getValue() instanceof String) || e.getValue().toString().isBlank())
                throw new IllegalArgumentException(e.getKey() + " must be a nonempty string");
        }
        if (source.toString().endsWith(".zip") && !opts.containsKey("zipEntry"))
            throw new IllegalArgumentException("Select a zipEntry from inspect_input first");
        if (opts.containsKey("reactionZipEntry") && !opts.containsKey("zipEntry"))
            throw new IllegalArgumentException("reactionZipEntry requires zipEntry");
        if (input.containsKey("launchGui") && !(input.get("launchGui") instanceof Boolean))
            throw new IllegalArgumentException("launchGui must be boolean");
        request.put("input", source.toString());
        request.put("importOptions", opts);
        request.put("threads", JsonFiles.integer(input, "threads", config.threads(), 1, 4096));
        request.put(
                "heap", ServerConfig.validHeap((String) input.getOrDefault("heap", config.heap())));
        request.put(
                "guiHeap",
                ServerConfig.validHeap((String) input.getOrDefault("guiHeap", config.guiHeap())));
        request.put("launchGui", input.getOrDefault("launchGui", false));
        request.put("guiVersion", guiVersion(input.getOrDefault("guiVersion", "gui2")));
        return request;
    }

    static String guiVersion(Object value) {
        if (!Set.of("gui2", "legacy").contains(value))
            throw new IllegalArgumentException("guiVersion must be gui2 or legacy");
        return (String) value;
    }

    static List<String> javaCommand(ServerConfig c, String heap, String main) {
        return new ArrayList<>(
                List.of(c.javaExecutable(), "-Xmx" + heap, "-cp", c.cliJar().toString(), main));
    }

    @SuppressWarnings("unchecked")
    static List<String> importCommand(
            ServerConfig c, Map<String, Object> r, Path raw, Path report) {
        List<String> a =
                javaCommand(
                        c,
                        (String) r.get("heap"),
                        "com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceImportCLI");
        String format = (String) r.get("format");
        a.addAll(
                List.of(
                        "--format",
                        format,
                        format.equals("csv") ? "--directory" : "--input",
                        (String) r.get("input"),
                        "--spaceName",
                        (String) r.get("spaceName"),
                        "--mode",
                        "FragFp",
                        "--threads",
                        r.get("threads").toString(),
                        "--rawOut",
                        raw.toString(),
                        "--reportOut",
                        report.toString()));
        Map<String, Object> opts = (Map<String, Object>) r.get("importOptions");
        if (!format.equals("csv") && !opts.containsKey("maxSets"))
            a.addAll(List.of("--maxSets", "3"));
        for (var e : opts.entrySet()) {
            if (e.getKey().equals("metadata"))
                ((Map<String, String>) e.getValue())
                        .forEach((k, v) -> a.addAll(List.of("--metadata", k + "=" + v)));
            else a.addAll(List.of("--" + e.getKey(), e.getValue().toString()));
        }
        a.addAll(List.of("--metadata", "space.role=full"));
        return a;
    }

    static List<String> buildCommand(
            ServerConfig c,
            Map<String, Object> r,
            Path raw,
            Path index,
            Path similarity,
            Path report) {
        List<String> a =
                javaCommand(
                        c,
                        (String) r.get("heap"),
                        "com.idorsia.research.chem.hyperspace.cli.RawSynthonSpaceBuildCLI");
        a.addAll(
                List.of(
                        "--rawIn",
                        raw.toString(),
                        "--threads",
                        r.get("threads").toString(),
                        "--synthonOut",
                        index.toString(),
                        "--descriptor",
                        "FragFp",
                        "--bits",
                        "512",
                        "--reportOut",
                        report.toString()));
        if (((List<?>) r.get("searchModes")).contains("similarity"))
            a.addAll(
                    List.of(
                            "--similarityOut",
                            similarity.toString(),
                            "--similarityThreads",
                            r.get("threads").toString()));
        return a;
    }
}
