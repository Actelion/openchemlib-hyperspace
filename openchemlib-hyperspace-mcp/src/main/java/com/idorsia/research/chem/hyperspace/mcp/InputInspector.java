package com.idorsia.research.chem.hyperspace.mcp;

import com.fasterxml.jackson.core.*;

import org.apache.commons.csv.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.*;

final class InputInspector {
    static final int LIMIT = 1024 * 1024;

    static Map<String, Object> inspect(Path input, String entry) throws Exception {
        Path p = input.toAbsolutePath().normalize();
        if (!Files.isReadable(p)) throw new IllegalArgumentException("Input is not readable: " + p);
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("input", p.toString());
        result.put("sampleOnly", true);
        result.put("note", "Samples identify layout, not chemical validity or whole-file counts.");
        if (Files.isDirectory(p)) {
            List<Object> samples = new ArrayList<>();
            try (var files = Files.list(p)) {
                for (Path f :
                        files.filter(Files::isRegularFile)
                                .filter(f -> f.toString().endsWith(".csv"))
                                .limit(5)
                                .toList())
                    try (InputStream in = Files.newInputStream(f)) {
                        samples.add(Map.of("file", f.getFileName().toString(), "table", table(in)));
                    }
            }
            result.put("samples", samples);
            result.put("formatCandidates", List.of("csv"));
            result.put(
                    "note",
                    "Per-reaction CSV: confirm reaction/set filename conventions and column"
                        + " mappings using RAWSPACE_IMPORTERS.");
        } else if (p.toString().endsWith(".zip")) {
            try (ZipFile zip = new ZipFile(p.toFile())) {
                List<String> entries =
                        zip.stream()
                                .filter(e -> !e.isDirectory())
                                .map(ZipEntry::getName)
                                .limit(100)
                                .toList();
                result.put("entries", entries);
                result.put("entriesTruncated", zip.size() > 100);
                if (entry != null) {
                    ZipEntry e = zip.getEntry(entry);
                    if (e == null || e.isDirectory())
                        throw new IllegalArgumentException("ZIP entry not found: " + entry);
                    try (InputStream in = zip.getInputStream(e)) {
                        result.put("table", table(in));
                    }
                    result.put("formatCandidates", List.of("enamine"));
                } else
                    result.put(
                            "required",
                            "Choose zipEntry; optionally reactionZipEntry. No archive contents are"
                                + " extracted.");
            }
        } else if (p.toString().endsWith(".rawspace") || p.toString().endsWith(".rawspace.gz")) {
            result.put("formatCandidates", List.of("rawspace"));
            try (InputStream file = Files.newInputStream(p);
                    InputStream in =
                            p.toString().endsWith(".gz") ? new GZIPInputStream(file) : file) {
                byte[] prefix = in.readNBytes(LIMIT);
                Map<String, Object> info = new LinkedHashMap<>();
                try (JsonParser parser = JsonFiles.JSON.getFactory().createParser(prefix)) {
                    if (parser.nextToken() != JsonToken.START_OBJECT)
                        throw new IllegalArgumentException(
                                "Rawspace must start with a JSON object");
                    while (parser.nextToken() == JsonToken.FIELD_NAME) {
                        String name = parser.currentName();
                        parser.nextToken();
                        if (name.equals("reactions")) break;
                        if (Set.of("name", "version", "metadata").contains(name))
                            info.put(name, JsonFiles.JSON.readValue(parser, Object.class));
                        else parser.skipChildren();
                    }
                } catch (com.fasterxml.jackson.core.io.JsonEOFException e) {
                    info.put("sampleTruncated", true);
                }
                result.put("rawspace", info);
            }
        } else
            try (InputStream in = Files.newInputStream(p)) {
                Map<String, Object> t = table(in);
                result.put("table", t);
                List<?> header = (List<?>) t.get("header");
                List<String> lower =
                        header.stream().map(v -> v.toString().toLowerCase(Locale.ROOT)).toList();
                List<String> formats = new ArrayList<>();
                if (lower.contains("synton_role")) formats.add("xtalpi");
                if (lower.contains("reaction_id")
                        && (lower.contains("synthon#") || lower.contains("synton#")))
                    formats.add("enamine");
                result.put("formatCandidates", formats);
                result.put("requiresExplicitFormat", true);
            }
        return result;
    }

    static Map<String, Object> table(InputStream in) throws Exception {
        byte[] bytes = in.readNBytes(LIMIT);
        String text = new String(bytes, StandardCharsets.UTF_8);
        // Drop an incomplete last record at the sampling boundary; quoted multiline fields may
        // still be truncated.
        if (bytes.length == LIMIT) {
            int end = text.lastIndexOf('\n');
            text = end < 0 ? "" : text.substring(0, end + 1);
        }
        String first = text.lines().findFirst().orElse("");
        char delimiter = first.contains("\t") ? '\t' : ',';
        List<List<String>> records = new ArrayList<>();
        String warning = "";
        try (CSVParser parser =
                CSVFormat.DEFAULT
                        .builder()
                        .setDelimiter(delimiter)
                        .setIgnoreEmptyLines(true)
                        .get()
                        .parse(new StringReader(text))) {
            Iterator<CSVRecord> it = parser.iterator();
            while (records.size() < 21 && it.hasNext()) {
                List<String> row = new ArrayList<>();
                it.next().forEach(row::add);
                records.add(row);
            }
        } catch (UncheckedIOException e) {
            warning =
                    "Sample ends inside a quoted record or contains malformed CSV: "
                            + e.getMessage();
        }
        return Map.of(
                "delimiter",
                delimiter == '\t' ? "tab" : "comma",
                "header",
                records.isEmpty() ? List.of() : records.get(0),
                "records",
                records.size() < 2 ? List.of() : records.subList(1, records.size()),
                "byteLimit",
                LIMIT,
                "limitReached",
                bytes.length == LIMIT,
                "warning",
                warning);
    }
}
