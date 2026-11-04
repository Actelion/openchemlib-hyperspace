package com.idorsia.research.chem.hyperspace.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.nio.file.*;
import java.util.*;

final class RawSpaceCompletionReport {
    static void write(String output, RawSynthonSpace raw, Map<String, String> outputs)
            throws Exception {
        if (output == null) return;
        Map<String, Object> report = new LinkedHashMap<>();
        report.put("success", true);
        report.put("spaceName", raw.getName());
        report.put("reactionCount", raw.getReactions().size());
        long count =
                raw.getReactions().values().stream()
                        .flatMap(r -> r.getRawFragmentSets().values().stream())
                        .mapToLong(List::size)
                        .sum();
        report.put("synthonCount", count);
        Map<String, String> files = new LinkedHashMap<>();
        outputs.forEach((k, v) -> files.put(k, Path.of(v).toAbsolutePath().toString()));
        report.put("outputs", files);
        Map<String, String> diagnostics = new LinkedHashMap<>();
        raw.getMetadata()
                .forEach(
                        (k, v) -> {
                            if (k.startsWith("parser.")) diagnostics.put(k, v);
                        });
        report.put("importDiagnostics", diagnostics);
        List<String> warnings = new ArrayList<>();
        if (diagnostics.containsKey("parser.inputReactionCount")) {
            int input = Integer.parseInt(diagnostics.get("parser.inputReactionCount"));
            if (input > raw.getReactions().size())
                warnings.add(
                        (input - raw.getReactions().size())
                                + " input reactions were excluded by import filtering/validation;"
                                + " see import log.");
        }
        report.put("warnings", warnings);
        Path target = Path.of(output).toAbsolutePath();
        Files.createDirectories(target.getParent());
        Path tmp = Files.createTempFile(target.getParent(), "completion-", ".json");
        try {
            new ObjectMapper().writerWithDefaultPrettyPrinter().writeValue(tmp.toFile(), report);
            try {
                Files.move(
                        tmp,
                        target,
                        StandardCopyOption.ATOMIC_MOVE,
                        StandardCopyOption.REPLACE_EXISTING);
            } catch (AtomicMoveNotSupportedException e) {
                Files.move(tmp, target, StandardCopyOption.REPLACE_EXISTING);
            }
        } finally {
            Files.deleteIfExists(tmp);
        }
    }
}
