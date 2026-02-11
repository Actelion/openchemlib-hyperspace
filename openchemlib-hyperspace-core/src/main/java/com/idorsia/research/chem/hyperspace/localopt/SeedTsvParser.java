package com.idorsia.research.chem.hyperspace.localopt;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SeedTsvParser {

    public static List<SeedAssembly> parse(Path path) throws IOException {
        List<SeedAssembly> seeds = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            if (header == null) {
                return seeds;
            }
            String[] columns = header.split("	");
            Map<String, Integer> index = new HashMap<>();
            for (int i = 0; i < columns.length; i++) {
                index.put(columns[i], i);
            }
            Integer rxnIdx = index.get("rxnId");
            Integer fragIdx = index.get("fragIds");
            Integer scoreIdx = index.get("phesaSimilarity");
            if (rxnIdx == null || fragIdx == null) {
                throw new IllegalArgumentException("Seed file must contain rxnId and fragIds columns");
            }
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                String[] parts = line.split("	");
                String reactionId = parts[rxnIdx];
                List<String> fragIds = parseFragments(parts[fragIdx]);
                double score = Double.NaN;
                if (scoreIdx != null && scoreIdx < parts.length) {
                    try {
                        score = Double.parseDouble(parts[scoreIdx]);
                    } catch (NumberFormatException ignored) {
                    }
                }
                seeds.add(new SeedAssembly(reactionId, fragIds, score));
            }
        }
        return seeds;
    }

    private static List<String> parseFragments(String column) {
        List<String> fragments = new ArrayList<>();
        if (column == null || column.isEmpty()) {
            return fragments;
        }
        String[] splits = column.split(";\\s*");
        for (String s : splits) {
            String trimmed = s.trim();
            if (!trimmed.isEmpty()) {
                fragments.add(trimmed);
            }
        }
        return fragments;
    }
}
