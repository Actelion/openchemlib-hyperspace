package com.idorsia.research.chem.hyperspace.gui2;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Startup configuration shared with the MCP and legacy GUI. */
public final class Gui2Configuration {
    public record Space(Path file, String name, int threads) {}

    public record Configuration(List<Space> spaces, List<String> warnings) {}

    private Gui2Configuration() {}

    public static Configuration read(Path input) throws IOException {
        Path config = input.toAbsolutePath().normalize();
        List<Space> spaces = new ArrayList<>();
        List<String> warnings = new ArrayList<>();
        if (config.getFileName().toString().endsWith(".json")) {
            JsonNode providers =
                    new ObjectMapper().readTree(config.toFile()).path("ServiceProviders");
            if (!providers.isArray()) throw new IOException("Missing ServiceProviders array");
            for (JsonNode provider : providers) {
                if (!"HyperspaceSSS".equals(provider.path("ServiceProvider").asText())) {
                    warnings.add(
                            "Skipped "
                                    + provider.path("ServiceName").asText("provider")
                                    + ": GUI 2 supports local substructure search only; use the"
                                    + " legacy GUI for similarity.");
                    continue;
                }
                JsonNode data = provider.path("Config");
                String file = data.path("File").asText("");
                if (file.isBlank()) throw new IOException("Missing index File");
                JsonNode limit = data.get("MaxNumberOfThreads");
                if (limit != null
                        && (!limit.isIntegralNumber()
                                || !limit.canConvertToInt()
                                || limit.intValue() < 1))
                    throw new IOException("MaxNumberOfThreads must be positive");
                int threads =
                        limit == null
                                ? Runtime.getRuntime().availableProcessors()
                                : limit.intValue();
                String name =
                        provider.path("ServiceName").asText(data.path("SpaceName").asText(file));
                spaces.add(new Space(resolve(config, file), name, threads));
            }
        } else {
            for (String line : Files.readAllLines(config)) {
                if (line.isBlank() || line.stripLeading().startsWith("#")) continue;
                Path file = resolve(config, line.trim());
                spaces.add(
                        new Space(
                                file,
                                file.getFileName().toString(),
                                Runtime.getRuntime().availableProcessors()));
            }
        }
        if (spaces.isEmpty())
            throw new IOException(
                    "No local substructure spaces configured. Use guiVersion: legacy for similarity"
                        + " search.");
        return new Configuration(List.copyOf(spaces), List.copyOf(warnings));
    }

    private static Path resolve(Path config, String value) throws IOException {
        Path file = config.getParent().resolve(value).normalize();
        if (!Files.isRegularFile(file) || !Files.isReadable(file))
            throw new IOException("Index is not readable: " + file);
        return file;
    }
}
