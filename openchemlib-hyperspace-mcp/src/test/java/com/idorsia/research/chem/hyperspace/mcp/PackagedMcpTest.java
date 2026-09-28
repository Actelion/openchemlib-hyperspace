package com.idorsia.research.chem.hyperspace.mcp;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import io.modelcontextprotocol.client.*;
import io.modelcontextprotocol.client.transport.*;
import io.modelcontextprotocol.json.McpJsonDefaults;
import io.modelcontextprotocol.spec.McpSchema;

import org.junit.jupiter.api.Test;

import java.nio.file.*;
import java.time.Duration;
import java.util.*;

class PackagedMcpTest {
    @Test
    void fingerprintPreparationSurvivesClientReconnect() throws Exception {
        String jar = System.getProperty("hyperspace3d.distributionJar");
        assumeTrue(jar != null && System.getProperty("hyperspace.mcpJar") != null);
        Path root = Files.createTempDirectory("mcp fingerprint ");
        Path input = Files.writeString(root.resolve("input.tsv"), "smiles\tid\nCCCCCC\tone\nc1ccccc1\ttwo\n");
        Path config = root.resolve("server.json");
        JsonFiles.write(config, Map.of("workspace", root.resolve("workspace").toString(),
                "cliJar", System.getProperty("hyperspace.cliJar", jar), "fingerprintJar", jar,
                "fingerprintLibDirectory", Path.of(jar).getParent().resolve("lib").toString(),
                "modelBundle", System.getProperty("hyperspace3d.modelBundle"),
                "compactBundle", System.getProperty("hyperspace3d.compactBundle"), "threads", 2, "heap", "1G"));
        String id;
        try (var client = connect(config)) {
            assertTrue(call(client, "get_workflow_help", Map.of("topic", "fingerprints")).get("text").toString().contains("scratch"));
            assertTrue(call(client, "inspect_molecule_library", Map.of("input", input.toString())).containsKey("columns"));
            var started = call(client, "prepare_molecule_fingerprints", Map.of("input", input.toString(), "sourceRowsPerShard", 1));
            assertEquals(true, started.get("accepted")); id = started.get("jobId").toString();
        }
        try (var client = connect(config)) {
            Map<String, Object> status;
            long end = System.nanoTime() + Duration.ofMinutes(2).toNanos();
            do {
                status = call(client, "get_job_status", Map.of("jobId", id));
                if (!Jobs.ACTIVE.contains(status.get("state"))) break;
                Thread.sleep(100);
            } while (System.nanoTime() < end);
            assertEquals("succeeded", status.get("state"), status.toString());
            assertEquals(1, ((List<?>) call(client, "list_molecule_libraries", Map.of()).get("libraries")).size());
            assertTrue(((List<?>) call(client, "list_spaces", Map.of()).get("spaces")).isEmpty());
            assertTrue(client.callTool(new McpSchema.CallToolRequest("configure_gui", Map.of("jobIds", List.of(id)))).isError());
            var artifacts = (Map<?, ?>) status.get("artifacts");
            assertTrue(Files.isRegularFile(Path.of(artifacts.get("manifest").toString())));
            assertTrue(call(client, "get_job_log", Map.of("jobId", id, "log", "fingerprint")).get("text").toString().contains("Molecule index complete"));
        }
    }

    McpSyncClient connect(Path config) {
        ServerParameters parameters =
                ServerParameters.builder(
                                Path.of(System.getProperty("java.home"), "bin", "java").toString())
                        .args(
                                "-jar",
                                System.getProperty("hyperspace.mcpJar"),
                                "--config",
                                config.toString())
                        .addEnvVar("DISPLAY", "")
                        .build();
        StdioClientTransport transport =
                new StdioClientTransport(parameters, McpJsonDefaults.getMapper());
        transport.setStdErrorHandler(System.err::println);
        McpSyncClient client =
                McpClient.sync(transport).requestTimeout(Duration.ofSeconds(30)).build();
        client.initialize();
        return client;
    }

    @SuppressWarnings("unchecked")
    Map<String, Object> call(McpSyncClient c, String name, Map<String, Object> args) {
        var result = c.callTool(new McpSchema.CallToolRequest(name, args));
        assertFalse(Boolean.TRUE.equals(result.isError()), result.toString());
        return (Map<String, Object>) result.structuredContent();
    }

    @Test
    void toyPreparationSurvivesClientReconnect() throws Exception {
        assumeTrue(
                System.getProperty("hyperspace.cliJar") != null
                        && System.getProperty("hyperspace.mcpJar") != null);
        Path root = Files.createTempDirectory("mcp packaged ");
        Path config = root.resolve("server.json");
        JsonFiles.write(
                config,
                Map.of(
                        "workspace",
                        root.resolve("workspace").toString(),
                        "cliJar",
                        System.getProperty("hyperspace.cliJar"),
                        "threads",
                        2,
                        "heap",
                        "2G",
                        "guiHeap",
                        "2G"));
        String id;
        try (McpSyncClient client = connect(config)) {
            assertEquals(15, client.listTools().tools().size());
            assertEquals(4, client.listResources().resources().size());
            assertEquals(1, client.listPrompts().prompts().size());
            assertTrue(
                    call(
                                    client,
                                    "inspect_input",
                                    Map.of("input", System.getProperty("hyperspace.toy")))
                            .containsKey("table"));
            assertTrue(
                    call(client, "get_workflow_help", Map.of("topic", "workflow"))
                            .get("text")
                            .toString()
                            .contains("prepare_space"));
            assertFalse(
                    client.readResource(client.listResources().resources().get(0))
                            .contents()
                            .isEmpty());
            assertTrue(
                    client.callTool(
                                    new McpSchema.CallToolRequest(
                                            "get_job_status", Map.of("jobId", "invalid")))
                            .isError());
            Map<String, Object> started =
                    call(
                            client,
                            "prepare_space",
                            Map.of(
                                    "input",
                                    System.getProperty("hyperspace.toy"),
                                    "format",
                                    "enamine",
                                    "spaceName",
                                    "Toy reconnect",
                                    "launchGui",
                                    false));
            assertEquals(true, started.get("accepted"));
            id = started.get("jobId").toString();
            assertTrue(
                    Jobs.ACTIVE.contains(
                            call(client, "get_job_status", Map.of("jobId", id)).get("state")));
        }
        try (McpSyncClient client = connect(config)) {
            Map<String, Object> status = Map.of();
            long deadline = System.nanoTime() + Duration.ofMinutes(3).toNanos();
            do {
                status = call(client, "get_job_status", Map.of("jobId", id));
                if (!Jobs.ACTIVE.contains(status.get("state"))) break;
                Thread.sleep(500);
            } while (System.nanoTime() < deadline);
            assertEquals("succeeded", status.get("state"), status.toString());
            Map<?, ?> artifacts = (Map<?, ?>) status.get("artifacts");
            assertTrue(Files.size(Path.of(artifacts.get("substructure").toString())) > 0);
            assertTrue(Files.size(Path.of(artifacts.get("similarity").toString())) > 0);
            assertFalse(((List<?>) call(client, "list_spaces", Map.of()).get("spaces")).isEmpty());
            Path gui =
                    Path.of(
                            call(client, "configure_gui", Map.of("jobIds", List.of(id)))
                                    .get("configPath")
                                    .toString());
            assertEquals(2, ((List<?>) JsonFiles.read(gui).get("ServiceProviders")).size());
            Map<?, ?> imported = (Map<?, ?>) status.get("importSummary");
            assertEquals(5, ((Number) imported.get("reactionCount")).intValue());
            assertFalse(((List<?>) imported.get("warnings")).isEmpty());
            String cancelId =
                    call(
                                    client,
                                    "prepare_space",
                                    Map.of(
                                            "input",
                                            artifacts.get("rawspace"),
                                            "format",
                                            "rawspace",
                                            "spaceName",
                                            "Cancel"))
                            .get("jobId")
                            .toString();
            call(client, "cancel_job", Map.of("jobId", cancelId));
            long cancelDeadline = System.nanoTime() + Duration.ofSeconds(30).toNanos();
            Map<String, Object> cancelled;
            do {
                cancelled = call(client, "get_job_status", Map.of("jobId", cancelId));
                if (!Jobs.ACTIVE.contains(cancelled.get("state"))) break;
                Thread.sleep(100);
            } while (System.nanoTime() < cancelDeadline);
            assertEquals("cancelled", cancelled.get("state"));
            Map<String,Object> launch = call(client, "launch_gui", Map.of("configPath", gui.toString()));
            assertEquals("not_launched", launch.get("state"));
            assertTrue(((List<?>)launch.get("command")).contains(gui.toString()));
            assertEquals("gui2", launch.get("guiVersion"));
            assertTrue(((List<?>) launch.get("command")).contains("GUI2"));
            assertFalse(((List<?>) launch.get("warnings")).isEmpty());
            Map<String,Object> legacy = call(client, "launch_gui",
                    Map.of("configPath", gui.toString(), "guiVersion", "legacy"));
            assertTrue(((List<?>) legacy.get("command")).contains("GUI"));
            System.out.println("MCP_TOY_GUI=" + gui);
            String save = System.getProperty("hyperspace.testManifest");
            if (save != null) JsonFiles.write(Path.of(save), status);
        }
    }
}
