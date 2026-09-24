package com.idorsia.research.chem.hyperspace.mcp;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

import java.nio.file.*;
import java.util.*;

class PreparationTest {
    @Test
    void guiVersionsAreExplicit() {
        assertEquals("gui2", Preparation.guiVersion("gui2"));
        assertEquals("legacy", Preparation.guiVersion("legacy"));
        assertThrows(IllegalArgumentException.class, () -> Preparation.guiVersion("other"));
    }

    @Test
    void commandsPreservePathsAndRejectUnsupportedOptions() throws Exception {
        Path root = Files.createTempDirectory("mcp spaces ");
        Path input =
                Files.writeString(root.resolve("input with spaces.txt"), "SMILES\tid\tset\trx\n");
        ServerConfig c = new ServerConfig(root, root.resolve("cli.jar"), "java", 2, "2G", "2G");
        Map<String, Object> r =
                Preparation.normalize(
                        Map.of(
                                "input",
                                input.toString(),
                                "format",
                                "enamine",
                                "spaceName",
                                "Test; $(not a shell)"),
                        c);
        assertEquals("gui2", r.get("guiVersion"));
        Map<String,Object> legacyInput = new java.util.HashMap<>(r);
        legacyInput.put("guiVersion", "legacy");
        legacyInput.put("launchGui", true);
        assertEquals("legacy", Preparation.normalize(legacyInput, c).get("guiVersion"));
        List<String> cmd =
                Preparation.importCommand(
                        c, r, root.resolve("raw.gz"), root.resolve("report.json"));
        assertTrue(cmd.contains(input.toString()));
        assertTrue(cmd.contains("Test; $(not a shell)"));
        assertEquals("3", cmd.get(cmd.indexOf("--maxSets") + 1));
        assertThrows(
                IllegalArgumentException.class,
                () ->
                        Preparation.normalize(
                                Map.of(
                                        "input",
                                        input.toString(),
                                        "format",
                                        "enamine",
                                        "spaceName",
                                        "Test",
                                        "importOptions",
                                        Map.of("idColumn", "id")),
                                c));
        assertThrows(IllegalArgumentException.class, () -> ServerConfig.validHeap("8G -evil"));
        assertThrows(
                IllegalArgumentException.class,
                () ->
                        Preparation.normalize(
                                Map.of(
                                        "input",
                                        input.toString(),
                                        "format",
                                        "enamine",
                                        "spaceName",
                                        "Test",
                                        "threads",
                                        0),
                                c));
    }

    @Test
    void mapsExistingImporterOptionsAndRawspace() throws Exception {
        Path root = Files.createTempDirectory("mcp formats");
        Path input =
                Files.writeString(
                        root.resolve("input.csv"), "SMILES,synton_id,synton_role,reaction_id\n");
        ServerConfig c = new ServerConfig(root, root.resolve("cli.jar"), "java", 2, "2G", "2G");
        Map<String, Object> x =
                Preparation.normalize(
                        Map.of(
                                "input",
                                input.toString(),
                                "format",
                                "xtalpi",
                                "spaceName",
                                "X",
                                "importOptions",
                                Map.of("idColumn", "supplier", "maxSets", 0)),
                        c);
        assertTrue(
                Preparation.importCommand(c, x, root.resolve("raw"), root.resolve("report"))
                        .contains("supplier"));
        Map<String, Object> raw =
                Preparation.normalize(
                        Map.of(
                                "input",
                                input.toString(),
                                "format",
                                "rawspace",
                                "spaceName",
                                "X",
                                "searchModes",
                                List.of("substructure")),
                        c);
        List<String> build =
                Preparation.buildCommand(
                        c,
                        raw,
                        input,
                        root.resolve("idx"),
                        root.resolve("sim"),
                        root.resolve("report"));
        assertFalse(build.contains("--similarityOut"));
        assertEquals(raw.get("threads").toString(), build.get(build.indexOf("--threads") + 1));
        assertEquals("512", build.get(build.indexOf("--bits") + 1));
    }

    @Test
    void inspectionIsBoundedAndHandlesQuotedCsv() throws Exception {
        Path f = Files.createTempFile("mcp-inspect-", ".csv");
        String text = "SMILES,synthon_id,synthon#,reaction_id\n";
        for (int i = 0; i < 40; i++) text += "C,\"id," + i + "\",0,rx\n";
        Files.writeString(f, text);
        Map<?, ?> table = (Map<?, ?>) InputInspector.inspect(f, null).get("table");
        assertEquals(20, ((List<?>) table.get("records")).size());
        assertEquals("id,0", ((List<?>) ((List<?>) table.get("records")).get(0)).get(1));
    }

    @Test
    void zipInspectionRequiresExplicitEntryAndReadsIt() throws Exception {
        Path f = Files.createTempFile("mcp-zip-", ".zip");
        try (var zip = new java.util.zip.ZipOutputStream(Files.newOutputStream(f))) {
            zip.putNextEntry(new java.util.zip.ZipEntry("nested/synthons.txt"));
            zip.write("SMILES\tsynthon_id\tsynthon#\treaction_id\nC\tx\t0\tr\n".getBytes());
            zip.closeEntry();
        }
        assertTrue(InputInspector.inspect(f, null).containsKey("required"));
        assertTrue(InputInspector.inspect(f, "nested/synthons.txt").containsKey("table"));
        assertThrows(IllegalArgumentException.class, () -> InputInspector.inspect(f, "missing"));
    }

    @Test
    void staleWorkerIsNotReportedAsSuccess() throws Exception {
        Path root = Files.createTempDirectory("mcp-status");
        Jobs jobs =
                new Jobs(new ServerConfig(root, root.resolve("cli.jar"), "java", 2, "2G", "2G"));
        String id = UUID.randomUUID().toString();
        Path d = root.resolve("jobs").resolve(id);
        JsonFiles.write(
                d.resolve("status.json"),
                Map.of(
                        "jobId",
                        id,
                        "state",
                        "running",
                        "createdAt",
                        "2020-01-01T00:00:00Z",
                        "workerPid",
                        Long.MAX_VALUE,
                        "workerStartedAt",
                        "2020-01-01T00:00:00Z"));
        assertEquals("interrupted", jobs.status(id).get("state"));
        assertThrows(IllegalArgumentException.class, () -> jobs.dir("../outside"));
    }

    @Test
    void completedGuiConfigurationKeepsBothProviders() throws Exception {
        List<Object> providers = new ArrayList<>();
        Jobs.addProviders(
                providers,
                Map.of(
                        "spaceName",
                        "Toy",
                        "threads",
                        2,
                        "searchModes",
                        List.of("substructure", "similarity")),
                Map.of("substructure", "/tmp/sss.data", "similarity", "/tmp/sim.data"),
                "123");
        assertEquals(2, providers.size());
        assertEquals("HyperspaceSimilarity", ((Map<?, ?>) providers.get(1)).get("ServiceProvider"));
    }
}
