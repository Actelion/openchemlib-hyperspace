package com.idorsia.research.chem.hyperspace.mcp;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;

class WorkerTest {
    static Worker worker(Path root) throws Exception {
        String id = UUID.randomUUID().toString();
        Path d = root.resolve("jobs").resolve(id);
        JsonFiles.write(
                d.resolve("server.json"),
                Map.of("workspace", root.toString(), "cliJar", root.resolve("cli.jar").toString()));
        JsonFiles.write(d.resolve("request.json"), Map.of());
        JsonFiles.write(
                d.resolve("status.json"),
                Map.of("jobId", id, "state", "pending", "createdAt", Instant.now().toString()));
        return new Worker(d);
    }

    static List<String> fixture(String mode) {
        return List.of(
                Path.of(System.getProperty("java.home"), "bin", "java").toString(),
                "-cp",
                System.getProperty("java.class.path"),
                Fixture.class.getName(),
                mode);
    }

    public static class Fixture {
        public static void main(String[] args) throws Exception {
            if (args[0].equals("fail")) {
                System.err.println("intentional failure");
                System.exit(7);
            }
            if (args[0].equals("sleep")) {
                System.out.println("running");
                Thread.sleep(60000);
            }
        }
    }

    @Test
    void cancellationStopsTheComputation() throws Exception {
        Worker w = worker(Files.createTempDirectory("worker-cancel"));
        ExecutorService executor = Executors.newSingleThreadExecutor();
        Future<?> task =
                executor.submit(
                        () -> {
                            assertThrows(
                                    InterruptedException.class,
                                    () ->
                                            w.stage(
                                                    "import",
                                                    fixture("sleep"),
                                                    w.dir.resolve("report.json")));
                        });
        try {
            long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(10);
            Map<String, Object> state = Map.of();
            while (System.nanoTime() < deadline) {
                state = JsonFiles.read(w.dir.resolve("status.json"));
                if (state.containsKey("childPid")) break;
                Thread.sleep(20);
            }
            assertTrue(state.containsKey("childPid"));
            long pid = ((Number) state.get("childPid")).longValue();
            Files.writeString(w.dir.resolve("cancel"), "cancel");
            task.get(15, TimeUnit.SECONDS);
            assertFalse(ProcessHandle.of(pid).map(ProcessHandle::isAlive).orElse(false));
        } finally {
            executor.shutdownNow();
        }
    }

    @Test
    void failureAndMissingReportsCannotPublishSuccess() throws Exception {
        Worker w = worker(Files.createTempDirectory("worker-fail"));
        Exception e =
                assertThrows(
                        IllegalStateException.class,
                        () -> w.stage("import", fixture("fail"), w.dir.resolve("report.json")));
        assertTrue(e.getMessage().contains("7"));
        assertTrue(w.artifacts.isEmpty());
        assertThrows(
                java.io.IOException.class,
                () -> w.stage("build", fixture("success"), w.dir.resolve("missing.json")));
        assertTrue(w.artifacts.isEmpty());
    }

    @Test
    void duplicateSubmissionReturnsActiveJob() throws Exception {
        Path root = Files.createTempDirectory("worker-busy");
        Worker w = worker(root);
        Worker.identity(w.status, "worker", ProcessHandle.current());
        w.status.put("state", "running");
        w.save();
        Path jar = Files.writeString(root.resolve("cli.jar"), "unused");
        Path input = Files.writeString(root.resolve("input.txt"), "unused");
        Jobs jobs = new Jobs(new ServerConfig(root, jar, "java", 1, "1G", "1G"));
        Map<String, Object> result =
                jobs.submit(
                        Map.of(
                                "input",
                                input.toString(),
                                "format",
                                "enamine",
                                "spaceName",
                                "Test"));
        assertEquals(false, result.get("accepted"));
        assertEquals("workspace_busy", result.get("reason"));
    }
}
