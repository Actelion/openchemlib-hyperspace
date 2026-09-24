package com.idorsia.research.chem.hyperspace.mcp;

import java.io.File;
import java.nio.channels.*;
import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.TimeUnit;

final class Worker {
    final Path dir;
    final ServerConfig config;
    final Map<String, Object> request;
    final Map<String, Object> status;
    final Map<String, Object> artifacts = new LinkedHashMap<>();

    Worker(Path dir) throws Exception {
        this.dir = dir;
        config = ServerConfig.load(dir.resolve("server.json"));
        request = JsonFiles.read(dir.resolve("request.json"));
        status = JsonFiles.read(dir.resolve("status.json"));
    }

    void save() throws Exception {
        status.put("heartbeatAt", Instant.now().toString());
        JsonFiles.write(dir.resolve("status.json"), status);
    }

    static void identity(Map<String, Object> s, String prefix, ProcessHandle process) {
        s.put(prefix + "Pid", process.pid());
        s.put(prefix + "StartedAt", process.info().startInstant().orElseThrow().toString());
    }

    void run() throws Exception {
        identity(status, "worker", ProcessHandle.current());
        status.put("state", "running");
        status.put("artifacts", artifacts);
        save();
        try (FileChannel channel =
                        FileChannel.open(
                                config.workspace().resolve("worker.lock"),
                                StandardOpenOption.CREATE,
                                StandardOpenOption.WRITE);
                FileLock lock = channel.tryLock()) {
            if (lock == null)
                throw new IllegalStateException("Another worker holds the workspace lock");
            Path staging = dir.resolve("staging");
            Files.createDirectories(staging);
            Path raw;
            if (request.get("format").equals("rawspace"))
                raw = Path.of((String) request.get("input"));
            else {
                raw = staging.resolve("space.rawspace.gz");
                stage(
                        "import",
                        Preparation.importCommand(
                                config, request, raw, dir.resolve("import-report.json")),
                        dir.resolve("import-report.json"));
                Path published = dir.resolve("space.rawspace.gz");
                Files.move(raw, published);
                raw = published;
                updateOutput("import", "rawspace", published);
            }
            artifacts.put("rawspace", raw.toString());
            save();
            Path index = staging.resolve("space_FragFp.data"),
                    similarity = staging.resolve("space_FragFp_similarity3.data");
            stage(
                    "build",
                    Preparation.buildCommand(
                            config,
                            request,
                            raw,
                            index,
                            similarity,
                            dir.resolve("build-report.json")),
                    dir.resolve("build-report.json"));
            Path published = dir.resolve(index.getFileName());
            Files.move(index, published);
            artifacts.put("substructure", published.toString());
            updateOutput("build", "substructure", published);
            if (((List<?>) request.get("searchModes")).contains("similarity")) {
                Path sim = dir.resolve(similarity.getFileName());
                Files.move(similarity, sim);
                artifacts.put("similarity", sim.toString());
                updateOutput("build", "similarity", sim);
            }
            if (Files.exists(dir.resolve("cancel"))) throw new InterruptedException("Cancelled");
            status.put("stage", "configure_gui");
            save();
            List<Object> providers = new ArrayList<>();
            Jobs.addProviders(
                    providers, request, artifacts, dir.getFileName().toString().substring(0, 8));
            Path gui = dir.resolve("gui.json");
            JsonFiles.write(gui, Map.of("ServiceProviders", providers));
            artifacts.put("guiConfig", gui.toString());
            if (Boolean.TRUE.equals(request.get("launchGui"))) {
                try {
                    status.put(
                            "gui",
                            Jobs.launch(
                                    config,
                                    gui,
                                    (String) request.get("guiHeap"),
                                    dir.resolve("gui.log"),
                                    Preparation.guiVersion(request.getOrDefault("guiVersion", "gui2"))));
                } catch (Exception e) {
                    status.put("gui", Map.of("state", "launch_failed", "error", e.toString()));
                }
            }
            status.put("state", "succeeded");
            status.put("stage", "complete");
            status.remove("progress");
        } catch (InterruptedException e) {
            status.put("state", "cancelled");
            status.put("error", e.getMessage());
        } catch (Exception e) {
            status.put("state", "failed");
            status.put("error", e.toString());
            e.printStackTrace();
        } finally {
            status.put("finishedAt", Instant.now().toString());
            save();
        }
    }

    @SuppressWarnings("unchecked")
    private void updateOutput(String stage, String key, Path published) throws Exception {
        Map<String, Object> summary = (Map<String, Object>) status.get(stage + "Summary");
        ((Map<String, Object>) summary.get("outputs")).put(key, published.toString());
        JsonFiles.write(dir.resolve(stage + "-report.json"), summary);
    }

    @SuppressWarnings("unchecked")
    void stage(String stage, List<String> command, Path report) throws Exception {
        if (Files.exists(dir.resolve("cancel")))
            throw new InterruptedException("Cancelled before " + stage);
        status.put("stage", stage);
        status.put("progress", "See stage log; no estimated percentage available");
        JsonFiles.write(dir.resolve(stage + "-command.json"), command);
        save();
        Process process =
                new ProcessBuilder(command)
                        .directory(dir.toFile())
                        .redirectInput(new File("/dev/null"))
                        .redirectErrorStream(true)
                        .redirectOutput(dir.resolve(stage + ".log").toFile())
                        .start();
        try {
            identity(status, "child", process.toHandle());
            save();
            while (!process.waitFor(1, TimeUnit.SECONDS)) {
                if (Files.exists(dir.resolve("cancel")))
                    throw new InterruptedException("Cancelled during " + stage);
                save();
            }
            status.put(stage + "ExitCode", process.exitValue());
            if (Files.exists(dir.resolve("cancel")))
                throw new InterruptedException("Cancelled during " + stage);
            if (process.exitValue() != 0)
                throw new IllegalStateException(
                        stage
                                + " failed with exit code "
                                + process.exitValue()
                                + "; see "
                                + stage
                                + ".log");
            Map<String, Object> summary = JsonFiles.read(report);
            if (!Boolean.TRUE.equals(summary.get("success")))
                throw new IllegalStateException(
                        "Missing successful completion report for " + stage);
            if (((Number) summary.get("reactionCount")).longValue() == 0)
                throw new IllegalStateException(
                        "No valid reactions were produced; see filtering diagnostics");
            for (Object path : ((Map<String, Object>) summary.get("outputs")).values()) {
                if (!Files.isRegularFile(Path.of(path.toString()))
                        || Files.size(Path.of(path.toString())) == 0)
                    throw new IllegalStateException("Missing/empty output: " + path);
            }
            status.put(stage + "Summary", summary);
        } finally {
            if (process.isAlive()) Jobs.terminate(process.toHandle());
            status.remove("childPid");
            status.remove("childStartedAt");
            save();
        }
    }
}
