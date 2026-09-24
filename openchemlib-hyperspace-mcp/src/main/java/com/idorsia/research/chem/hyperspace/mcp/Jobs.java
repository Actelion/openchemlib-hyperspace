package com.idorsia.research.chem.hyperspace.mcp;

import java.io.*;
import java.nio.channels.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.TimeUnit;

final class Jobs {
    static final Set<String> ACTIVE = Set.of("pending", "running");
    final ServerConfig config;

    Jobs(ServerConfig config) throws IOException {
        this.config = config;
        Files.createDirectories(config.workspace().resolve("jobs"));
    }

    Path dir(String id) {
        if (!id.matches("[a-f0-9-]{36}")) throw new IllegalArgumentException("Invalid job ID");
        Path p = config.workspace().resolve("jobs").resolve(id);
        if (!Files.isDirectory(p)) throw new IllegalArgumentException("Unknown job: " + id);
        return p;
    }

    static String self() throws Exception {
        return Path.of(
                        HyperspaceMcp.class
                                .getProtectionDomain()
                                .getCodeSource()
                                .getLocation()
                                .toURI())
                .toString();
    }

    synchronized Map<String, Object> submit(Map<String, Object> input) throws Exception {
        Map<String, Object> r = Preparation.normalize(input, config);
        if (!Files.isRegularFile(config.cliJar()))
            throw new IllegalArgumentException("CLI JAR missing: " + config.cliJar());
        String setsid = ServerConfig.setsid();
        try (FileChannel channel =
                        FileChannel.open(
                                config.workspace().resolve("submit.lock"),
                                StandardOpenOption.CREATE,
                                StandardOpenOption.WRITE);
                FileLock lock = channel.lock()) {
            for (Map<String, Object> s : list())
                if (ACTIVE.contains(s.get("state")))
                    return Map.of("accepted", false, "reason", "workspace_busy", "activeJob", s);
            String id = UUID.randomUUID().toString();
            Path d = config.workspace().resolve("jobs").resolve(id);
            Files.createDirectories(d);
            JsonFiles.write(d.resolve("request.json"), r);
            JsonFiles.write(d.resolve("server.json"), config.asMap());
            Map<String, Object> state = new LinkedHashMap<>();
            state.put("jobId", id);
            state.put("state", "pending");
            state.put("stage", "starting");
            state.put("createdAt", Instant.now().toString());
            state.put("artifacts", new LinkedHashMap<>());
            state.put("cliJar", config.cliJar().toString());
            state.put("cliJarSha256", sha256(config.cliJar()));
            state.put("mcpVersion", "0.2.5");
            state.put("javaVersion", System.getProperty("java.version"));
            try (var jar = new java.util.jar.JarFile(config.cliJar().toFile())) {
                var manifest = jar.getManifest();
                state.put("cliJarVersion", manifest == null ? "unknown" : Objects.toString(
                        manifest.getMainAttributes().getValue("Implementation-Version"), "unknown"));
            }
            JsonFiles.write(d.resolve("status.json"), state);
            try {
                new ProcessBuilder(
                                setsid,
                                config.javaExecutable(),
                                "-Xmx256M",
                                "-cp",
                                self(),
                                HyperspaceMcp.class.getName(),
                                "--worker",
                                d.toString())
                        .redirectInput(new File("/dev/null"))
                        .redirectErrorStream(true)
                        .redirectOutput(d.resolve("worker.log").toFile())
                        .start();
            } catch (Exception e) {
                state.put("state", "failed");
                state.put("error", e.toString());
                JsonFiles.write(d.resolve("status.json"), state);
                throw e;
            }
            return Map.of("accepted", true, "jobId", id, "status", state);
        }
    }

    static String sha256(Path p) throws Exception {
        MessageDigest hash = MessageDigest.getInstance("SHA-256");
        try (InputStream in = Files.newInputStream(p)) {
            byte[] b = new byte[65536];
            int n;
            while ((n = in.read(b)) != -1) hash.update(b, 0, n);
        }
        return HexFormat.of().formatHex(hash.digest());
    }

    Map<String, Object> status(String id) throws Exception {
        Map<String, Object> s = JsonFiles.read(dir(id).resolve("status.json"));
        if (ACTIVE.contains(s.get("state"))) {
            boolean alive = s.containsKey("workerPid") && matchingProcess(s, "worker").isPresent();
            boolean old =
                    Instant.parse((String) s.get("createdAt"))
                            .isBefore(Instant.now().minusSeconds(60));
            if (!alive && (s.containsKey("workerPid") || old)) {
                s.put("state", "interrupted");
                s.put(
                        "error",
                        "Worker is no longer running; completed artifacts can be reused in a new"
                            + " job.");
                if (matchingProcess(s, "child").isPresent()) {
                    s.put("state", "running");
                    s.put(
                            "error",
                            "Worker stopped but recorded computation is still running; cancel this"
                                + " job before resubmitting.");
                }
            }
        }
        s.put(
                "elapsedSeconds",
                java.time.Duration.between(
                                Instant.parse((String) s.get("createdAt")),
                                Instant.parse(
                                        (String)
                                                s.getOrDefault(
                                                        "finishedAt", Instant.now().toString())))
                        .getSeconds());
        return s;
    }

    static Optional<ProcessHandle> matchingProcess(Map<String, Object> s, String prefix) {
        if (!(s.get(prefix + "Pid") instanceof Number)) return Optional.empty();
        return ProcessHandle.of(((Number) s.get(prefix + "Pid")).longValue())
                .filter(
                        p ->
                                p.isAlive()
                                        && p.info()
                                                .startInstant()
                                                .map(
                                                        t ->
                                                                t.toString()
                                                                        .equals(
                                                                                s.get(
                                                                                        prefix
                                                                                                + "StartedAt")))
                                                .orElse(false));
    }

    List<Map<String, Object>> list() throws Exception {
        List<Map<String, Object>> all = new ArrayList<>();
        try (var paths = Files.list(config.workspace().resolve("jobs"))) {
            for (Path p : paths.sorted().toList())
                if (Files.isRegularFile(p.resolve("status.json")))
                    all.add(status(p.getFileName().toString()));
        }
        return all;
    }

    Map<String, Object> cancel(String id) throws Exception {
        Path d = dir(id);
        Map<String, Object> s = status(id);
        if (ACTIVE.contains(s.get("state"))) {
            Files.writeString(d.resolve("cancel"), "cancel requested");
            if (matchingProcess(s, "worker").isEmpty())
                matchingProcess(s, "child").ifPresent(Jobs::terminate);
        }
        return Map.of("jobId", id, "cancellationRequested", ACTIVE.contains(s.get("state")));
    }

    static void terminate(ProcessHandle p) {
        List<ProcessHandle> children = p.descendants().toList();
        children.forEach(ProcessHandle::destroy);
        p.destroy();
        try {
            p.onExit().get(5, TimeUnit.SECONDS);
        } catch (Exception ignored) {
            p.destroyForcibly();
        }
        children.stream().filter(ProcessHandle::isAlive).forEach(ProcessHandle::destroyForcibly);
    }

    Map<String, Object> log(String id, String log, long offset, int limit) throws Exception {
        if (!Set.of("worker", "import", "build", "gui").contains(log))
            throw new IllegalArgumentException("log must be worker, import, build, or gui");
        if (offset < 0 || limit < 1 || limit > 65536)
            throw new IllegalArgumentException("Invalid offset or limit (maximum 65536 bytes)");
        Path f = dir(id).resolve(log + ".log");
        if (!Files.exists(f)) return Map.of("text", "", "nextOffset", 0, "eof", true);
        try (RandomAccessFile in = new RandomAccessFile(f.toFile(), "r")) {
            in.seek(Math.min(offset, in.length()));
            byte[] b = new byte[limit];
            int n = in.read(b);
            return Map.of(
                    "text",
                    new String(b, 0, Math.max(0, n), StandardCharsets.UTF_8),
                    "nextOffset",
                    in.getFilePointer(),
                    "eof",
                    in.getFilePointer() >= in.length());
        }
    }

    static boolean displayAvailable() {
        return System.getenv("DISPLAY") != null && !System.getenv("DISPLAY").isBlank();
    }

    static Map<String, Object> launch(ServerConfig c, Path gui, String heap, Path log)
            throws Exception {
        return launch(c, gui, heap, log, "gui2");
    }

    static Map<String, Object> launch(ServerConfig c, Path gui, String heap, Path log, String version)
            throws Exception {
        version = Preparation.guiVersion(version);
        List<String> warnings = new ArrayList<>();
        if (version.equals("gui2")) {
            Object value = JsonFiles.read(gui).get("ServiceProviders");
            if (!(value instanceof List<?>)) throw new IllegalArgumentException("Missing ServiceProviders");
            boolean substructure = false;
            for (Object item : (List<?>) value) {
                if (item instanceof Map<?, ?> provider
                        && "HyperspaceSSS".equals(provider.get("ServiceProvider"))) substructure = true;
                else warnings.add("GUI 2 supports local substructure search only. Use guiVersion: legacy for similarity.");
            }
            if (!substructure) throw new IllegalArgumentException(
                    "No local substructure providers. Use guiVersion: legacy for similarity search.");
        }
        List<String> cmd =
                List.of(
                        c.javaExecutable(),
                        "-Xmx" + ServerConfig.validHeap(heap),
                        "-jar",
                        c.cliJar().toString(),
                        version.equals("gui2") ? "GUI2" : "GUI",
                        gui.toString());
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("config", gui.toString());
        result.put("guiVersion", version);
        result.put("warnings", warnings.stream().distinct().toList());
        result.put("command", cmd);
        result.put(
                "shellCommand",
                cmd.stream()
                        .map(v -> "'" + v.replace("'", "'\"'\"'") + "'")
                        .collect(java.util.stream.Collectors.joining(" ")));
        if (!displayAvailable()) {
            result.put("state", "not_launched");
            result.put("reason", "No DISPLAY available; run the command in a desktop session.");
            return result;
        }
        List<String> detached = new ArrayList<>();
        detached.add(ServerConfig.setsid());
        detached.addAll(cmd);
        Process p =
                new ProcessBuilder(detached)
                        .redirectInput(new File("/dev/null"))
                        .redirectErrorStream(true)
                        .redirectOutput(ProcessBuilder.Redirect.appendTo(log.toFile()))
                        .start();
        result.put("pid", p.pid());
        result.put("log", log.toString());
        result.put("state", p.waitFor(500, TimeUnit.MILLISECONDS) ? "launch_failed" : "launched");
        result.put(
                "note",
                "Process launch does not mean indexes are ready. Check the GUI loading status.");
        return result;
    }

    @SuppressWarnings("unchecked")
    Path configure(List<String> ids, Path out) throws Exception {
        if (ids.isEmpty()) throw new IllegalArgumentException("Select at least one space/job ID");
        List<Object> providers = new ArrayList<>();
        for (String id : ids) {
            Map<String, Object> s = status(id);
            if (!s.get("state").equals("succeeded"))
                throw new IllegalArgumentException("Space is not complete: " + id);
            Map<String, Object> r = JsonFiles.read(dir(id).resolve("request.json"));
            Map<String, Object> a = (Map<String, Object>) s.get("artifacts");
            addProviders(providers, r, a, id.substring(0, 8));
        }
        JsonFiles.write(out, Map.of("ServiceProviders", providers));
        return out;
    }

    static void addProviders(
            List<Object> providers,
            Map<String, Object> r,
            Map<String, Object> artifacts,
            String suffix) {
        for (Object mode : (List<?>) r.get("searchModes")) {
            String key = mode.equals("similarity") ? "similarity" : "substructure";
            String name = r.get("spaceName") + " - " + key + " [" + suffix + "]";
            providers.add(
                    Map.of(
                            "ServiceName",
                            name,
                            "ServiceProvider",
                            key.equals("similarity") ? "HyperspaceSimilarity" : "HyperspaceSSS",
                            "Config",
                            Map.of(
                                    "SpaceName",
                                    r.get("spaceName"),
                                    "File",
                                    artifacts.get(key),
                                    "MaxNumberOfThreads",
                                    r.get("threads"))));
        }
    }
}
