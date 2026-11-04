package com.idorsia.research.chem.hyperspace.mcp;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

/** Generates files only. No cluster connection or shell evaluation occurs during export. */
final class SlurmFingerprintBundle {
    private static final String CLI = "com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintWorkflowCLI";

    @SuppressWarnings("unchecked")
    static Map<String, Object> export(Path profile, Path output) throws Exception {
        Map<String, Object> document = JsonFiles.read(profile);
        if (!(document.get("workflow") instanceof Map<?, ?>) || !(document.get("slurm") instanceof Map<?, ?>))
            throw new IllegalArgumentException("profile needs workflow and slurm objects");
        Map<String, Object> w = new LinkedHashMap<>((Map<String, Object>) document.get("workflow"));
        Map<String, Object> s = new LinkedHashMap<>((Map<String, Object>) document.get("slurm"));
        allowed(w, Set.of("formatVersion", "library", "smilesColumn", "idColumn", "sharedDirectory", "scratchRoot",
                "applicationJar", "libraryDirectory", "modelBundle", "compactBundle", "javaExecutable", "heap",
                "rowsPerPartition", "maxArrayTasks", "runtime"));
        allowed(s, Set.of("partition", "cpuPartition", "account", "qos", "constraint", "gpuRequest", "cpus",
                "memory", "time", "cpuMemory", "cpuTime", "maxConcurrentTasks", "environmentSetup"));
        w.putIfAbsent("formatVersion", 1);
        if (JsonFiles.integer(w, "formatVersion", 1, 1, 1) != 1) throw new IllegalArgumentException("unsupported format");
        w.putIfAbsent("scratchRoot", "/scratch");
        for (String key : List.of("library", "sharedDirectory", "scratchRoot", "applicationJar", "libraryDirectory", "modelBundle", "compactBundle"))
            absolute(w, key);
        w.putIfAbsent("javaExecutable", "java");
        safeText(JsonFiles.string(w, "javaExecutable"));
        w.put("heap", ServerConfig.validHeap((String) w.getOrDefault("heap", "8G")));
        w.put("rowsPerPartition", JsonFiles.integer(w, "rowsPerPartition", 1_000_000, 1, 8_000_000));
        w.put("maxArrayTasks", JsonFiles.integer(w, "maxArrayTasks", 1000, 1, 100000));
        w.putIfAbsent("smilesColumn", "smiles"); w.putIfAbsent("idColumn", "id");
        for (String key : List.of("smilesColumn", "idColumn")) {
            String value = JsonFiles.string(w, key); safeText(value);
            if (value.contains("\t")) throw new IllegalArgumentException("column names cannot contain tabs");
        }
        if (w.get("smilesColumn").equals(w.get("idColumn"))) throw new IllegalArgumentException("columns must differ");
        int cpus = JsonFiles.integer(s, "cpus", 8, 1, 4096);
        var runtime = new LinkedHashMap<>((Map<String, Object>) w.getOrDefault("runtime", Map.of()));
        allowed(runtime, Set.of("device", "cudaDeviceId", "cpuWorkers", "encoderBatchSize", "queueCapacity", "progressIntervalSeconds"));
        if (!Objects.toString(runtime.getOrDefault("device", "CUDA")).equalsIgnoreCase("CUDA"))
            throw new IllegalArgumentException("Slurm GPU bundle requires runtime.device CUDA");
        runtime.put("device", "CUDA");
        runtime.put("cudaDeviceId", JsonFiles.integer(runtime, "cudaDeviceId", 0, 0, 0));
        runtime.put("cpuWorkers", JsonFiles.integer(runtime, "cpuWorkers", Math.max(1, cpus - 2), 1, cpus));
        runtime.put("encoderBatchSize", JsonFiles.integer(runtime, "encoderBatchSize", 512, 1, 65536));
        runtime.put("queueCapacity", JsonFiles.integer(runtime, "queueCapacity", 4, 1, 1024));
        runtime.put("progressIntervalSeconds", JsonFiles.integer(runtime, "progressIntervalSeconds", 30, 0, 3600));
        w.put("runtime", runtime);
        String partition = directive(s, "partition", null);
        String cpuPartition = directive(s, "cpuPartition", partition);
        String gpu = directive(s, "gpuRequest", "gpu:1");
        if (!gpu.matches("gpu(?::[A-Za-z0-9_]+)?:1")) throw new IllegalArgumentException("request exactly one GPU per task");
        String memory = memory(s, "memory", "24G");
        String cpuMemory = memory(s, "cpuMemory", "8G");
        String time = time(s, "time", "24:00:00"), cpuTime = time(s, "cpuTime", "24:00:00");
        int concurrent = JsonFiles.integer(s, "maxConcurrentTasks", 4, 1, 100000);
        String setup = "";
        if (s.containsKey("environmentSetup")) {
            absolute(s, "environmentSetup");
            w.put("environmentSetup", s.get("environmentSetup"));
            setup = "source " + quote((String) s.get("environmentSetup")) + "\n";
        }
        String common = "";
        for (String k : List.of("account", "qos")) if (s.containsKey(k)) common += "#SBATCH --" + k + "=" + directive(s, k, null) + "\n";
        String constraint = s.containsKey("constraint") ? "#SBATCH --constraint=" + directive(s, "constraint", null) + "\n" : "";
        String cp = w.get("applicationJar") + ":" + w.get("libraryDirectory") + "/*";
        String java = quote((String) w.get("javaExecutable"));
        String launch = java + " -Xmx4G -cp " + quote(cp) + " " + CLI;
        Map<String, String> scripts = new LinkedHashMap<>();
        for (String phase : List.of("prepare", "smoke", "compute", "finalize")) {
            boolean gpuPhase = phase.equals("smoke") || phase.equals("compute");
            String header = "#!/usr/bin/env bash\n#SBATCH --job-name=hs-fp-" + phase + "\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n"
                    + common + "#SBATCH --partition=" + (gpuPhase ? partition : cpuPartition) + "\n"
                    + "#SBATCH --cpus-per-task=" + (gpuPhase ? cpus : 2) + "\n"
                    + "#SBATCH --mem=" + (gpuPhase ? memory : cpuMemory) + "\n"
                    + "#SBATCH --time=" + (gpuPhase ? time : cpuTime) + "\n"
                    + (gpuPhase ? "#SBATCH --gres=" + gpu + "\n" + constraint : "");
            String operation = phase.equals("compute") ? "worker" : phase;
            String args = phase.equals("compute") ? " --task \"${SLURM_ARRAY_TASK_ID:?}\" --tasks \"${2:?array task count}\"" : "";
            scripts.put(phase + ".sbatch", header + "set -euo pipefail\nCONFIG=${1:?absolute workflow config required}\n"
                    + "cd -- \"$(dirname -- \"$CONFIG\")\"\n" + setup
                    + "srun " + launch + " " + operation + " --config \"$CONFIG\"" + args + "\n");
            String helper = "#!/usr/bin/env bash\nset -euo pipefail\n"
                    + "BUNDLE=$(cd -- \"$(dirname -- \"${BASH_SOURCE[0]}\")\" && pwd)\n"
                    + "CONFIG=\"$BUNDLE/workflow-config.json\"\nmkdir -p -- \"$BUNDLE/logs\"\n" + setup;
            if (phase.equals("compute")) helper += "TASKS=$(" + launch + " array-tasks --config \"$CONFIG\")\n"
                    + "[[ $TASKS =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid array size' >&2; exit 1; }\n";
            helper += "sbatch --parsable --chdir=\"$BUNDLE\" --output=\"$BUNDLE/logs/%x-%A_%a.out\" "
                    + (phase.equals("compute") ? "--array=\"0-$((TASKS-1))%" + concurrent + "\" " : "")
                    + "\"$BUNDLE/" + phase + ".sbatch\" \"$CONFIG\"" + (phase.equals("compute") ? " \"$TASKS\"" : "") + "\n";
            scripts.put("submit_" + phase + ".sh", helper);
        }
        output = output.toAbsolutePath().normalize();
        if (Files.exists(output)) try (var entries = Files.list(output)) {
            if (entries.findAny().isPresent()) throw new IOException("bundle output must be empty");
        }
        Files.createDirectories(output);
        JsonFiles.write(output.resolve("workflow-config.json"), w);
        JsonFiles.write(output.resolve("cluster-profile.json"), s);
        for (var script : scripts.entrySet()) Files.writeString(output.resolve(script.getKey()), script.getValue(), StandardOpenOption.CREATE_NEW);
        Files.writeString(output.resolve("README.md"), "# Fingerprint job bundle\n\n"
                + "Copy this bundle to cluster shared storage. Verify workflow-config.json, application/model paths and environment setup.\n\n"
                + "Run `bash submit_prepare.sh`; wait for success. Then run `bash submit_smoke.sh`; wait for success.\n"
                + "Run `bash submit_compute.sh`; wait for every array task. Finally run `bash submit_finalize.sh`.\n\n"
                + "The helpers submit jobs only when YOU execute them. Exporting this bundle submits nothing.\n"
                + "Use squeue/sacct to monitor the printed IDs; logs are in logs/.\n"
                + "Re-submit compute to skip verified completed partitions. Partial scratch work is recomputed.\n"
                + "Re-submit finalize to reuse verified cache shards. Do not edit inputs, models or configuration after preparation.\n"
                + "Final cache: " + w.get("sharedDirectory") + "/cache\n\n"
                + "Shared filesystem must support atomic rename and file locks. Scratch is task-private and cleaned on normal exit.\n"
                + "Hard-linked published payloads are immutable. Copy fallback requires additional storage. No automatic input cleanup.\n"
                + "CUDA/cuDNN availability and GPU support are cluster acceptance checks. No driver installation is performed.\n");
        return Map.of("bundleDirectory", output.toString(), "submitted", false,
                "scripts", scripts.keySet(), "workflowConfig", output.resolve("workflow-config.json").toString());
    }

    private static void allowed(Map<String, Object> m, Set<String> keys) {
        for (String key : m.keySet()) if (!keys.contains(key)) throw new IllegalArgumentException("Unknown profile key: " + key);
    }
    private static void safeText(String s) {
        if (s.chars().anyMatch(c -> c == 0 || c == '\n' || c == '\r')) throw new IllegalArgumentException("control characters are not allowed");
    }
    private static void absolute(Map<String, Object> m, String key) {
        String value = JsonFiles.string(m, key); safeText(value);
        if (!Path.of(value).isAbsolute()) throw new IllegalArgumentException(key + " must be an absolute cluster path");
    }
    private static String directive(Map<String, Object> s, String key, String fallback) {
        String value = (String) s.getOrDefault(key, fallback);
        if (value == null || !value.matches("[A-Za-z0-9_.:-]+")) throw new IllegalArgumentException("invalid Slurm " + key);
        return value;
    }
    private static String memory(Map<String, Object> s, String key, String fallback) {
        return ServerConfig.validHeap((String) s.getOrDefault(key, fallback));
    }
    private static String time(Map<String, Object> s, String key, String fallback) {
        String value = (String) s.getOrDefault(key, fallback);
        if (!value.matches("[0-9]+(?:-[0-9]+)?:[0-5][0-9]:[0-5][0-9]")) throw new IllegalArgumentException("invalid time limit");
        return value;
    }
    static String quote(String text) { safeText(text); return "'" + text.replace("'", "'\"'\"'") + "'"; }
}
