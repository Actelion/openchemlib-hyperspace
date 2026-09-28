package com.idorsia.research.chem.hyperspace.mcp;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.nio.file.*;
import java.util.*;
import static org.junit.jupiter.api.Assertions.*;

class FingerprintToolsTest {
    @TempDir Path temp;

    @Test void legacyConfigAndBoundedInspection() throws Exception {
        Path config = temp.resolve("server.json");
        JsonFiles.write(config, Map.of("workspace", temp.toString(), "cliJar", "cli.jar"));
        var old = ServerConfig.load(config);
        assertNull(old.fingerprintJar());
        assertThrows(IllegalArgumentException.class, () -> FingerprintJobs.normalize(Map.of(), old));
        Path input = temp.resolve("library.tsv");
        Files.writeString(input, "smiles\tid\n" + "CCO\tx\n".repeat(30000));
        var info = FingerprintJobs.inspect(input);
        assertEquals(List.of("smiles", "id"), info.get("columns"));
        assertEquals(5, ((List<?>) info.get("sampleLines")).size());
        Path models = Files.createDirectories(temp.resolve("models"));
        var current = new ServerConfig(temp, temp.resolve("cli.jar"), "java", 2, "2G", "2G",
                temp.resolve("fp.jar"), models, models, models);
        var request = FingerprintJobs.normalize(Map.of("input", input.toString()), current);
        assertEquals("CPU", request.get("device"));
        assertEquals("fingerprint", request.get("jobKind"));
        assertTrue(FingerprintJobs.command(current, request, config).contains("-cp"));
        assertThrows(IllegalArgumentException.class, () -> FingerprintJobs.normalize(
                Map.of("input", input.toString(), "threads", 0), current));
    }

    Map<String, Object> profile() {
        Map<String, Object> w = new LinkedHashMap<>();
        for (String key : List.of("library", "sharedDirectory", "applicationJar", "libraryDirectory", "modelBundle", "compactBundle"))
            w.put(key, "/shared/with space/'quoted/" + key);
        return new LinkedHashMap<>(Map.of("workflow", w,
                "slurm", new LinkedHashMap<>(Map.of("partition", "GPU-01", "cpuPartition", "CPU-01", "gpuRequest", "gpu:rtx5080:1"))));
    }

    @SuppressWarnings("unchecked")
    @Test void exporterQuotesPathsAndHelpersOnlySubmitWhenInvoked() throws Exception {
        var p = profile();
        Path bin = Files.createDirectories(temp.resolve("bin"));
        Path java = bin.resolve("java-fixture");
        Files.writeString(java, "#!/bin/bash\n[[ ${FAIL_SMOKE:-0} == 0 ]] || exit 1\nprintf '3\\n'\n");
        java.toFile().setExecutable(true);
        ((Map<String, Object>) p.get("workflow")).put("javaExecutable", java.toString());
        Path submit = bin.resolve("sbatch");
        Files.writeString(submit, "#!/bin/bash\nprintf '%s\\n' \"$@\" > \"$SUBMIT_LOG\"\nprintf '12345\\n'\n");
        submit.toFile().setExecutable(true);
        Path profile = temp.resolve("profile.json"), bundle = temp.resolve("export with spaces");
        JsonFiles.write(profile, p);
        var exported = SlurmFingerprintBundle.export(profile, bundle);
        assertEquals(false, exported.get("submitted"));
        assertFalse(Files.exists(temp.resolve("submitted")));
        try (var paths = Files.list(bundle)) {
            for (Path script : paths.filter(f -> f.toString().endsWith(".sh") || f.toString().endsWith(".sbatch")).toList()) {
                Process check = new ProcessBuilder("bash", "-n", script.toString()).inheritIO().start();
                assertEquals(0, check.waitFor(), script.toString());
            }
        }
        var command = new ProcessBuilder("bash", bundle.resolve("submit_compute.sh").toString());
        command.environment().put("PATH", bin + ":" + System.getenv("PATH"));
        command.environment().put("SUBMIT_LOG", temp.resolve("submitted").toString());
        assertEquals(0, command.inheritIO().start().waitFor());
        assertTrue(Files.readString(temp.resolve("submitted")).contains("--array=0-2%4"));
        Files.delete(temp.resolve("submitted"));
        command.environment().put("FAIL_SMOKE", "1");
        assertNotEquals(0, command.inheritIO().start().waitFor());
        assertFalse(Files.exists(temp.resolve("submitted")));
        assertTrue(Files.readString(bundle.resolve("compute.sbatch")).contains("gpu:rtx5080:1"));
        assertThrows(java.io.IOException.class, () -> SlurmFingerprintBundle.export(profile, bundle));
    }

    @SuppressWarnings("unchecked")
    @Test void rejectsUnsafeSlurmDirectives() throws Exception {
        var p = profile();
        ((Map<String, Object>) p.get("slurm")).put("partition", "GPU-01\n#SBATCH --exclusive");
        Path profile = temp.resolve("profile.json"); JsonFiles.write(profile, p);
        assertThrows(IllegalArgumentException.class, () -> SlurmFingerprintBundle.export(profile, temp.resolve("out")));
        assertFalse(Files.exists(temp.resolve("out")));
    }

    @Test void fingerprintCancellationUsesExistingProcessLifecycle() throws Exception {
        Worker w = WorkerTest.worker(temp);
        Files.writeString(w.dir.resolve("cancel"), "cancel");
        assertThrows(InterruptedException.class, () -> w.stage("fingerprint", WorkerTest.fixture("sleep"), temp.resolve("manifest.json")));
        assertFalse(w.status.containsKey("childPid"));
    }
}
