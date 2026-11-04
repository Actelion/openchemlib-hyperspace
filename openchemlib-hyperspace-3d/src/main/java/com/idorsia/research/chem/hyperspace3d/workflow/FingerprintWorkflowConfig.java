package com.idorsia.research.chem.hyperspace3d.workflow;

import com.idorsia.research.chem.hyperspace3d.index.build.MoleculeFingerprintIndexBuildConfig;
import java.nio.file.Path;

/** Cluster paths are resolved on the machine executing the workflow, not by MCP. */
public final class FingerprintWorkflowConfig {
    public int formatVersion = 1;
    public String library;
    public String smilesColumn = "smiles";
    public String idColumn = "id";
    public String sharedDirectory;
    public String scratchRoot = "/scratch";
    public String applicationJar;
    public String libraryDirectory;
    public String modelBundle;
    public String compactBundle;
    public String environmentSetup;
    public String javaExecutable = "java";
    public String heap = "8G";
    public int rowsPerPartition = 1_000_000;
    public int maxArrayTasks = 1000;
    public MoleculeFingerprintIndexBuildConfig.Runtime runtime = new MoleculeFingerprintIndexBuildConfig.Runtime();

    public FingerprintWorkflowConfig() {
        runtime.cpuWorkers = 6;
        runtime.encoderBatchSize = 512;
    }

    public void resolve(Path configFile) {
        Path base = configFile.toAbsolutePath().getParent();
        library = resolve(base, library);
        sharedDirectory = resolve(base, sharedDirectory);
        scratchRoot = resolve(base, scratchRoot);
        applicationJar = resolve(base, applicationJar);
        libraryDirectory = resolve(base, libraryDirectory);
        modelBundle = resolve(base, modelBundle);
        compactBundle = resolve(base, compactBundle);
        if (environmentSetup != null) environmentSetup = resolve(base, environmentSetup);
        if (javaExecutable.contains("/")) javaExecutable = resolve(base, javaExecutable);
        if (formatVersion != 1 || rowsPerPartition < 1 || rowsPerPartition > 8_000_000
                || maxArrayTasks < 1 || maxArrayTasks > 100000
                || !heap.matches("[1-9][0-9]*[mMgG]") || javaExecutable.isBlank())
            throw new IllegalArgumentException("invalid workflow version, partition size, array limit or JVM settings");
        build(Path.of(library), Path.of(sharedDirectory).resolve("cache"),
                Path.of(modelBundle), Path.of(compactBundle)).validate();
    }

    public MoleculeFingerprintIndexBuildConfig build(Path input, Path output, Path base, Path compact) {
        var c = new MoleculeFingerprintIndexBuildConfig();
        c.inputs.library = input.toString();
        c.inputs.smilesColumn = smilesColumn;
        c.inputs.idColumn = idColumn;
        c.inputs.modelBundle = base.toString();
        c.inputs.compactBundle = compact.toString();
        c.output.directory = output.toString();
        c.output.sourceRowsPerShard = rowsPerPartition;
        c.runtime = runtime;
        return c;
    }

    private static String resolve(Path base, String text) {
        if (text == null || text.isBlank()) throw new IllegalArgumentException("required workflow path is missing");
        return base.resolve(text).normalize().toAbsolutePath().toString();
    }
}
