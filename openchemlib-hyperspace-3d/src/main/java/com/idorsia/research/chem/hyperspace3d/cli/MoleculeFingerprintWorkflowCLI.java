package com.idorsia.research.chem.hyperspace3d.cli;

import com.idorsia.research.chem.hyperspace3d.workflow.*;
import java.nio.file.Path;

/** Shared-filesystem orchestration; GPU computation and application data live on task-local scratch. */
public final class MoleculeFingerprintWorkflowCLI {
    private MoleculeFingerprintWorkflowCLI() {}
    public static void main(String[] args) throws Exception {
        if (args.length < 3 || !args[1].equals("--config"))
            throw new IllegalArgumentException("Usage: MoleculeFingerprintWorkflowCLI prepare|smoke|worker|finalize|array-tasks --config workflow.json [--task N --tasks N]");
        var workflow = new FingerprintWorkflow(Path.of(args[2]));
        if (!args[0].equals("worker") && args.length != 3) throw new IllegalArgumentException("unexpected arguments");
        switch (args[0]) {
            case "prepare" -> workflow.prepare();
            case "smoke" -> workflow.smoke();
            case "finalize" -> System.out.println(FingerprintFinalizer.finish(workflow));
            case "array-tasks" -> System.out.println(workflow.arrayTasks());
            case "worker" -> {
                if (args.length != 7 || !args[3].equals("--task") || !args[5].equals("--tasks"))
                    throw new IllegalArgumentException("worker requires --task N --tasks N");
                workflow.worker(Integer.parseInt(args[4]), Integer.parseInt(args[6]));
            }
            default -> throw new IllegalArgumentException("unknown workflow command");
        }
    }
}
