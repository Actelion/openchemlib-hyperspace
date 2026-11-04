package com.idorsia.research.chem.hyperspace3d;

import com.idorsia.research.chem.hyperspace3d.workflow.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.cli.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;
import static com.idorsia.research.chem.hyperspace3d.workflow.WorkflowFiles.*;

class FingerprintWorkflowOnnxTest {
    @TempDir Path temp;
    @Test void packagedCpuWorkersMatchSingleBuildAndScreen() throws Exception {
        String jar = System.getProperty("hyperspace3d.distributionJar");
        assumeTrue(jar != null, "package the 3D distribution and set hyperspace3d.distributionJar");
        Path input = Files.writeString(temp.resolve("input.tsv"), FingerprintWorkflowTest.INPUT.replace("CCO", "CCCCCC").replace("CCN", "CCCCCN"));
        var c = new FingerprintWorkflowConfig();
        c.library = input.toString(); c.sharedDirectory = temp.resolve("shared").toString();
        c.scratchRoot = temp.resolve("scratch").toString();
        c.applicationJar = jar; c.libraryDirectory = Path.of(jar).getParent().resolve("lib").toString();
        c.modelBundle = Path.of(System.getProperty("hyperspace3d.modelBundle", "model-bundles/deepspace7-v1")).toAbsolutePath().toString();
        c.compactBundle = Path.of(System.getProperty("hyperspace3d.compactBundle", "model-bundles/deepspace7-skelspheres16")).toAbsolutePath().toString();
        c.rowsPerPartition = 2; c.runtime.device = "CPU"; c.runtime.cpuWorkers = 2;
        c.runtime.encoderBatchSize = 2; c.heap = "1G";
        c.javaExecutable = Path.of(System.getProperty("java.home"), "bin", "java").toString();
        Path config = temp.resolve("workflow.json"); json(config, c);
        var workflow = new FingerprintWorkflow(config);
        workflow.prepare(); workflow.smoke();
        assertEquals(3, workflow.arrayTasks());
        workflow.worker(0, 2); workflow.worker(1, 2);
        workflow.worker(0, 2); // completed partitions survive task resubmission
        Path cache = FingerprintFinalizer.finish(workflow).getParent();
        Path baseline = temp.resolve("baseline");
        Path build = temp.resolve("build.json");
        json(build, c.build(input, baseline, Path.of(c.modelBundle), Path.of(c.compactBundle)));
        MoleculeFingerprintIndexBuilderCLI.main(new String[]{"--config", build.toString()});
        var expected = records(baseline); var actual = records(cache);
        assertEquals(expected.size(), actual.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i).sourceRow(), actual.get(i).sourceRow());
            assertEquals(expected.get(i).moleculeId(), actual.get(i).moleculeId());
            assertEquals(expected.get(i).smiles(), actual.get(i).smiles());
            assertArrayEquals(expected.get(i).base128(), actual.get(i).base128(), 0.001f);
            assertArrayEquals(expected.get(i).compact16(), actual.get(i).compact16(), 0.001f);
        }
        assertEquals(screen(baseline, c, "baseline-search"), screen(cache, c, "partition-search"));
        try (var files = Files.list(Path.of(c.scratchRoot))) { assertEquals(0, files.count()); }
    }
    static List<MoleculeFingerprintRecord> records(Path cache) throws Exception {
        var manifest = MoleculeFingerprintIndexReader.loadManifest(cache);
        List<MoleculeFingerprintRecord> result = new ArrayList<>();
        for (int i = 0; i < manifest.shards.size(); i++)
            try (var reader = new MoleculeFingerprintIndexReader(cache, manifest, i)) { result.addAll(reader.readBatch(100)); }
        return result;
    }
    String screen(Path cache, FingerprintWorkflowConfig build, String name) throws Exception {
        var c = new MoleculeSimilaritySearchConfig();
        c.inputs.index = cache.toString(); c.inputs.modelBundle = build.modelBundle;
        c.query.structure = "CCCCCC"; c.query.identifier = "query";
        c.screening.resultTopK = 3; c.output.reportTopK = 3;
        c.runtime.device = "CPU"; c.runtime.comparatorBatchSize = 4;
        c.output.hitsTsv = temp.resolve(name + ".tsv").toString();
        c.output.hitsSdf = temp.resolve(name + ".sdf").toString();
        c.output.summaryMarkdown = temp.resolve(name + ".md").toString();
        c.output.runManifest = temp.resolve(name + "-manifest.json").toString();
        Path file = temp.resolve(name + ".json"); json(file, c);
        MoleculeSimilaritySearchCLI.main(new String[]{"--config", file.toString()});
        return Files.readString(Path.of(c.output.hitsTsv));
    }
}
