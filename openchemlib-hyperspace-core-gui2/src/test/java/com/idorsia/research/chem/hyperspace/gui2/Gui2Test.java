package com.idorsia.research.chem.hyperspace.gui2;

import static org.junit.Assert.*;

import com.actelion.research.chem.*;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.*;
import com.idorsia.research.chem.hyperspace.gui2.task.*;

import org.junit.Test;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.zip.GZIPOutputStream;

import javax.swing.*;

public class Gui2Test {
    private Map<String, Object> provider(String type, String file) {
        return Map.of(
                "ServiceProvider",
                type,
                "ServiceName",
                "Supplier test",
                "Config",
                Map.of("File", file, "MaxNumberOfThreads", 2));
    }

    @Test
    public void configSupportsRelativePathsAndWarnsAboutSimilarity() throws Exception {
        Path root = Files.createTempDirectory("gui2-config");
        Files.write(root.resolve("index.data"), new byte[] {1});
        Path json = root.resolve("gui.json");
        new ObjectMapper()
                .writeValue(
                        json.toFile(),
                        Map.of(
                                "ServiceProviders",
                                List.of(
                                        provider("HyperspaceSSS", "index.data"),
                                        provider("HyperspaceSimilarity", "unused.data"))));
        var config = Gui2Configuration.read(json);
        assertEquals(1, config.spaces().size());
        assertEquals(root.resolve("index.data"), config.spaces().get(0).file());
        assertEquals("Supplier test", config.spaces().get(0).name());
        assertEquals(2, config.spaces().get(0).threads());
        assertEquals(1, config.warnings().size());
        Path lines = root.resolve("spaces.conf");
        Files.writeString(lines, "# spaces\nindex.data\n");
        assertEquals(
                root.resolve("index.data"), Gui2Configuration.read(lines).spaces().get(0).file());
        new ObjectMapper()
                .writeValue(
                        json.toFile(),
                        Map.of(
                                "ServiceProviders",
                                List.of(provider("HyperspaceSimilarity", "unused.data"))));
        try {
            Gui2Configuration.read(json);
            fail("Expected similarity-only rejection");
        } catch (IOException ex) {
            assertTrue(ex.getMessage().contains("legacy"));
        }
        new ObjectMapper()
                .writeValue(
                        json.toFile(),
                        Map.of(
                                "ServiceProviders",
                                List.of(provider("HyperspaceSSS", "missing.data"))));
        try {
            Gui2Configuration.read(json);
            fail("Expected missing file rejection");
        } catch (IOException ex) {
            assertTrue(ex.getMessage().contains("readable"));
        }
    }

    private Path toyIndex() throws Exception {
        SynthonSpace space = new SynthonSpace();
        space.setFP(SynthonSpace.resolveDescriptorHandlerFromName("FragFp"), 512);
        Map<Integer, List<Object>> sets = new HashMap<>();
        Map<String, String> ids = new HashMap<>();
        int i = 0;
        for (String smiles : List.of("C[U]", "N[U]")) {
            StereoMolecule molecule = new StereoMolecule();
            new SmilesParser().parse(molecule, smiles);
            String idcode = molecule.getIDCode();
            sets.put(i, List.of(idcode));
            ids.put(idcode, "synthon-" + i++);
        }
        space.addReaction("toy", sets, ids, null);
        space.initAfterJavaDeserialization();
        space.reinitHelperMaps();
        space.reinitBitTree();
        Path file = Files.createTempFile("gui2-toy", ".data");
        try (ObjectOutputStream out =
                new ObjectOutputStream(new GZIPOutputStream(Files.newOutputStream(file)))) {
            out.writeObject(space);
        }
        return file;
    }

    private void waitFor(java.util.function.BooleanSupplier condition) throws Exception {
        long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(10);
        while (System.nanoTime() < deadline) {
            boolean[] ok = {false};
            SwingUtilities.invokeAndWait(() -> ok[0] = condition.getAsBoolean());
            if (ok[0]) return;
            Thread.sleep(30);
        }
        fail("Timed out waiting for EDT update");
    }

    @Test
    public void loadsAndSearchesWithConfiguredThreads() throws Exception {
        LeetHyperspaceModel model = new LeetHyperspaceModel();
        LoadSynthonSpaceTask load =
                new LoadSynthonSpaceTask(model, toyIndex().toString(), "Toy", 2);
        load.execute();
        SynthonSpace space = load.get(30, TimeUnit.SECONDS);
        waitFor(() -> !model.getSynthonSpaces().isEmpty());
        assertEquals(100, load.getProgress());
        assertEquals("Toy", model.getSynthonSpaces().get(0).getName());
        assertEquals(2, model.getSynthonSpaces().get(0).getThreads());
        assertFalse(
                Thread.getAllStackTraces().keySet().stream()
                        .anyMatch(t -> t.isAlive() && t.getName().equals("gui2-load-progress")));
        StereoMolecule query = new StereoMolecule();
        new SmilesParser().parse(query, "CN");
        query.setFragment(true);
        CombinatorialSearchResultModel results = new CombinatorialSearchResultModel(query);
        SubstructureSearchTask search = new SubstructureSearchTask(space, query, results, 2);
        search.execute();
        search.get(30, TimeUnit.SECONDS);
        waitFor(() -> !results.getHits().isEmpty());
        CombinatorialSearchResultModel batches = new CombinatorialSearchResultModel(query);
        class BatchTask extends SubstructureSearchTask {
            BatchTask() {
                super(space, query, batches, 2);
            }

            void feed(List<SynthonSpace.CombinatorialHit> hits) {
                process(List.of(hits));
            }
        }
        BatchTask task = new BatchTask();
        var hits = results.getHits();
        SwingUtilities.invokeAndWait(
                () -> {
                    task.feed(hits);
                    task.feed(hits);
                });
        assertEquals(2 * hits.size(), batches.getHits().size());
    }

    @Test
    public void failedLoadDoesNotPublishSpace() throws Exception {
        org.junit.Assume.assumeTrue(java.awt.GraphicsEnvironment.isHeadless());
        Path file = Files.createTempFile("gui2-bad", ".data");
        Files.writeString(file, "not an index");
        LeetHyperspaceModel model = new LeetHyperspaceModel();
        LoadSynthonSpaceTask load = new LoadSynthonSpaceTask(model, file.toString());
        load.execute();
        try {
            load.get(10, TimeUnit.SECONDS);
            fail("Expected failure");
        } catch (ExecutionException expected) {
        }
        waitFor(() -> load.getName().startsWith("Failed:"));
        assertTrue(model.getSynthonSpaces().isEmpty());
        assertTrue(load.getProgress() < 100);
    }

    @Test
    public void renderGuiWhenRequested() throws Exception {
        org.junit.Assume.assumeTrue(Boolean.getBoolean("hyperspace.gui2Smoke"));
        var configuration =
                new Gui2Configuration.Configuration(
                        List.of(
                                new Gui2Configuration.Space(
                                        toyIndex(), "Toy substructure space", 2)),
                        List.of());
        JFrame[] frame = new JFrame[1];
        try {
            SwingUtilities.invokeAndWait(() -> frame[0] = LeetHyperspaceMain.open(configuration));
            Thread.sleep(2000);
            SwingUtilities.invokeAndWait(
                    () -> {
                        try {
                            java.awt.image.BufferedImage image =
                                    new java.awt.image.BufferedImage(
                                            frame[0].getWidth(),
                                            frame[0].getHeight(),
                                            java.awt.image.BufferedImage.TYPE_INT_RGB);
                            java.awt.Graphics2D graphics = image.createGraphics();
                            frame[0].printAll(graphics);
                            graphics.dispose();
                            javax.imageio.ImageIO.write(
                                    image, "png", new File("/tmp/hyperspace-gui2.png"));
                        } catch (IOException ex) {
                            throw new RuntimeException(ex);
                        }
                    });
        } finally {
            SwingUtilities.invokeAndWait(
                    () -> {
                        if (frame[0] != null) frame[0].dispose();
                    });
        }
    }
}
