package com.idorsia.research.chem.hyperspace.gui;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.gui.process.*;
import com.idorsia.research.chem.hyperspace.gui.search.*;

import org.json.JSONObject;
import org.junit.jupiter.api.Test;

import java.nio.file.*;
import java.util.*;

import javax.swing.*;

class PreparedSpaceSmokeTest {
    @Test
    void preparedArtifactsLoadAndSearch() throws Exception {
        String manifest = System.getProperty("hyperspace.testManifest");
        assumeTrue(manifest != null);
        JSONObject artifacts =
                new JSONObject(Files.readString(Path.of(manifest))).getJSONObject("artifacts");
        var raw = HyperspaceIOUtils.loadRawSynthonSpace(artifacts.getString("rawspace"));
        var reaction = raw.getReactions().values().iterator().next();
        List<StereoMolecule> parts = new ArrayList<>();
        reaction.getRawFragmentSets().entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .forEach(
                        e ->
                                parts.add(
                                        new IDCodeParser()
                                                .getCompactMolecule(
                                                        e.getValue().get(0).getIdcode())));
        StereoMolecule query = SynthonAssembler.assembleSynthons_faster(parts);
        var providers =
                HyperspaceInit.loadSearchProviderInitFile(
                        null, Path.of(artifacts.getString("guiConfig")).toFile());
        for (AbstractSearchProvider provider : providers.values()) {
            AbstractHyperspaceProcess load = provider.startInitialization();
            assertTrue(load.waitUntilDoneOrFailed(30000));
            SwingUtilities.invokeAndWait(() -> {});
            assertEquals(AbstractSearchProvider.SearchProviderStatus.READY, provider.getStatus());
            AbstractSearchProvider.SearchConfiguration config =
                    provider instanceof HyperspaceSimilaritySearch
                            ? new HyperspaceSimilaritySearch.SimilaritySearchConfiguration(
                                    query, 3, 3, 100, 2)
                            : new HyperspaceSubstructureSearch.SubstructureSearchConfiguration(
                                    query, true, 3, 20, 2);
            AbstractHyperspaceSearchProcess search = provider.runSearch(config);
            assertTrue(search.waitUntilDoneOrFailed(60000));
            assertEquals(AbstractHyperspaceProcess.ProcessStatus.DONE, search.getProcessStatus());
            assertFalse(search.getSearchResults().isEmpty());
        }
    }

    @Test
    void guiOpensPreparedConfiguration() throws Exception {
        String manifest = System.getProperty("hyperspace.testManifest");
        assumeTrue(manifest != null && Boolean.getBoolean("hyperspace.guiSmoke"));
        assumeTrue(!java.awt.GraphicsEnvironment.isHeadless());
        String config =
                new JSONObject(Files.readString(Path.of(manifest)))
                        .getJSONObject("artifacts")
                        .getString("guiConfig");
        HyperspaceSearchGUI.main(new String[] {config});
        try {
            long until = System.nanoTime() + java.util.concurrent.TimeUnit.SECONDS.toNanos(30);
            boolean found = false;
            while (System.nanoTime() < until) {
                SwingUtilities.invokeAndWait(() -> {});
                for (java.awt.Frame f : java.awt.Frame.getFrames())
                    if (f.isVisible() && f.getTitle().contains(config)) found = true;
                if (found) break;
                Thread.sleep(100);
            }
            assertTrue(found);
            Thread.sleep(2000);
            for (java.awt.Frame f : java.awt.Frame.getFrames())
                if (f.isVisible() && f.getTitle().contains(config)) {
                    var image =
                            new java.awt.image.BufferedImage(
                                    f.getWidth(),
                                    f.getHeight(),
                                    java.awt.image.BufferedImage.TYPE_INT_RGB);
                    SwingUtilities.invokeAndWait(
                            () -> {
                                var graphics = image.createGraphics();
                                try {
                                    f.printAll(graphics);
                                } finally {
                                    graphics.dispose();
                                }
                            });
                    javax.imageio.ImageIO.write(
                            image, "png", Path.of("/tmp/hyperspace-mcp-gui.png").toFile());
                    if (Boolean.getBoolean("hyperspace.screenCapture")) {
                        javax.imageio.ImageIO.write(new java.awt.Robot().createScreenCapture(f.getBounds()),
                                "png", Path.of("/tmp/hyperspace-mcp-gui-screen.png").toFile());
                    }
                }
        } finally {
            SwingUtilities.invokeAndWait(
                    () -> {
                        for (java.awt.Window w : java.awt.Window.getWindows()) w.dispose();
                    });
        }
    }
}
