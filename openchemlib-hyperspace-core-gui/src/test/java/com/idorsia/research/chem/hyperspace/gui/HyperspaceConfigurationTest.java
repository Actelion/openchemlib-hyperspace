package com.idorsia.research.chem.hyperspace.gui;

import static org.junit.jupiter.api.Assertions.*;

import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.gui.search.*;

import org.json.*;
import org.junit.jupiter.api.Test;

import java.nio.file.*;
import java.util.*;

import javax.swing.SwingUtilities;

class HyperspaceConfigurationTest {
    @Test
    void savesBothModesAndPreservesNamesAndRelativePaths() throws Exception {
        Path root = Files.createTempDirectory("gui-config");
        Files.writeString(root.resolve("sss.data"), "not loaded");
        Files.writeString(root.resolve("sim.data"), "not loaded");
        JSONArray providers = new JSONArray();
        for (String type : List.of("HyperspaceSSS", "HyperspaceSimilarity"))
            providers.put(
                    new JSONObject()
                            .put("ServiceName", "Named " + type)
                            .put("ServiceProvider", type)
                            .put(
                                    "Config",
                                    new JSONObject()
                                            .put("SpaceName", "Toy")
                                            .put(
                                                    "File",
                                                    type.equals("HyperspaceSSS")
                                                            ? "sss.data"
                                                            : "sim.data")
                                            .put("MaxNumberOfThreads", 2)));
        Path config = root.resolve("gui.json");
        Files.writeString(config, new JSONObject().put("ServiceProviders", providers).toString());
        Map<String, AbstractSearchProvider> loaded =
                HyperspaceInit.loadSearchProviderInitFile(null, config.toFile());
        assertEquals(2, loaded.size());
        Path saved = root.resolve("saved.json");
        HyperspaceInit.saveSearchProviders(new ArrayList<>(loaded.values()), saved.toFile());
        Map<String, AbstractSearchProvider> roundtrip =
                HyperspaceInit.loadSearchProviderInitFile(null, saved.toFile());
        assertEquals(loaded.keySet(), roundtrip.keySet());
        assertTrue(
                roundtrip.get("Named HyperspaceSimilarity") instanceof HyperspaceSimilaritySearch);
        JSONObject c =
                new JSONObject(
                        roundtrip
                                .get("Named HyperspaceSimilarity")
                                .getSearchProviderConfiguration()
                                .serializeToJSON());
        assertEquals(root.resolve("sim.data").toString(), c.getString("File"));
        Files.delete(root.resolve("sim.data"));
        assertThrows(
                java.io.IOException.class,
                () -> HyperspaceInit.loadSearchProviderInitFile(null, saved.toFile()));
    }

    @Test
    void corruptIndexPropagatesFailedStateAndMessage() throws Exception {
        Path bad = Files.createTempFile("bad-space", ".data");
        Files.writeString(bad, "corrupt");
        for (AbstractSearchProvider provider :
                List.of(new HyperspaceSubstructureSearch(), new HyperspaceSimilaritySearch())) {
            provider.setConfigurationAndGUI(
                    new HyperspaceSubstructureSearch.InitializationConfig(
                            "Broken", bad.toString(), 1),
                    null);
            AbstractHyperspaceProcess process = provider.startInitialization();
            assertTrue(process.waitUntilDoneOrFailed(5000));
            SwingUtilities.invokeAndWait(() -> {});
            assertEquals(
                    AbstractHyperspaceProcess.ProcessStatus.FAILED, process.getProcessStatus());
            assertEquals(AbstractSearchProvider.SearchProviderStatus.ERROR, provider.getStatus());
            assertTrue(process.getProcessStatusMessage().contains(bad.toString()));
        }
    }

    @Test
    void alreadyFinishedProcessDoesNotWaitUntilTimeout() throws Exception {
        class Finished extends AbstractHyperspaceProcess {
            public String getName() {
                return "test";
            }

            void done() {
                setProcessStatus(ProcessStatus.DONE);
            }
        }
        Finished p = new Finished();
        p.done();
        SwingUtilities.invokeAndWait(() -> {});
        assertTrue(p.waitUntilDoneOrFailed(1));
    }
}
