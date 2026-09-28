package com.idorsia.research.chem.hyperspace.gui2.model;

import static org.junit.Assert.*;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.view.*;

import org.junit.Test;

import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.*;
import javax.swing.event.TableModelEvent;

public class ExpansionLifecycleTest {
    private SynthonSpace.CombinatorialHit hit() {
        StereoMolecule molecule = new StereoMolecule();
        molecule.addAtom(6);
        return new SynthonSpace.CombinatorialHit(
                "test", Map.of(new SynthonSpace.FragType("test", 0), List.of(
                        new SynthonSpace.FragId("test", 0, molecule, "fixture", new BitSet(), new BitSet()))), null, Map.of());
    }

    private SynthonAssembler.ExpandedCombinatorialHit molecule(String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, smiles);
        return new SynthonAssembler.ExpandedCombinatorialHit(molecule.getIDCode(), List.of(
                new SynthonSpace.FragId("test", 0, molecule, "example", new BitSet(), new BitSet())));
    }

    private void add(CombinatorialSearchResultModel source) throws Exception {
        SwingUtilities.invokeAndWait(() -> source.addResults(List.of(hit())));
    }

    private Runnable take(BlockingQueue<Runnable> queue) throws Exception {
        Runnable update = queue.poll(10, TimeUnit.SECONDS);
        assertNotNull("Expected queued publication", update);
        return update;
    }

    private void close(RealTimeExpandingSearchResultModel model) throws Exception {
        SwingUtilities.invokeAndWait(model::dispose);
        assertTrue(
                "Expansion executors must terminate", model.awaitTermination(10, TimeUnit.SECONDS));
    }

    @Test
    public void replacedSearchRejectsQueuedOldPublication() throws Exception {
        var sourceA = new CombinatorialSearchResultModel(null);
        var sourceB = new CombinatorialSearchResultModel(null);
        var a = molecule("CC");
        var b = molecule("NN");
        BlockingQueue<Runnable> updatesA = new LinkedBlockingQueue<>();
        BlockingQueue<Runnable> updatesB = new LinkedBlockingQueue<>();
        var modelA =
                new RealTimeExpandingSearchResultModel(sourceA, 10, h -> List.of(a), updatesA::add);
        var modelB =
                new RealTimeExpandingSearchResultModel(sourceB, 10, h -> List.of(b), updatesB::add);
        try {
            add(sourceA);
            Runnable stale = take(updatesA);
            SwingUtilities.invokeAndWait(
                    () -> {
                        RealTimeExpandingHitsView view = new RealTimeExpandingHitsView(modelA);
                        view.setModel(modelB);
                        assertTrue(modelA.isDisposed());
                    });
            add(sourceB);
            SwingUtilities.invokeAndWait(take(updatesB));
            SwingUtilities.invokeAndWait(stale);
            SwingUtilities.invokeAndWait(
                    () -> {
                        assertEquals(0, modelA.getTableModel().getRowCount());
                        assertEquals(1, modelB.getTableModel().getRowCount());
                        assertSame(b, modelB.getTableModel().getValueAt(0, 1));
                    });
            add(sourceA);
            modelA.dispose();
            assertTrue(modelA.awaitTermination(10, TimeUnit.SECONDS));
            assertTrue(updatesA.isEmpty());
        } finally {
            close(modelA);
            close(modelB);
        }
    }

    @Test
    public void optionChangesRejectStaleGenerationsAndResetAllRows() throws Exception {
        var source = new CombinatorialSearchResultModel(null);
        var expanded = molecule("CO");
        BlockingQueue<Runnable> updates = new LinkedBlockingQueue<>();
        var model =
                new RealTimeExpandingSearchResultModel(
                        source, 10, h -> List.of(expanded), updates::add);
        try {
            add(source);
            Runnable first = take(updates);
            SwingUtilities.invokeAndWait(() -> model.setStructurePostprocessOptions(false, true));
            Runnable second = take(updates);
            SwingUtilities.invokeAndWait(() -> model.setStructurePostprocessOptions(true, false));
            Runnable latest = take(updates);
            SwingUtilities.invokeAndWait(latest);
            SwingUtilities.invokeAndWait(second);
            SwingUtilities.invokeAndWait(first);
            SwingUtilities.invokeAndWait(
                    () -> assertEquals(1, model.getTableModel().getRowCount()));
            SwingUtilities.invokeAndWait(
                    () -> {
                        model.setStructurePostprocessOptions(false, false);
                        assertEquals(0, model.getTableModel().getRowCount());
                    });
            SwingUtilities.invokeAndWait(take(updates));
            SwingUtilities.invokeAndWait(
                    () -> assertSame(expanded, model.getTableModel().getValueAt(0, 1)));
        } finally {
            close(model);
        }
    }

    @Test
    public void runningOldWorkCannotPublishAfterDisposal() throws Exception {
        var source = new CombinatorialSearchResultModel(null);
        CountDownLatch started = new CountDownLatch(1), release = new CountDownLatch(1);
        AtomicInteger publications = new AtomicInteger();
        var expanded = molecule("C");
        var model =
                new RealTimeExpandingSearchResultModel(
                        source,
                        10,
                        h -> {
                            started.countDown();
                            boolean done = false;
                            while (!done) {
                                try {
                                    release.await();
                                    done = true;
                                } catch (InterruptedException ignored) {
                                }
                            }
                            return List.of(expanded);
                        },
                        update -> publications.incrementAndGet());
        try {
            add(source);
            assertTrue(started.await(10, TimeUnit.SECONDS));
            SwingUtilities.invokeAndWait(model::dispose);
            add(source);
            release.countDown();
            assertTrue(model.awaitTermination(10, TimeUnit.SECONDS));
            assertEquals(0, publications.get());
        } finally {
            release.countDown();
            close(model);
        }
    }

    @Test
    public void rowLimitsAndInsertionEventsAreExactAndOnEdt() throws Exception {
        var source = new CombinatorialSearchResultModel(null);
        var a = molecule("C");
        var b = molecule("N");
        BlockingQueue<Runnable> updates = new LinkedBlockingQueue<>();
        var model =
                new RealTimeExpandingSearchResultModel(source, 1, h -> List.of(a, b), updates::add);
        List<TableModelEvent> events = new ArrayList<>();
        try {
            SwingUtilities.invokeAndWait(
                    () ->
                            model.getTableModel()
                                    .addTableModelListener(
                                            event -> {
                                                assertTrue(SwingUtilities.isEventDispatchThread());
                                                events.add(event);
                                            }));
            add(source);
            SwingUtilities.invokeAndWait(take(updates));
            SwingUtilities.invokeAndWait(
                    () -> {
                        assertEquals(1, model.getTableModel().getRowCount());
                        assertSame(a, model.getTableModel().getValueAt(0, 1));
                        assertTrue(
                                model.getTableModel()
                                        .getStructureData(0)
                                        .startsWith(a.assembled_idcode + " "));
                        assertEquals(1, events.size());
                        assertEquals(0, events.get(0).getFirstRow());
                        assertEquals(0, events.get(0).getLastRow());
                    });
        } finally {
            close(model);
        }
    }

    @Test
    public void expansionFailureIsReportedOnEdt() throws Exception {
        var source = new CombinatorialSearchResultModel(null);
        CountDownLatch reported = new CountDownLatch(1);
        var model =
                new RealTimeExpandingSearchResultModel(
                        source,
                        10,
                        h -> {
                            throw new IllegalStateException("test failure");
                        },
                        SwingUtilities::invokeAndWait);
        try {
            model.addListener(
                    () -> {
                        assertTrue(SwingUtilities.isEventDispatchThread());
                        reported.countDown();
                    });
            add(source);
            assertTrue(reported.await(10, TimeUnit.SECONDS));
            assertTrue(model.getResultsInfoString().contains("test failure"));
        } finally {
            close(model);
        }
    }

    private RealTimeExpandingHitsView findExpanded(Container root) {
        if (root instanceof RealTimeExpandingHitsView view) return view;
        for (Component component : root.getComponents()) {
            if (component instanceof Container container) {
                var found = findExpanded(container);
                if (found != null) return found;
            }
        }
        return null;
    }

    @Test
    public void replacingSelectedPreviewAndTopLevelViewsDisposesModels() throws Exception {
        org.junit.Assume.assumeFalse(GraphicsEnvironment.isHeadless());
        List<RealTimeExpandingSearchResultModel> models = new ArrayList<>();
        SwingUtilities.invokeAndWait(
                () -> {
                    var source = new CombinatorialSearchResultModel(null);
                    var view = new CombinatorialHitsView(source);
                    for (int i = 0; i < 5; i++) {
                        var old = findExpanded(view).getModel();
                        models.add(old);
                        view.setExpandedHits(List.of());
                        assertTrue(old.isDisposed());
                    }
                    models.add(findExpanded(view).getModel());
                    var root = new LeetHyperspaceView(new LeetHyperspaceModel());
                    root.addResultsView(view);
                    root.addResultsView(
                            new CombinatorialHitsView(new CombinatorialSearchResultModel(null)));
                    assertTrue(models.get(models.size() - 1).isDisposed());
                    var expanded = new RealTimeExpandingSearchResultModel(source, 10);
                    models.add(expanded);
                    root.addResultsView(new RealTimeExpandingHitsView(expanded));
                    root.disposeResults();
                    root.disposeResults();
                    assertTrue(expanded.isDisposed());
                });
        for (var model : models) close(model);
    }

    @Test
    public void emptyBatchDoesNotEmitAnInsertion() throws Exception {
        var source = new CombinatorialSearchResultModel(null);
        BlockingQueue<Runnable> updates = new LinkedBlockingQueue<>();
        var model =
                new RealTimeExpandingSearchResultModel(source, 10, h -> List.of(), updates::add);
        AtomicInteger events = new AtomicInteger();
        try {
            SwingUtilities.invokeAndWait(
                    () ->
                            model.getTableModel()
                                    .addTableModelListener(e -> events.incrementAndGet()));
            add(source);
            SwingUtilities.invokeAndWait(take(updates));
            assertEquals(0, events.get());
        } finally {
            close(model);
        }
    }

    @Test
    public void desktopSwitchAndWindowCloseSmoke() throws Exception {
        org.junit.Assume.assumeTrue(
                Boolean.getBoolean("hyperspace.gui2Smoke") && !GraphicsEnvironment.isHeadless());
        JFrame[] frame = new JFrame[1];
        var sourceA = new CombinatorialSearchResultModel(null);
        var sourceB = new CombinatorialSearchResultModel(null);
        var a = molecule("c1ccccc1");
        var b = molecule("CC(=O)Nc1ccc(O)cc1");
        BlockingQueue<Runnable> updatesA = new LinkedBlockingQueue<>(),
                updatesB = new LinkedBlockingQueue<>();
        var modelA =
                new RealTimeExpandingSearchResultModel(sourceA, 10, h -> List.of(a), updatesA::add);
        var modelB =
                new RealTimeExpandingSearchResultModel(sourceB, 10, h -> List.of(b), updatesB::add);
        CountDownLatch closed = new CountDownLatch(1);
        try {
            SwingUtilities.invokeAndWait(
                    () -> {
                        frame[0] =
                                com.idorsia.research.chem.hyperspace.gui2.LeetHyperspaceMain.open(
                                        null);
                        frame[0].addWindowListener(
                                new java.awt.event.WindowAdapter() {
                                    @Override
                                    public void windowClosed(java.awt.event.WindowEvent e) {
                                        closed.countDown();
                                    }
                                });
                        var view = (LeetHyperspaceView) frame[0].getContentPane().getComponent(0);
                        view.addResultsView(new RealTimeExpandingHitsView(modelA));
                    });
            add(sourceA);
            Runnable stale = take(updatesA);
            SwingUtilities.invokeAndWait(
                    () -> {
                        var view = (LeetHyperspaceView) frame[0].getContentPane().getComponent(0);
                        view.addResultsView(new RealTimeExpandingHitsView(modelB));
                        selectEnumerated(view);
                    });
            add(sourceB);
            SwingUtilities.invokeAndWait(take(updatesB));
            SwingUtilities.invokeAndWait(stale);
            SwingUtilities.invokeAndWait(() -> assertTrue(clickOption(frame[0], "Highlight")));
            SwingUtilities.invokeAndWait(take(updatesB));
            SwingUtilities.invokeAndWait(() -> assertTrue(clickOption(frame[0], "Align")));
            SwingUtilities.invokeAndWait(take(updatesB));
            SwingUtilities.invokeAndWait(
                    () -> {
                        assertEquals(1, modelB.getTableModel().getRowCount());
                        assertSame(b, modelB.getTableModel().getValueAt(0, 1));
                        java.awt.image.BufferedImage image =
                                new java.awt.image.BufferedImage(
                                        frame[0].getWidth(),
                                        frame[0].getHeight(),
                                        java.awt.image.BufferedImage.TYPE_INT_RGB);
                        Graphics2D graphics = image.createGraphics();
                        frame[0].printAll(graphics);
                        graphics.dispose();
                        try {
                            javax.imageio.ImageIO.write(
                                    image,
                                    "png",
                                    new java.io.File("/tmp/gui2-hardened-results.png"));
                        } catch (java.io.IOException ex) {
                            throw new RuntimeException(ex);
                        }
                        frame[0].dispose();
                    });
            assertTrue(closed.await(10, TimeUnit.SECONDS));
            assertTrue(modelA.isDisposed());
            assertTrue(modelB.isDisposed());
        } finally {
            SwingUtilities.invokeAndWait(
                    () -> {
                        if (frame[0] != null) frame[0].dispose();
                    });
            close(modelA);
            close(modelB);
        }
    }

    private boolean clickOption(Container container, String name) {
        if (container instanceof JCheckBox checkbox && name.equals(checkbox.getText())) {
            checkbox.doClick();
            return true;
        }
        for (Component component : container.getComponents()) {
            if (component instanceof Container child && clickOption(child, name)) return true;
        }
        return false;
    }

    private void selectEnumerated(Container container) {
        if (container instanceof JTabbedPane tabs) {
            for (int i = 0; i < tabs.getTabCount(); i++) {
                if ("Enumerated".equals(tabs.getTitleAt(i))) tabs.setSelectedIndex(i);
            }
        }
        for (Component component : container.getComponents()) {
            if (component instanceof Container child) selectEnumerated(child);
        }
    }
}
