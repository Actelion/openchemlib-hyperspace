package com.idorsia.research.chem.hyperspace.gui2.model;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.coords.InventorTemplate;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.*;
import java.util.concurrent.*;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;

/** Owns bounded expansion work for one search. Disposal is permanent. */
public class RealTimeExpandingSearchResultModel implements AutoCloseable {
    private final CombinatorialSearchResultModel resultModel;
    private final StereoMolecule query;
    private final int maxExpandedHits;
    private final int workers =
            Math.max(1, Math.min(4, Runtime.getRuntime().availableProcessors()));
    private final ExecutorService coordinator =
            Executors.newSingleThreadExecutor(threadFactory("coordinator"));
    private final ThreadPoolExecutor expansionPool =
            (ThreadPoolExecutor) Executors.newFixedThreadPool(workers, threadFactory("worker"));
    private final Object lifecycle = new Object();
    private final Set<Future<?>> tasks = new HashSet<>();
    private final Expansion expansion;
    private final Publication publication;
    private final CombinatorialSearchResultModel.CombinatorialSearchResultModelListener
            sourceListener = this::requestExpansion;
    private final List<RealTimeExpandingSearchResultModelListener> listeners =
            new CopyOnWriteArrayList<>();
    private final RealTimeExpandingTableModel tableModel = new RealTimeExpandingTableModel();
    // Only the EDT accesses rows. Workers read the separately published count.
    private final List<Row> rows = new ArrayList<>();
    private volatile int shown;
    private volatile long generation;
    private volatile boolean disposed;
    private volatile boolean highlightSubstructure = true;
    private volatile boolean alignSubstructure;
    private volatile String error;
    private boolean requested;
    private boolean draining;

    @FunctionalInterface
    interface Expansion {
        List<SynthonAssembler.ExpandedCombinatorialHit> expand(SynthonSpace.CombinatorialHit hit)
                throws Exception;
    }

    @FunctionalInterface
    interface Publication {
        void publish(Runnable update) throws Exception;
    }

    private record Settings(long generation, boolean highlight, boolean align) {}

    private record Row(String structure, SynthonAssembler.ExpandedCombinatorialHit hit) {}

    public RealTimeExpandingSearchResultModel(CombinatorialSearchResultModel model, int limit) {
        this(
                model,
                limit,
                hit -> SynthonAssembler.expandCombinatorialHit(hit, 1024),
                SwingUtilities::invokeAndWait);
    }

    // Test hooks allow deterministic interleaving without timing-dependent chemistry.
    RealTimeExpandingSearchResultModel(
            CombinatorialSearchResultModel model,
            int limit,
            Expansion expansion,
            Publication publication) {
        if (limit < 1) throw new IllegalArgumentException("Expansion limit must be positive");
        resultModel = Objects.requireNonNull(model);
        query = model.getQuery() == null ? null : new StereoMolecule(model.getQuery());
        maxExpandedHits = limit;
        this.expansion = expansion;
        this.publication = publication;
        model.addListener(sourceListener);
        requestExpansion();
    }

    private static ThreadFactory threadFactory(String role) {
        return runnable -> {
            Thread thread = new Thread(runnable, "gui2-expansion-" + role);
            thread.setPriority(Thread.MIN_PRIORITY);
            thread.setDaemon(true);
            return thread;
        };
    }

    private boolean current(long token) {
        return !disposed && generation == token;
    }

    private void requestExpansion() {
        synchronized (lifecycle) {
            if (disposed) return;
            requested = true;
            if (!draining) {
                draining = true;
                coordinator.execute(this::drain);
            }
        }
    }

    private void drain() {
        Set<SynthonSpace.CombinatorialHit> processed =
                Collections.newSetFromMap(new IdentityHashMap<>());
        long processedGeneration = -1;
        while (true) {
            Settings settings;
            synchronized (lifecycle) {
                if (disposed || !requested) {
                    draining = false;
                    return;
                }
                requested = false;
                settings = new Settings(generation, highlightSubstructure, alignSubstructure);
            }
            if (processedGeneration != settings.generation()) {
                processed.clear();
                processedGeneration = settings.generation();
            }
            try {
                expandAvailable(settings, processed);
            } catch (InterruptedException ex) {
                Thread.currentThread().interrupt();
                return;
            } catch (Exception ex) {
                if (current(settings.generation())) {
                    ex.printStackTrace(System.err);
                    SwingUtilities.invokeLater(
                            () -> {
                                if (!current(settings.generation())) return;
                                error = "Expansion failed: " + ex.getMessage();
                                fireResultsChanged();
                            });
                }
            }
        }
    }

    // At most 'workers' batches are submitted. Waiting for EDT publication also bounds UI
    // callbacks.
    private void expandAvailable(Settings settings, Set<SynthonSpace.CombinatorialHit> processed)
            throws Exception {
        Deque<Future<List<Row>>> pending = new ArrayDeque<>();
        try {
            for (SynthonSpace.CombinatorialHit hit : resultModel.getHits()) {
                if (!current(settings.generation()) || shown >= maxExpandedHits) break;
                if (!processed.add(hit)) continue;
                for (SynthonSpace.CombinatorialHit chunk : splitHit(hit)) {
                    if (!current(settings.generation()) || shown >= maxExpandedHits) break;
                    synchronized (lifecycle) {
                        if (!current(settings.generation())) break;
                        Future<List<Row>> task =
                                expansionPool.submit(() -> expandChunk(chunk, settings));
                        tasks.add(task);
                        pending.addLast(task);
                    }
                    if (pending.size() >= workers) publishNext(pending.removeFirst(), settings);
                }
            }
            while (!pending.isEmpty() && current(settings.generation()))
                publishNext(pending.removeFirst(), settings);
        } finally {
            synchronized (lifecycle) {
                for (Future<?> task : pending) {
                    task.cancel(true);
                    tasks.remove(task);
                }
                expansionPool.purge();
            }
        }
    }

    private List<Row> expandChunk(SynthonSpace.CombinatorialHit hit, Settings settings)
            throws Exception {
        if (!current(settings.generation()) || Thread.currentThread().isInterrupted())
            return List.of();
        List<Row> batch = new ArrayList<>();
        StereoMolecule localQuery = query == null ? null : new StereoMolecule(query);
        CoordinateInventor inventor = new CoordinateInventor();
        if (settings.align() && localQuery != null) {
            inventor.setCustomTemplateList(
                    Collections.singletonList(
                            new InventorTemplate(
                                    localQuery,
                                    new DescriptorHandlerLongFFP512().createDescriptor(localQuery),
                                    true)));
        }
        for (SynthonAssembler.ExpandedCombinatorialHit expanded : expansion.expand(hit)) {
            if (!current(settings.generation()) || Thread.currentThread().isInterrupted()) break;
            StereoMolecule molecule = HyperspaceUtils.parseIDCode(expanded.assembled_idcode);
            inventor.invent(molecule);
            if (settings.highlight() && localQuery != null)
                HyperspaceUtils.setHighlightedSubstructure(molecule, localQuery);
            Canonizer canonizer = new Canonizer(molecule);
            batch.add(
                    new Row(
                            canonizer.getIDCode() + " " + canonizer.getEncodedCoordinates(),
                            expanded));
        }
        return batch;
    }

    private void publishNext(Future<List<Row>> task, Settings settings) throws Exception {
        try {
            List<Row> batch = task.get();
            if (!current(settings.generation())) return;
            publication.publish(
                    () -> {
                        if (!current(settings.generation())) return;
                        int oldSize = rows.size();
                        int count = Math.min(batch.size(), maxExpandedHits - oldSize);
                        if (count == 0) return;
                        rows.addAll(batch.subList(0, count));
                        shown = rows.size();
                        tableModel.fireTableRowsInserted(oldSize, rows.size() - 1);
                        fireResultsChanged();
                    });
        } catch (CancellationException ex) {
            if (current(settings.generation())) throw ex;
        } finally {
            synchronized (lifecycle) {
                tasks.remove(task);
            }
        }
    }

    public void setStructurePostprocessOptions(boolean highlight, boolean align) {
        if (!SwingUtilities.isEventDispatchThread()) {
            SwingUtilities.invokeLater(() -> setStructurePostprocessOptions(highlight, align));
            return;
        }
        synchronized (lifecycle) {
            if (disposed || (highlight == highlightSubstructure && align == alignSubstructure))
                return;
            generation++;
            highlightSubstructure = highlight;
            alignSubstructure = align;
            cancelTasks();
        }
        rows.clear();
        shown = 0;
        error = null;
        tableModel.fireTableDataChanged();
        fireResultsChanged();
        requestExpansion();
    }

    private void cancelTasks() {
        for (Future<?> task : tasks) task.cancel(true);
        tasks.clear();
        expansionPool.purge();
    }

    public void dispose() {
        synchronized (lifecycle) {
            if (disposed) return;
            disposed = true;
            generation++;
            resultModel.removeListener(sourceListener);
            cancelTasks();
            coordinator.shutdownNow();
            expansionPool.shutdownNow();
            listeners.clear();
        }
    }

    @Override
    public void close() {
        dispose();
    }

    public boolean isDisposed() {
        return disposed;
    }

    boolean isTerminated() {
        return coordinator.isTerminated() && expansionPool.isTerminated();
    }

    boolean awaitTermination(long timeout, TimeUnit unit) throws InterruptedException {
        long deadline = System.nanoTime() + unit.toNanos(timeout);
        return coordinator.awaitTermination(timeout, unit)
                && expansionPool.awaitTermination(
                        Math.max(0, deadline - System.nanoTime()), TimeUnit.NANOSECONDS);
    }

    public boolean isHighlightSubstructure() {
        return highlightSubstructure;
    }

    public boolean isAlignSubstructure() {
        return alignSubstructure;
    }

    public CombinatorialSearchResultModel getCombinatorialSearchResultModel() {
        return resultModel;
    }

    public RealTimeExpandingTableModel getTableModel() {
        return tableModel;
    }

    public void addListener(RealTimeExpandingSearchResultModelListener listener) {
        synchronized (lifecycle) {
            if (!disposed) listeners.add(listener);
        }
    }

    public boolean removeListener(RealTimeExpandingSearchResultModelListener listener) {
        return listeners.remove(listener);
    }

    private void fireResultsChanged() {
        for (var listener : listeners) listener.resultsChanged();
    }

    public String getResultsInfoString() {
        long total =
                resultModel.getHits().stream()
                        .mapToLong(
                                hit ->
                                        hit.hit_fragments.values().stream()
                                                .mapToLong(List::size)
                                                .reduce(1, (a, b) -> a * b))
                        .sum();
        return error == null ? String.format("Results: %6d  Showing: %6d", total, shown) : error;
    }

    public interface RealTimeExpandingSearchResultModelListener {
        void resultsChanged();
    }

    public class RealTimeExpandingTableModel extends AbstractTableModel {
        public String getStructureData(int row) {
            return rows.get(row).structure();
        }

        @Override
        public int getRowCount() {
            return rows.size();
        }

        @Override
        public int getColumnCount() {
            return 2;
        }

        @Override
        public String getColumnName(int column) {
            return column == 0 ? "Structure" : "Synthons";
        }

        @Override
        public Object getValueAt(int row, int column) {
            return column == 0 ? rows.get(row).structure() : rows.get(row).hit();
        }
    }

    private static Iterable<SynthonSpace.CombinatorialHit> splitHit(
            SynthonSpace.CombinatorialHit hit) {
        List<SynthonSpace.FragType> types = new ArrayList<>(hit.hit_fragments.keySet());
        types.sort(Comparator.comparingInt(type -> type.frag));
        long count = 1;
        SynthonSpace.FragType longest = null;
        for (var type : types) {
            int size = hit.hit_fragments.get(type).size();
            count = size == 0 ? 0 : count > Long.MAX_VALUE / size ? Long.MAX_VALUE : count * size;
            if (longest == null || size > hit.hit_fragments.get(longest).size()) longest = type;
        }
        if (count <= 200) return List.of(hit);
        var axis = longest;
        var fragments = hit.hit_fragments.get(axis);
        int chunkSize = Math.max(1, (int) (fragments.size() / (count / 100.0)));
        return () ->
                new Iterator<>() {
                    private int start;

                    public boolean hasNext() {
                        return start < fragments.size();
                    }

                    public SynthonSpace.CombinatorialHit next() {
                        if (!hasNext()) throw new NoSuchElementException();
                        int end = Math.min(fragments.size(), start + chunkSize);
                        Map<SynthonSpace.FragType, List<SynthonSpace.FragId>> sets =
                                new HashMap<>(hit.hit_fragments);
                        sets.put(axis, new ArrayList<>(fragments.subList(start, end)));
                        start = end;
                        return new SynthonSpace.CombinatorialHit(
                                hit.rxn, sets, hit.sri, hit.mapping);
                    }
                };
    }
}
