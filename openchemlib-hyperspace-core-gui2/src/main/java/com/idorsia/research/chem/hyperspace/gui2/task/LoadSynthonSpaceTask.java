package com.idorsia.research.chem.hyperspace.gui2.task;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.model.LoadedSynthonSpace;

import org.apache.commons.io.input.CountingInputStream;

import java.awt.GraphicsEnvironment;
import java.io.*;
import java.nio.file.*;
import java.util.concurrent.*;
import java.util.zip.GZIPInputStream;

import javax.swing.*;

public class LoadSynthonSpaceTask extends SwingWorker<SynthonSpace, Double>
        implements HyperspaceTask {
    private final LeetHyperspaceModel model;
    private final String filepath;
    private final String name;
    private final int threads;
    private volatile String failure;

    public LoadSynthonSpaceTask(LeetHyperspaceModel model, String filepath) {
        this(
                model,
                filepath,
                Path.of(filepath).getFileName().toString(),
                Runtime.getRuntime().availableProcessors());
    }

    public LoadSynthonSpaceTask(
            LeetHyperspaceModel model, String filepath, String name, int threads) {
        if (threads < 1) throw new IllegalArgumentException("threads must be positive");
        this.model = model;
        this.filepath = filepath;
        this.name = name;
        this.threads = threads;
    }

    @Override
    public String getName() {
        return failure == null ? "Load space " + name : "Failed: " + name + " - " + failure;
    }

    public SwingWorker getThisWorker() {
        return this;
    }

    @Override
    protected SynthonSpace doInBackground() throws Exception {
        long size = Files.size(Path.of(filepath));
        try (CountingInputStream input =
                        new CountingInputStream(
                                new BufferedInputStream(Files.newInputStream(Path.of(filepath))));
                ObjectInputStream objects = new ObjectInputStream(new GZIPInputStream(input))) {
            ScheduledExecutorService reporter =
                    Executors.newSingleThreadScheduledExecutor(
                            r -> new Thread(r, "gui2-load-progress"));
            try {
                reporter.scheduleAtFixedRate(
                        () ->
                                setProgress(
                                        (int)
                                                Math.min(
                                                        80,
                                                        80.0
                                                                * input.getByteCount()
                                                                / Math.max(1, size))),
                        0,
                        250,
                        TimeUnit.MILLISECONDS);
                Object value = objects.readObject();
                if (!(value instanceof SynthonSpace))
                    throw new IOException("Not a substructure SynthonSpace index");
                SynthonSpace space = (SynthonSpace) value;
                space.initAfterJavaDeserialization();
                return space;
            } finally {
                reporter.shutdownNow();
                boolean interrupted = Thread.interrupted();
                while (!reporter.isTerminated()) {
                    try {
                        reporter.awaitTermination(1, TimeUnit.SECONDS);
                    } catch (InterruptedException ex) {
                        interrupted = true;
                    }
                }
                if (interrupted) Thread.currentThread().interrupt();
            }
        }
    }

    @Override
    protected void done() {
        try {
            model.addSynthonSpace(new LoadedSynthonSpace(get(), name, threads));
            setProgress(100);
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();
            reportFailure(ex);
        } catch (ExecutionException | CancellationException ex) {
            reportFailure(ex instanceof ExecutionException ? ex.getCause() : ex);
        }
    }

    private void reportFailure(Throwable error) {
        failure = error.toString();
        System.err.println("Cannot load " + filepath + ": " + failure);
        firePropertyChange("loadError", null, failure);
        if (!GraphicsEnvironment.isHeadless())
            JOptionPane.showMessageDialog(
                    null,
                    "Cannot load " + name + "\n" + failure,
                    "Space loading failed",
                    JOptionPane.ERROR_MESSAGE);
    }
}
