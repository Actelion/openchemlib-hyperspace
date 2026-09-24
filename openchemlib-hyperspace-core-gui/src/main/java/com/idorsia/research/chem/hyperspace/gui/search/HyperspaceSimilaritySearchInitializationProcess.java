package com.idorsia.research.chem.hyperspace.gui.search;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSimilaritySpace3;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import org.apache.commons.io.input.CountingInputStream;
import java.io.*;
import java.nio.file.*;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class HyperspaceSimilaritySearchInitializationProcess extends AbstractHyperspaceProcess implements AbstractHyperspaceProcess.HasProgress {
    private final String file;
    private volatile long fileSize;
    private volatile CountingInputStream counting;
    private volatile SynthonSimilaritySpace3 space;
    public HyperspaceSimilaritySearchInitializationProcess(AbstractSearchProvider provider, String file) { this.file = file; }
    public void startInitializationAsync() {
        setProcessStatus(ProcessStatus.COMPUTING);
        new Thread(() -> {
            try {
                fileSize = Files.size(Path.of(file));
                counting = new CountingInputStream(new BufferedInputStream(Files.newInputStream(Path.of(file))));
                try (CountingInputStream bytes = counting; ObjectInputStream in = new ObjectInputStream(new GZIPInputStream(bytes))) {
                    space = (SynthonSimilaritySpace3) in.readObject();
                    space.initAfterJavaDeserialization();
                }
                setProcessStatus(ProcessStatus.DONE);
            } catch (Exception | LinkageError e) {
                setProcessStatusMessage(file + ": " + e.toString());
                setProcessStatus(ProcessStatus.FAILED);
            } catch (OutOfMemoryError e) {
                setProcessStatusMessage("Insufficient Java heap while loading " + file + ". Restart the GUI with a larger -Xmx setting.");
                setProcessStatus(ProcessStatus.FAILED);
            }
        }, "hyperspace-load").start();
    }
    public double getStatusOfInitialization() {
        if (getProcessStatus() == ProcessStatus.DONE) return 1.0;
        return counting == null || fileSize <= 0 ? 0.0 : Math.min(0.95, 0.95 * counting.getByteCount() / fileSize);
    }
    @Override public String getName() { return "Load similarity space: " + Path.of(file).getFileName(); }
    public SynthonSimilaritySpace3 getSpace() { return space; }
    public List<StereoMolecule> getQueryStructures() { return List.of(); }
    @Override public double getProgress() { return getStatusOfInitialization(); }
}
