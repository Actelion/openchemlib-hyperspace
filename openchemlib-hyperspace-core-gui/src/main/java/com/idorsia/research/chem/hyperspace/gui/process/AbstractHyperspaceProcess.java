package com.idorsia.research.chem.hyperspace.gui.process;

import com.actelion.research.chem.StereoMolecule;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

/**
 * Idorsia Pharmaceuticals Ltd. 2021
 * Thomas Liphardt
 *
 * Hyperspace-GUI
 */




public abstract class AbstractHyperspaceProcess {

    public static enum ProcessStatus { WAITING, COMPUTING , DONE , FAILED }

    public abstract String getName();

    private volatile ProcessStatus status       = ProcessStatus.WAITING;
    private volatile String        status_msg   = "";

    public ProcessStatus getProcessStatus() {
        return this.status;
    }

    public String getProcessStatusMessage() {
        return this.status_msg;
    }

    //public abstract List<StereoMolecule> getQueryStructures();


    private List<HyperspaceProcessListener> listeners = new java.util.concurrent.CopyOnWriteArrayList<>();

    public static interface HyperspaceProcessListener {
        public void processStatusChanged();
    }

    protected void setProcessStatus(ProcessStatus ps) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                status = ps;
                fireProcessStatusChanged();
            }
        });
    }

    protected void setProcessStatusMessage(String status_message) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                status_msg = status_message;
                fireProcessStatusChanged();
            }
        });
    }

    protected void fireProcessStatusChanged() {
        for(HyperspaceProcessListener li : listeners) {
            li.processStatusChanged();
        }
    }

    public void addSearchProviderListener(HyperspaceProcessListener li) {
        this.listeners.add(li);
    }
    public void removeSearchProviderListener(HyperspaceProcessListener li) {
        this.listeners.remove(li);
    }


    /** Waits for completion, including a process that completed before this call. */
    public boolean waitUntilDoneOrFailed(int timeout_ms) throws InterruptedException {
        if (SwingUtilities.isEventDispatchThread()) throw new IllegalStateException("Do not wait for processes on the Swing event thread");
        long deadline = System.nanoTime() + java.util.concurrent.TimeUnit.MILLISECONDS.toNanos(timeout_ms);
        do {
            if (status == ProcessStatus.DONE || status == ProcessStatus.FAILED) return true;
            Thread.sleep(20);
        } while (System.nanoTime() < deadline);
        return status == ProcessStatus.DONE || status == ProcessStatus.FAILED;
    }

    public static interface HasProgress {
        public double getProgress();
    }

}
