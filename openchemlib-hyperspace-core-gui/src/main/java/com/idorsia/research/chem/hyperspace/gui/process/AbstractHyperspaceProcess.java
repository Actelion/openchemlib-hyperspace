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

    private ProcessStatus status       = ProcessStatus.WAITING;
    private String        status_msg   = "";

    public ProcessStatus getProcessStatus() {
        return this.status;
    }

    public String getProcessStatusMessage() {
        return this.status_msg;
    }

    //public abstract List<StereoMolecule> getQueryStructures();


    private List<HyperspaceProcessListener> listeners = new ArrayList<>();

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


    private boolean waitUntilDoneOrFailed_Flag_done        = false;
    private boolean waitUntilDoneOrFailed_Flag_interrupted = false;

    /**
     * Returns true if the process reached DONE or FAILED, returns
     * false if the timeout was hit.
     *
     * @param timeout_ms
     * @return
     */
    public boolean waitUntilDoneOrFailed(int timeout_ms) throws InterruptedException {
        long ts_a = System.currentTimeMillis();

        //final AbstractHyperspaceProcess thisObject = this;
        this.addSearchProviderListener(new HyperspaceProcessListener() {
            @Override
            public void processStatusChanged() {
                if(getProcessStatus() == ProcessStatus.DONE || getProcessStatus() == ProcessStatus.FAILED) {
                    waitUntilDoneOrFailed_Flag_done = true;
                    return;
                }
                else {
                    //System.out.println("[waitUntilDoneOrFailed] :: continue waiting");
                }
            }
        });

        while( System.currentTimeMillis() < (ts_a + timeout_ms) ) {
            Thread.sleep(20);
            if( waitUntilDoneOrFailed_Flag_done || waitUntilDoneOrFailed_Flag_interrupted ) {
                break;
            }
        }

        if( System.currentTimeMillis() >= (ts_a + timeout_ms) ) {
            waitUntilDoneOrFailed_Flag_interrupted = true;
        }

        return waitUntilDoneOrFailed_Flag_done;
    }

    public static interface HasProgress {
        public double getProgress();
    }

}
