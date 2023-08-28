package com.idorsia.research.chem.hyperspace.gui2.task;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.*;

public class LongRunningTask extends SwingWorker<Void, Integer> {


    private String name;


    public LongRunningTask(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    @Override
    protected Void doInBackground() throws Exception {
        // Perform the long-running task here
        int totalProgress = 100;
        for (int progress = 0; progress <= totalProgress; progress++) {
            // Update the progress
            //setProgress(progress);
            publish(progress);

            // Simulate work by sleeping for a while
            Thread.sleep(500);

            // Check if the task is cancelled
            if (isCancelled()) {
                break;
            }
        }
        return null;
    }

    @Override
    protected void process(java.util.List<Integer> chunks) {
        // Update the GUI with the progress values received
        int latestProgress = chunks.get(chunks.size() - 1);
        this.setProgress(latestProgress);
        //outputTextArea.append("Progress: " + latestProgress + "%\n");
    }

    @Override
    protected void done() {
        // Executed on the Event Dispatch Thread (EDT) when the task is complete
        //outputTextArea.append("Task completed!\n");
    }
}
