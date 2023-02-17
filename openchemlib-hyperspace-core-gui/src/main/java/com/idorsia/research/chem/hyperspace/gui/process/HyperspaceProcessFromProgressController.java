package com.idorsia.research.chem.hyperspace.gui.process;

import com.actelion.research.calc.ProgressController;

public class HyperspaceProcessFromProgressController extends AbstractHyperspaceProcess implements ProgressController , AbstractHyperspaceProcess.HasProgress {

    public HyperspaceProcessFromProgressController() {

    }

    @Override
    public String getName() {
        return null;
    }

    private int target = 100;
    private int progress = 0;

    @Override
    public void startProgress(String s, int i, int i1) {
        this.target = i1;
        this.setProgress(i);
        this.setProcessStatus(ProcessStatus.COMPUTING);
        this.setProcessStatusMessage(s);
    }

    @Override
    public void updateProgress(int i) {
        this.setProcessStatus(ProcessStatus.COMPUTING);
        this.setProgress(i);
    }

    public void setProgress(int progress) {
        if(progress >= target) {
            this.progress = progress;
            this.setProcessStatus(ProcessStatus.DONE);
            return;
        }
        this.progress = progress;

    }

    @Override
    public double getProgress() {
        return (1.0*progress) / (1.0*target);
    }

    @Override
    public void updateProgress(int i, String s) {
        this.setProcessStatusMessage(s);
        this.setProgress(i);
    }

    @Override
    public void stopProgress() {

    }

    @Override
    public void showErrorMessage(String s) {

    }

    @Override
    public boolean threadMustDie() {
        return false;
    }
}
