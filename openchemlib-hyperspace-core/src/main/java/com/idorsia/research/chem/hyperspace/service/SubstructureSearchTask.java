package com.idorsia.research.chem.hyperspace.service;

import com.actelion.research.chem.StereoMolecule;

public class SubstructureSearchTask {


    private String queryIDCode = "";
    private boolean screenBuildingBlocks = true;
    private String spaceName = "";

    private SubstructureSearchTaskConfiguration configuration = new SubstructureSearchTaskConfiguration();

    public SubstructureSearchTask() {}

    public SubstructureSearchTask(String queryIDCode, boolean screenBuildingBlocks, String spaceName, SubstructureSearchTaskConfiguration configuration) {
        this.queryIDCode = queryIDCode;
        this.screenBuildingBlocks = screenBuildingBlocks;
        this.spaceName = spaceName;
        this.configuration = configuration;
    }

    public String getQueryIDCode() {
        return queryIDCode;
    }

    public void setQueryIDCode(String queryIDCode) {
        this.queryIDCode = queryIDCode;
    }

    public boolean isScreenBuildingBlocks() {
        return screenBuildingBlocks;
    }

    public void setScreenBuildingBlocks(boolean screenBuildingBlocks) {
        this.screenBuildingBlocks = screenBuildingBlocks;
    }

    public String getSpaceName() {
        return spaceName;
    }

    public void setSpaceName(String spaceName) {
        this.spaceName = spaceName;
    }

    public SubstructureSearchTaskConfiguration getConfiguration() {
        return configuration;
    }

    public void setConfiguration(SubstructureSearchTaskConfiguration configuration) {
        this.configuration = configuration;
    }
}
