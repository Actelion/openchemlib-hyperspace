package com.idorsia.research.chem.hyperspace.importer;

import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class ConcreteMultiStepSynthonReaction implements MultiStepSynthonReaction {

    private ImporterSynthonReaction originalSynthonReaction;
    private ImporterSynthonReaction connectorRemappedSynthonReaction;

    public ConcreteMultiStepSynthonReaction(ImporterSynthonReaction originalSynthonReaction, ConnectorRemap remap) {
        if(remap==null) {remap = new ConnectorRemap(new HashMap<>());}
        this.originalSynthonReaction = originalSynthonReaction;
        this.connectorRemappedSynthonReaction = ConnectorRemap.remap(remap,originalSynthonReaction);
    }

    @Override
    public boolean isConcrete() {
        return true;
    }

    @Override
    public Set<Integer> getOutConnectors() {
        return connectorRemappedSynthonReaction.getOutConnectors();
    }

    @Override
    public Set<Integer> getAllUsedConnectors() {
        return connectorRemappedSynthonReaction.getAllUsedConnectors();
    }

    @Override
    public List<ImporterSynthonSet> getAllSynthonSets() {
        return this.connectorRemappedSynthonReaction.getSynthonSets();
    }

    @Override
    public void remapConnectors(ConnectorRemap remap) {
        this.connectorRemappedSynthonReaction = ConnectorRemap.remap(remap,this.connectorRemappedSynthonReaction);
    }
}
