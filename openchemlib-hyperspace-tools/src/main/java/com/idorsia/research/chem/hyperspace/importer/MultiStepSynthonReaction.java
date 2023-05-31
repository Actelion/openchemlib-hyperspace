package com.idorsia.research.chem.hyperspace.importer;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Set;

public interface MultiStepSynthonReaction {

    public boolean isConcrete();
    public Set<Integer> getOutConnectors();
    public Set<Integer> getAllUsedConnectors();
    public List<ImporterSynthonSet> getAllSynthonSets();

    /**
     * The MultiStepSynthonReaction supports remapping of connectors.
     * This function must always update the current state of remapping, i.e.
     * if we currently have some remapping in place, then this additional
     * remap is based on the existing remap and should update the
     * remap accordingly.
     *
     * @param rempa
     */
    public void remapConnectors(ConnectorRemap rempa);
}
