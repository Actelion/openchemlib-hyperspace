package com.idorsia.research.chem.hyperspace.importer;

import java.util.ArrayList;
import java.util.Set;

public class ImporterSynthonSet {

    public final String filePath;
    public final ArrayList<ImporterSynthon> synthons;

    public Set<Integer> getConnectors() {
        return synthons.get(0).getUsedConnectors();
    }

    public ImporterSynthonSet(String filePath, ArrayList<ImporterSynthon> synthons) {
        this.filePath = filePath;
        this.synthons = synthons;
    }

}
