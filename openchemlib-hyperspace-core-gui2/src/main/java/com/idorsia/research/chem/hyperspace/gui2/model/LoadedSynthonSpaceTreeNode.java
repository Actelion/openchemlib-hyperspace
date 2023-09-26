package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.tree.SynthonSpaceTreeNode;

public class LoadedSynthonSpaceTreeNode extends SynthonSpaceTreeNode {

    private LoadedSynthonSpace loadedSynthonSpace;

    public LoadedSynthonSpaceTreeNode(LoadedSynthonSpace space) {
        super(space.getSpace());
        this.loadedSynthonSpace = space;
    }

    @Override
    public String toString() {
        return this.loadedSynthonSpace.getName();
    }
}
