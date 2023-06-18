package com.idorsia.research.chem.hyperspace.tree;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.tree.TreeNode;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.stream.Collectors;

public class SynthonSpaceTreeNode implements TreeNode {

    private SynthonSpace space;

    public SynthonSpaceTreeNode(SynthonSpace space) {
        this.space = space;
    }

    public SynthonSpace getSynthonSpace() {
        return this.space;
    }

    @Override
    public TreeNode getChildAt(int childIndex) {
        return new SynthonRxnTreeNode( this, this.space, this.space.getRxnIds().get(childIndex) );
    }

    @Override
    public int getChildCount() {
        return this.space.getRxnIds().size();
    }

    @Override
    public TreeNode getParent() {
        return null;
    }

    @Override
    public int getIndex(TreeNode node) {
        if(node instanceof SynthonRxnTreeNode) {
            SynthonSpace.RxnId rxn_id = ((SynthonRxnTreeNode)node).getRxnId();
            return this.space.getRxnIds().indexOf(rxn_id);
        }

        return -1;
    }

    @Override
    public boolean getAllowsChildren() {
        return true;
    }

    @Override
    public boolean isLeaf() {
        return false;
    }

    @Override
    public Enumeration<? extends TreeNode> children() {
        List<SynthonRxnTreeNode> nodes = this.space.getRxnIds().stream().map(xi -> new SynthonRxnTreeNode(this,this.space,xi)).collect(Collectors.toList());
        return Collections.enumeration(nodes);
    }

    @Override
    public String toString() {
        return "SynthonSpace";
    }

}
