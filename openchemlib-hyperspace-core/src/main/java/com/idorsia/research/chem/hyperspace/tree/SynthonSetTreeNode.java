package com.idorsia.research.chem.hyperspace.tree;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.tree.TreeNode;
import java.util.Collections;
import java.util.Enumeration;
import java.util.List;
import java.util.stream.Collectors;

public class SynthonSetTreeNode implements TreeNode {

    private SynthonRxnTreeNode parentNode;
    private SynthonSpace space;
    private SynthonSpace.FragType fragType;

    public SynthonSetTreeNode(SynthonRxnTreeNode parentNode, SynthonSpace space, SynthonSpace.FragType fragType) {
        this.parentNode = parentNode;
        this.space = space;
        this.fragType = fragType;
    }

    public SynthonSpace.FragType getFragType() {
        return this.fragType;
    }

    @Override
    public TreeNode getChildAt(int childIndex) {
        return new SynthonTreeNode(this, space.getSynthonSet(fragType.rxn_id,fragType.frag).get(childIndex) );
    }

    @Override
    public int getChildCount() {
        return space.getSynthonSet(fragType.rxn_id,fragType.frag).size();
    }

    @Override
    public TreeNode getParent() {
        return parentNode;
    }

    @Override
    public int getIndex(TreeNode node) {
        if(node instanceof SynthonTreeNode) {
            SynthonSpace.FragId synthon = ((SynthonTreeNode) node).getSynthon();
            return space.getSynthonSet(fragType.rxn_id, fragType.frag).indexOf( synthon );
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

    private SynthonSetTreeNode getThis() {return this;}

    @Override
    public Enumeration<? extends TreeNode> children() {
        List<SynthonTreeNode> nodes = this.space.getSynthonSet(this.fragType.rxn_id,this.fragType.frag).stream().map(xi -> new SynthonTreeNode(getThis(),xi)).collect(Collectors.toList());
        return Collections.enumeration(nodes);
    }

    @Override
    public String toString() {
        try {
            String sa = "F:" + this.fragType.frag + "(n=" + space.getSynthonSet(fragType.rxn_id, fragType.frag).size() + ")";
            return sa;
        }
        catch(Exception ex) {ex.printStackTrace();}
        return "F:error";
    }
}
