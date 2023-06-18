package com.idorsia.research.chem.hyperspace.tree;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.tree.TreeNode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;

public class SynthonTreeNode implements TreeNode {

    private SynthonSetTreeNode parentNode;
    private SynthonSpace.FragId synthon;

    public SynthonTreeNode(SynthonSetTreeNode parentNode, SynthonSpace.FragId synthon) {
        this.parentNode = parentNode;
        this.synthon = synthon;
    }

    public SynthonSpace.FragId getSynthon() {
        return this.synthon;
    }

    @Override
    public TreeNode getChildAt(int childIndex) {
        return null;
    }

    @Override
    public int getChildCount() {
        return 0;
    }

    @Override
    public TreeNode getParent() {
        return this.parentNode;
    }

    @Override
    public int getIndex(TreeNode node) {
        return 0;
    }

    @Override
    public boolean getAllowsChildren() {
        return false;
    }

    @Override
    public boolean isLeaf() {
        return true;
    }

    @Override
    public Enumeration<? extends TreeNode> children() {
        return Collections.enumeration(new ArrayList());
    }

    @Override
    public String toString() {
        return this.synthon.fragment_id;
    }
}
