package com.idorsia.research.chem.hyperspace.tree;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.tree.TreeNode;
import java.util.Enumeration;
import java.util.List;
import java.util.stream.Collectors;

public class SynthonRxnTreeNode implements TreeNode {

    private SynthonSpaceTreeNode parent;
    private SynthonSpace space;
    private SynthonSpace.RxnId rxn;

    public SynthonRxnTreeNode(SynthonSpaceTreeNode parent, SynthonSpace space, String rxn) {
        this.parent = parent;
        this.space = space;
        this.rxn = space.getRxn(rxn);
    }

    public SynthonSpace.RxnId getRxnId() {
        return this.rxn;
    }

    @Override
    public TreeNode getChildAt(int childIndex) {
        List<Integer> fids = space.getFragTypes(rxn.id).keySet().stream().sorted().collect(Collectors.toList());
        return new SynthonSetTreeNode(this,space,this.space.getFragTypes(this.rxn.id).get(fids.get(childIndex)));
    }

    @Override
    public int getChildCount() {
        return this.space.getFragTypes(this.rxn.id).size();
    }

    @Override
    public TreeNode getParent() {
        return this.parent;
    }

    @Override
    public int getIndex(TreeNode node) {
        if(node instanceof SynthonSetTreeNode) {
            List<Integer> fids = space.getFragTypes(rxn.id).keySet().stream().sorted().collect(Collectors.toList());
            return fids.indexOf( ((SynthonSetTreeNode)node).getFragType().frag );
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
        return null;
    }

    @Override
    public String toString() {
        return this.rxn.id;
    }

}
