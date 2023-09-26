package com.idorsia.research.chem.hyperspace.gui2.model;

import com.idorsia.research.chem.hyperspace.tree.SynthonSpaceTreeNode;

import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreeNode;
import java.util.Enumeration;

public class SynthonSpacesTreeModel implements LeetHyperspaceModel.LeetHyperspaceModelListener {

    private LeetHyperspaceModel model;
    private DefaultTreeModel    treeModel;

    private SynthonSpaceListRootNode root;

    private class SynthonSpaceListRootNode implements TreeNode {
        @Override
        public TreeNode getChildAt(int childIndex) {
            return new LoadedSynthonSpaceTreeNode(model.getSynthonSpaces().get(childIndex));
        }

        @Override
        public int getChildCount() {
            return model.getSynthonSpaces().size();
        }

        @Override
        public TreeNode getParent() {
            return null;
        }

        @Override
        public int getIndex(TreeNode node) {
            if(node instanceof SynthonSpaceTreeNode) {
                return model.getSynthonSpaces().indexOf(((SynthonSpaceTreeNode)node).getSynthonSpace());
            }
            return -1;
        }

        @Override
        public boolean getAllowsChildren() {
            return true;
        }

        @Override
        public boolean isLeaf() {
            return model.getSynthonSpaces().isEmpty();
        }

        @Override
        public Enumeration<? extends TreeNode> children() {
            return null;
        }
    }

    public SynthonSpacesTreeModel(LeetHyperspaceModel model) {
        this.model = model;
        this.root = new SynthonSpaceListRootNode();
        this.treeModel = new DefaultTreeModel(this.root);
        model.addListener(this);
    }

    public TreeModel getTreeModel() {
        return this.treeModel;
    }

    @Override
    public void synthonSpacesChanged() {
        treeModel.nodeStructureChanged(root);
    }
}
