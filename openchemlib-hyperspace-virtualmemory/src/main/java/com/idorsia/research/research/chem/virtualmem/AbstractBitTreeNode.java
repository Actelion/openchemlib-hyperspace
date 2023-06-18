package com.idorsia.research.research.chem.virtualmem;

import java.util.BitSet;
import java.util.Iterator;

public abstract class AbstractBitTreeNode {

    protected int bits;
    protected int splitBit;
    //protected AbstractBitTreeNode left;
    //protected AbstractBitTreeNode right;

    public abstract Iterator<BitSet> getBitSetIterator();
    public abstract void initNodeData(Iterator<BitSet> data, String dbPath, int numBitSetsInMemory, int numBitSetsInLeaf);

    public abstract boolean isLeaf();
    public abstract AbstractBitTreeNode getLeftChild();
    public abstract AbstractBitTreeNode getRightChild();

}
