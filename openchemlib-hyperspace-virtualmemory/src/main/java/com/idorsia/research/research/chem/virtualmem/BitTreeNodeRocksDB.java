package com.idorsia.research.research.chem.virtualmem;

import org.rocksdb.RocksDB;
import org.rocksdb.RocksDBException;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

public class BitTreeNodeRocksDB extends AbstractBitTreeNode {

    private RocksDB db;

    private boolean isLeaf;

    /**
     * Splitting strategy parameters..
     */
    private int    splitStrategy_MaxTries    = 20;
    private double splitStrategy_acceptRatio = 0.35;

    public BitTreeNodeRocksDB(RocksDB db) {
        this.db = db;
    }

    @Override
    public boolean isLeaf() {
        return isLeaf;
    }

    @Override
    public AbstractBitTreeNode getLeftChild() {
        return null;
    }

    @Override
    public AbstractBitTreeNode getRightChild() {
        return null;
    }

    public void initNodeData(Iterator<BitSet> data, String dbPath, int numBitSetsInMemory, int numBitSetsInLeaf) {
        try {
            RocksDB.loadLibrary();
            db = RocksDB.open(dbPath);
        } catch (RocksDBException e) {
            e.printStackTrace();
            return;
        }

        // Temp storage in case that this will be a leaf
        //List<BitSet> bitSetList = new ArrayList<>();
        //data.forEachRemaining(bitSetList::add);
        // Decide if this node should be a leaf

        /**
        if (bitSetList.size() <= numBitSetsInLeaf) {
            isLeaf = true;
            bitSets = bitSetList;
            return;
        }

        isLeaf = false;
        // Determine the splitting bit
        this.splitBit = findSplittingBit.apply(bitSetList.iterator());

        // Collect BitSets for left and right children
        List<BitSet> leftBitSets = new ArrayList<>();
        List<BitSet> rightBitSets = new ArrayList<>();

        for (BitSet bitSet : bitSetList) {
            if (bitSet.get(splitBit)) {
                rightBitSets.add(bitSet);
            } else {
                leftBitSets.add(bitSet);
            }
        }

        // Initialize children nodes
        this.left = new RocksDBBitTreeNode();
        this.left.initNodeData(leftBitSets.iterator(), dbPath + "0", numBitSetsInMemory, numBitSetsInLeaf, findSplittingBit);

        this.right = new RocksDBBitTreeNode();
        this.right.initNodeData(rightBitSets.iterator(), dbPath + "1", numBitSetsInMemory, numBitSetsInLeaf, findSplittingBit);

         **/
    }

    public Iterator<BitSet> getBitSetIterator() {
        if (isLeaf) {
            return null;//bitSets.iterator();
        } else {
            Iterator<BitSet> leftIterator = getLeftChild().getBitSetIterator();
            Iterator<BitSet> rightIterator = getRightChild().getBitSetIterator();
            return new Iterator<BitSet>() {
                @Override
                public boolean hasNext() {
                    return leftIterator.hasNext() || rightIterator.hasNext();
                }

                @Override
                public BitSet next() {
                    if (leftIterator.hasNext()) {
                        return leftIterator.next();
                    } else {
                        return rightIterator.next();
                    }
                }
            };
        }
    }



}
