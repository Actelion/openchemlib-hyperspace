package com.idorsia.research.research.chem.virtualmem;

import org.rocksdb.Options;
import org.rocksdb.RocksDB;
import org.rocksdb.RocksDBException;
import org.rocksdb.RocksIterator;

import java.util.BitSet;
import java.util.Iterator;
import java.util.Random;

public class RocksDBBenchmark {

    public static void main(String[] args) throws RocksDBException {
        // Open the RocksDB database
        RocksDB.loadLibrary();
        Options options = new Options().setCreateIfMissing(true);
        RocksDB db = RocksDB.open(options, "C:\\Temp\\rocksdb\\test_a.xdb");

        int bitSetLength = 512; // Length of each BitSet
        int numBitSets = 10000000; // Number of BitSet objects to store
        int batchSize = 100; // Number of BitSets to generate in each batch

        // Create an iterator to generate the BitSet objects
        BitSetIterator bitSetIterator = new BitSetIterator(bitSetLength, numBitSets);

        // Store the BitSet objects in RocksDB
        System.out.println("Store bitsets:\n");
        while (bitSetIterator.hasNext()) {
            BitSet bitSet = bitSetIterator.next();
            byte[] key = Integer.toString(bitSetIterator.getCurrentIndex()).getBytes();
            byte[] value = bitSet.toByteArray();
            db.put(key, value);
        }

        System.out.println("Iterate over bitsets:");
        RocksIterator iterator = db.newIterator();
        for (iterator.seekToFirst(); iterator.isValid(); iterator.next()) {
            byte[] key = iterator.key();
            byte[] value = iterator.value();
            BitSet bitSet = BitSet.valueOf(value);
            //processBitSet(bitSet);
        }

        // Close the RocksDB database
        db.close();

        // Close the RocksDB database
        db.close();
    }

    // Iterator to generate random BitSets
    private static class BitSetIterator implements Iterator<BitSet> {
        private final int bitSetLength;
        private int currentIndex;
        private int remainingBitSets;

        private Random r = new Random(123);
        public BitSetIterator(int bitSetLength, int numBitSets) {
            this.bitSetLength = bitSetLength;
            this.currentIndex = 0;
            this.remainingBitSets = numBitSets;
        }

        @Override
        public boolean hasNext() {
            return remainingBitSets > 0;
        }

        @Override
        public BitSet next() {
            BitSet bitSet = generateRandomBitSet();
            return bitSet;
        }

        public int getCurrentIndex() {
            return currentIndex;
        }

        private BitSet generateRandomBitSet() {
            BitSet bitSet = new BitSet(bitSetLength);
            for (int i = 0; i < bitSetLength; i++) {
                if (r.nextDouble() < 0.5) {
                    bitSet.set(i);
                }
            }
            return bitSet;
        }
    }

}
