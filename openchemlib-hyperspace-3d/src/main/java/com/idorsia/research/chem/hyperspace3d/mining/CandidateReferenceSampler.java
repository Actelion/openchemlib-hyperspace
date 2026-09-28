package com.idorsia.research.chem.hyperspace3d.mining;

import com.idorsia.research.chem.hyperspace3d.index.MoleculeVectorReference;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Bounded-memory bottom-hash sample of stable source vector references. */
public final class CandidateReferenceSampler {
    private final long seed;
    private final long[] priorities;
    private final long[] references;
    private int size;

    public CandidateReferenceSampler(int capacity, long seed) {
        if (capacity < 1) throw new IllegalArgumentException("capacity must be positive");
        this.seed = seed;
        this.priorities = new long[capacity];
        this.references = new long[capacity];
    }

    public void offer(int shard, long localRow) {
        if (shard < 0 || shard > 0xffff || localRow < 0 || localRow > 0x0000_ffff_ffff_ffffL) {
            throw new IllegalArgumentException("source reference cannot be packed");
        }
        long reference = ((long) shard << 48) | localRow;
        long priority = PredictedRankMiner.mix(seed, shard, localRow);
        if (size < priorities.length) {
            priorities[size] = priority; references[size] = reference; siftUp(size++);
        } else if (unsigned(priority, priorities[0]) < 0
                || (priority == priorities[0] && reference < references[0])) {
            priorities[0] = priority; references[0] = reference; siftDown(0);
        }
    }

    public List<SampledReference> results() {
        List<SampledReference> result = new ArrayList<>(size);
        for (int index = 0; index < size; index++) {
            long packed = references[index];
            result.add(new SampledReference(new MoleculeVectorReference(
                    (int) (packed >>> 48), packed & 0x0000_ffff_ffff_ffffL), priorities[index]));
        }
        result.sort((left, right) -> {
            int value = unsigned(left.priority, right.priority);
            if (value != 0) return value;
            value = Integer.compare(left.reference.shardIndex(), right.reference.shardIndex());
            return value != 0 ? value : Long.compare(left.reference.localRow(), right.reference.localRow());
        });
        return result;
    }

    private void siftUp(int index) {
        while (index > 0) {
            int parent = (index - 1) >>> 1;
            if (worseOrEqual(parent, index)) return;
            swap(parent, index); index = parent;
        }
    }
    private void siftDown(int index) {
        while (true) {
            int left = index * 2 + 1;
            if (left >= size) return;
            int worst = left;
            int right = left + 1;
            if (right < size && !worseOrEqual(left, right)) worst = right;
            if (worseOrEqual(index, worst)) return;
            swap(index, worst); index = worst;
        }
    }
    private boolean worseOrEqual(int left, int right) {
        int value = unsigned(priorities[left], priorities[right]);
        return value > 0 || (value == 0 && references[left] >= references[right]);
    }
    private void swap(int left, int right) {
        long value = priorities[left]; priorities[left] = priorities[right]; priorities[right] = value;
        value = references[left]; references[left] = references[right]; references[right] = value;
    }
    private static int unsigned(long left, long right) {
        return Long.compare(left ^ Long.MIN_VALUE, right ^ Long.MIN_VALUE);
    }

    public record SampledReference(MoleculeVectorReference reference, long priority) {}
}
