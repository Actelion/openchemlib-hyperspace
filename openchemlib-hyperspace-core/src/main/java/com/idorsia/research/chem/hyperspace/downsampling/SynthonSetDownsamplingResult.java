package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

/**
 * Result container for a single reaction fragment set.
 */
public class SynthonSetDownsamplingResult implements Serializable {

    private final SynthonSpace.FragType fragType;
    private final List<SynthonSpace.FragId> representatives;
    private final List<ClusterInfo> clusterInfos;
    private final int originalSize;

    public SynthonSetDownsamplingResult(SynthonSpace.FragType fragType,
                                        List<SynthonSpace.FragId> representatives,
                                        List<ClusterInfo> clusterInfos,
                                        int originalSize) {
        this.fragType = Objects.requireNonNull(fragType, "fragType");
        this.representatives = Collections.unmodifiableList(new ArrayList<>(representatives));
        this.clusterInfos = Collections.unmodifiableList(new ArrayList<>(clusterInfos));
        this.originalSize = originalSize;
    }

    public SynthonSpace.FragType getFragType() {
        return fragType;
    }

    public List<SynthonSpace.FragId> getRepresentatives() {
        return representatives;
    }

    public List<ClusterInfo> getClusterInfos() {
        return clusterInfos;
    }

    public int getOriginalSize() {
        return originalSize;
    }

    public int getRetainedSize() {
        return representatives.size();
    }

    public int getDiscardedSize() {
        return Math.max(0, originalSize - getRetainedSize());
    }

    public static final class ClusterInfo implements Serializable {
        private final SynthonSpace.FragId representative;
        private final int members;
        private final double minSimilarityWithinCluster;

        public ClusterInfo(SynthonSpace.FragId representative, int members, double minSimilarityWithinCluster) {
            this.representative = Objects.requireNonNull(representative, "representative");
            this.members = members;
            this.minSimilarityWithinCluster = minSimilarityWithinCluster;
        }

        public SynthonSpace.FragId getRepresentative() {
            return representative;
        }

        public int getMembers() {
            return members;
        }

        public double getMinSimilarityWithinCluster() {
            return minSimilarityWithinCluster;
        }
    }
}
