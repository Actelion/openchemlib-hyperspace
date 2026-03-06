package com.idorsia.research.synthonbrowser;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import javax.swing.tree.DefaultMutableTreeNode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

final class SynthonTreeBuilder {

    private static final DecimalFormat MILLION_FORMAT = new DecimalFormat("#,##0.0");

    enum NodeType {
        ROOT,
        REACTION,
        SET
    }

    record NodeData(
            NodeType type,
            String reactionId,
            Integer setIndex,
            int synthonCount,
            String label
    ) {
        @Override
        public String toString() {
            return label;
        }
    }

    private SynthonTreeBuilder() {}

    static DefaultMutableTreeNode buildTree(RawSynthonSpace space, boolean useDownsampledSets) {
        int totalReactions = space.getReactions().size();
        int totalSynthons = 0;
        for (RawSynthonSpace.ReactionData rd : space.getReactions().values()) {
            Map<Integer, List<RawSynthon>> sets = activeSets(rd, useDownsampledSets);
            totalSynthons += countSynthons(sets);
        }

        String rootLabel = String.format(
                "%s (v%s) - %d reactions, %d synthons%s",
                safe(space.getName()),
                safe(space.getVersion()),
                totalReactions,
                totalSynthons,
                useDownsampledSets ? " [downsampled]" : ""
        );
        DefaultMutableTreeNode root = new DefaultMutableTreeNode(
                new NodeData(NodeType.ROOT, null, null, totalSynthons, rootLabel)
        );

        List<String> reactionIds = new ArrayList<>(space.getReactions().keySet());
        reactionIds.sort(Comparator.naturalOrder());
        for (String reactionId : reactionIds) {
            RawSynthonSpace.ReactionData rd = space.getReactions().get(reactionId);
            if (rd == null) {
                continue;
            }
            Map<Integer, List<RawSynthon>> sets = activeSets(rd, useDownsampledSets);
            int reactionSynthons = countSynthons(sets);
            long estimatedCompounds = estimateCompoundCount(sets);
            String reactionLabel = String.format(
                    "%s (%d sets, %d synthons, %sM compounds)",
                    reactionId,
                    sets.size(),
                    reactionSynthons,
                    formatMillions(estimatedCompounds)
            );
            DefaultMutableTreeNode reactionNode = new DefaultMutableTreeNode(
                    new NodeData(NodeType.REACTION, reactionId, null, reactionSynthons, reactionLabel)
            );
            root.add(reactionNode);

            List<Integer> setIndexes = new ArrayList<>(sets.keySet());
            setIndexes.sort(Comparator.naturalOrder());
            for (Integer setIndex : setIndexes) {
                List<RawSynthon> synthons = sets.get(setIndex);
                int setSize = synthons == null ? 0 : synthons.size();
                String setLabel = String.format("Set %d (%d)", setIndex, setSize);
                DefaultMutableTreeNode setNode = new DefaultMutableTreeNode(
                        new NodeData(NodeType.SET, reactionId, setIndex, setSize, setLabel)
                );
                reactionNode.add(setNode);
            }
        }
        return root;
    }

    static Map<Integer, List<RawSynthon>> activeSets(
            RawSynthonSpace.ReactionData reactionData,
            boolean useDownsampledSets
    ) {
        if (useDownsampledSets && reactionData.getRawDownsampledSets() != null && !reactionData.getRawDownsampledSets().isEmpty()) {
            return reactionData.getRawDownsampledSets();
        }
        return reactionData.getRawFragmentSets();
    }

    static int countSynthons(Map<Integer, List<RawSynthon>> sets) {
        int count = 0;
        for (List<RawSynthon> synthons : sets.values()) {
            if (synthons != null) {
                count += synthons.size();
            }
        }
        return count;
    }

    static long estimateCompoundCount(Map<Integer, List<RawSynthon>> sets) {
        long product = 1L;
        boolean hasAnySet = false;
        for (List<RawSynthon> synthons : sets.values()) {
            if (synthons == null || synthons.isEmpty()) {
                return 0L;
            }
            hasAnySet = true;
            int size = synthons.size();
            if (product > Long.MAX_VALUE / size) {
                return Long.MAX_VALUE;
            }
            product *= size;
        }
        return hasAnySet ? product : 0L;
    }

    static String formatMillions(long count) {
        return MILLION_FORMAT.format(count / 1_000_000.0);
    }

    private static String safe(String value) {
        return value == null ? "N/A" : value;
    }
}
