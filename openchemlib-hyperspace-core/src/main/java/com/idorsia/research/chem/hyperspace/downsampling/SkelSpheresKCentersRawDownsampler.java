package com.idorsia.research.chem.hyperspace.downsampling;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSet;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Raw-synthon variant of the SkelSpheres k-centers downsampler.
 */
public class SkelSpheresKCentersRawDownsampler implements RawSynthonDownsampler {

    private final RawSynthonDescriptorCache<int[]> descriptorCache;
    private final ThreadLocal<IDCodeParser> parsers;
    private final ThreadLocal<DescriptorHandlerBinarySkelSpheres> descriptorHandlers;

    public SkelSpheresKCentersRawDownsampler() {
        this.parsers = ThreadLocal.withInitial(IDCodeParser::new);
        this.descriptorHandlers = ThreadLocal.withInitial(DescriptorHandlerBinarySkelSpheres::new);
        this.descriptorCache = new RawSynthonDescriptorCache<>(this::createDescriptor);
    }

    @Override
    public String getName() {
        return "SkelSpheresKCentersRaw";
    }

    @Override
    public SynthonSetDownsamplingResult downsample(RawSynthonSet synthonSet,
                                                   SynthonDownsamplingRequest request) {
        try {
            List<RawSynthon> synthons = synthonSet.getSynthons();
            if (synthons.isEmpty()) {
                return new SynthonSetDownsamplingResult(synthonSet.toFragType(),
                        Collections.emptyList(),
                        Collections.emptyList(),
                        0);
            }

            List<RawSynthon> working = new ArrayList<>(synthons);
            Collections.shuffle(working, new Random(request.getRandomSeed()));

            int maxCenters = request.getEffectiveMaxCenters(working.size());
            boolean includeClusterMembers = request.isIncludeClusterMembers();
            List<CenterRecord> centers = new ArrayList<>();
            long startNanos = System.nanoTime();
            long intervalNanos = request.getProgressReportIntervalSeconds() * 1_000_000_000L;
            long lastReportNanos = startNanos;
            int skippedDescriptors = 0;
            reportProgress(synthonSet, request, 0, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
            for (int idx = 0; idx < working.size(); idx++) {
                RawSynthon candidate = working.get(idx);
                int processed = idx + 1;
                int[] descriptor = descriptorCache.getOrCompute(candidate);
                if (descriptor == null) {
                    skippedDescriptors++;
                    long now = System.nanoTime();
                    if (shouldReport(now, lastReportNanos, intervalNanos)) {
                        reportProgress(synthonSet, request, processed, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
                        lastReportNanos = now;
                    }
                    continue;
                }
                CenterWithSimilarity best = findBestCenter(candidate, descriptor, centers, request.isEnforceConnectorEquivalence());
                boolean belowThreshold = best == null || best.similarity < request.getMinSimilarity();
                boolean canCreateCenter = maxCenters <= 0 || centers.size() < maxCenters;
                if ((best == null || belowThreshold) && canCreateCenter) {
                    centers.add(new CenterRecord(candidate, descriptor, includeClusterMembers));
                } else if (best != null) {
                    best.center.registerAssignment(candidate, best.similarity);
                } else if (!centers.isEmpty()) {
                    CenterWithSimilarity fallback = findBestCenter(candidate, descriptor, centers, false);
                    if (fallback != null) {
                        fallback.center.registerAssignment(candidate, fallback.similarity);
                    }
                } else {
                    centers.add(new CenterRecord(candidate, descriptor, includeClusterMembers));
                }
                long now = System.nanoTime();
                if (shouldReport(now, lastReportNanos, intervalNanos)) {
                    reportProgress(synthonSet, request, processed, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
                    lastReportNanos = now;
                }
            }
            reportProgress(synthonSet, request, working.size(), working.size(), centers, maxCenters, startNanos, skippedDescriptors, true);

            List<SynthonSpace.FragId> representatives = new ArrayList<>();
            List<SynthonSetDownsamplingResult.ClusterInfo> clusters = new ArrayList<>();
            for (CenterRecord center : centers) {
                SynthonSpace.FragId frag = center.representative.toFragId();
                representatives.add(frag);
                clusters.add(center.toClusterInfo());
            }

            return new SynthonSetDownsamplingResult(synthonSet.toFragType(),
                    representatives,
                    clusters,
                    synthons.size());
        } finally {
            descriptorCache.clear();
        }
    }

    private boolean shouldReport(long nowNanos, long lastReportNanos, long intervalNanos) {
        return intervalNanos > 0 && nowNanos - lastReportNanos >= intervalNanos;
    }

    private void reportProgress(RawSynthonSet synthonSet,
                                SynthonDownsamplingRequest request,
                                int processed,
                                int total,
                                List<CenterRecord> centers,
                                int maxCenters,
                                long startNanos,
                                int skippedDescriptors,
                                boolean completed) {
        try {
            request.getProgressReporter().onProgress(new SynthonDownsamplingRequest.ProgressEvent(
                    synthonSet.getReactionId(),
                    synthonSet.getFragmentIndex(),
                    processed,
                    total,
                    centers.size(),
                    maxCenters,
                    (System.nanoTime() - startNanos) / 1_000_000L,
                    skippedDescriptors,
                    completed));
        } catch (RuntimeException ignored) {
            // Progress reporting is diagnostic only and must not affect downsampling results.
        }
    }

    private CenterWithSimilarity findBestCenter(RawSynthon candidate,
                                                int[] descriptor,
                                                List<CenterRecord> centers,
                                                boolean enforceConnectorEquivalence) {
        CenterRecord best = null;
        double bestSimilarity = -1.0;
        BitSet connectors = enforceConnectorEquivalence ? candidate.getConnectors() : null;
        for (CenterRecord center : centers) {
            if (enforceConnectorEquivalence && !equalConnectors(connectors, center.representative.getConnectors())) {
                continue;
            }
            double similarity = descriptorHandlers.get().getSimilarity(descriptor, center.descriptor);
            if (similarity > bestSimilarity) {
                bestSimilarity = similarity;
                best = center;
            }
        }
        if (best == null) {
            return null;
        }
        return new CenterWithSimilarity(best, bestSimilarity);
    }

    private boolean equalConnectors(BitSet a, BitSet b) {
        if (a == null || b == null) {
            return false;
        }
        return a.equals(b);
    }

    private int[] createDescriptor(RawSynthon synthon) {
        try {
            StereoMolecule mol = parsers.get().getCompactMolecule(synthon.getIdcode());
            mol.ensureHelperArrays(Molecule.cHelperCIP);
            return descriptorHandlers.get().createDescriptor(mol);
        } catch (Exception e) {
            return null;
        }
    }

    private static final class CenterRecord {
        private final RawSynthon representative;
        private final int[] descriptor;
        private final List<SynthonSetDownsamplingResult.ClusterMember> memberAssignments;
        private int members = 1;
        private double minSimilarityWithinCluster = 1.0;

        private CenterRecord(RawSynthon representative,
                             int[] descriptor,
                             boolean includeClusterMembers) {
            this.representative = representative;
            this.descriptor = descriptor;
            this.memberAssignments = includeClusterMembers ? new ArrayList<>() : null;
            if (this.memberAssignments != null) {
                this.memberAssignments.add(new SynthonSetDownsamplingResult.ClusterMember(representative.toFragId(), 1.0));
            }
        }

        private void registerAssignment(RawSynthon member, double similarity) {
            members++;
            minSimilarityWithinCluster = Math.min(minSimilarityWithinCluster, similarity);
            if (memberAssignments != null) {
                memberAssignments.add(new SynthonSetDownsamplingResult.ClusterMember(member.toFragId(), similarity));
            }
        }

        private SynthonSetDownsamplingResult.ClusterInfo toClusterInfo() {
            SynthonSpace.FragId frag = representative.toFragId();
            if (memberAssignments == null) {
                return new SynthonSetDownsamplingResult.ClusterInfo(frag, members, minSimilarityWithinCluster);
            }
            return new SynthonSetDownsamplingResult.ClusterInfo(frag,
                    members,
                    minSimilarityWithinCluster,
                    memberAssignments);
        }
    }

    private static final class CenterWithSimilarity {
        private final CenterRecord center;
        private final double similarity;

        private CenterWithSimilarity(CenterRecord center, double similarity) {
            this.center = center;
            this.similarity = similarity;
        }
    }
}
