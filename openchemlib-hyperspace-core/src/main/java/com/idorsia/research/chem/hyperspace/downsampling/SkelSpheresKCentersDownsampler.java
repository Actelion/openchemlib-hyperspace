package com.idorsia.research.chem.hyperspace.downsampling;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.CachedStereoMoleculeProvider;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.Random;

/**
 * Simple online-like k-centers variant that uses SkelSpheres to pick representatives.
 */
public class SkelSpheresKCentersDownsampler implements SynthonDownsampler {

    private final CachedStereoMoleculeProvider stereoMoleculeProvider = new CachedStereoMoleculeProvider();
    private final SynthonDescriptorCache<int[]> descriptorCache;
    private final ThreadLocal<DescriptorHandlerBinarySkelSpheres> descriptorHandler;

    public SkelSpheresKCentersDownsampler() {
        this.descriptorHandler = ThreadLocal.withInitial(DescriptorHandlerBinarySkelSpheres::new);
        this.descriptorCache = new SynthonDescriptorCache<>(this::createDescriptor);
    }

    @Override
    public String getName() {
        return "SkelSpheresKCenters";
    }

    @Override
    public SynthonSetDownsamplingResult downsample(SynthonSpace sourceSpace,
                                                   SynthonSpace.FragType fragType,
                                                   List<SynthonSpace.FragId> synthons,
                                                   SynthonDownsamplingRequest request) {
        Objects.requireNonNull(fragType, "fragType");
        Objects.requireNonNull(synthons, "synthons");
        Objects.requireNonNull(request, "request");

        if (synthons.isEmpty()) {
            return new SynthonSetDownsamplingResult(fragType, Collections.emptyList(), Collections.emptyList(), 0);
        }

        List<SynthonSpace.FragId> working = new ArrayList<>(synthons);
        Collections.shuffle(working, new Random(request.getRandomSeed()));

        int maxCenters = request.getEffectiveMaxCenters(working.size());
        boolean includeClusterMembers = request.isIncludeClusterMembers();
        List<CenterRecord> centers = new ArrayList<>();
        long startNanos = System.nanoTime();
        long intervalNanos = request.getProgressReportIntervalSeconds() * 1_000_000_000L;
        long lastReportNanos = startNanos;
        int skippedDescriptors = 0;
        reportProgress(fragType, request, 0, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
        for (int idx = 0; idx < working.size(); idx++) {
            SynthonSpace.FragId candidate = working.get(idx);
            int processed = idx + 1;
            int[] descriptor = descriptorCache.getOrCompute(candidate);
            if (descriptor == null) {
                skippedDescriptors++;
                long now = System.nanoTime();
                if (shouldReport(now, lastReportNanos, intervalNanos)) {
                    reportProgress(fragType, request, processed, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
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
                // all centers filtered by connector equivalence but limit reached -> assign to overall best ignoring connectors
                CenterWithSimilarity fallback = findBestCenter(candidate, descriptor, centers, false);
                if (fallback != null) {
                    fallback.center.registerAssignment(candidate, fallback.similarity);
                }
            } else {
                centers.add(new CenterRecord(candidate, descriptor, includeClusterMembers));
            }
            long now = System.nanoTime();
            if (shouldReport(now, lastReportNanos, intervalNanos)) {
                reportProgress(fragType, request, processed, working.size(), centers, maxCenters, startNanos, skippedDescriptors, false);
                lastReportNanos = now;
            }
        }
        reportProgress(fragType, request, working.size(), working.size(), centers, maxCenters, startNanos, skippedDescriptors, true);

        List<SynthonSpace.FragId> representatives = new ArrayList<>();
        List<SynthonSetDownsamplingResult.ClusterInfo> clusterInfos = new ArrayList<>();
        for (CenterRecord center : centers) {
            representatives.add(center.representative);
            clusterInfos.add(center.toClusterInfo());
        }

        return new SynthonSetDownsamplingResult(fragType, representatives, clusterInfos, synthons.size());
    }

    private boolean shouldReport(long nowNanos, long lastReportNanos, long intervalNanos) {
        return intervalNanos > 0 && nowNanos - lastReportNanos >= intervalNanos;
    }

    private void reportProgress(SynthonSpace.FragType fragType,
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
                    fragType.rxn_id,
                    fragType.frag,
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

    private CenterWithSimilarity findBestCenter(SynthonSpace.FragId candidate,
                                                int[] descriptor,
                                                List<CenterRecord> centers,
                                                boolean enforceConnectorEquivalence) {
        CenterRecord best = null;
        double bestSimilarity = -1.0;
        BitSet candidateConnectors = enforceConnectorEquivalence ? candidate.getConnectors() : null;
        for (CenterRecord center : centers) {
            if (enforceConnectorEquivalence && !connectorsEqual(candidateConnectors, center.representative.getConnectors())) {
                continue;
            }
            double similarity = descriptorHandler.get().getSimilarity(descriptor, center.descriptor);
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

    private boolean connectorsEqual(BitSet a, BitSet b) {
        if (a == null || b == null) {
            return false;
        }
        return a.equals(b);
    }

    private int[] createDescriptor(SynthonSpace.FragId fragId) {
        StereoMolecule molecule = stereoMoleculeProvider.parseIDCode(fragId.idcode);
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return descriptorHandler.get().createDescriptor(molecule);
    }

    private static final class CenterRecord {
        private final SynthonSpace.FragId representative;
        private final int[] descriptor;
        private final List<SynthonSetDownsamplingResult.ClusterMember> memberAssignments;
        private int members = 1;
        private double minSimilarityWithinCluster = 1.0;

        private CenterRecord(SynthonSpace.FragId representative,
                             int[] descriptor,
                             boolean includeClusterMembers) {
            this.representative = representative;
            this.descriptor = descriptor;
            this.memberAssignments = includeClusterMembers ? new ArrayList<>() : null;
            if (this.memberAssignments != null) {
                this.memberAssignments.add(new SynthonSetDownsamplingResult.ClusterMember(representative, 1.0));
            }
        }

        private void registerAssignment(SynthonSpace.FragId member, double similarity) {
            members++;
            minSimilarityWithinCluster = Math.min(minSimilarityWithinCluster, similarity);
            if (memberAssignments != null) {
                memberAssignments.add(new SynthonSetDownsamplingResult.ClusterMember(member, similarity));
            }
        }

        private SynthonSetDownsamplingResult.ClusterInfo toClusterInfo() {
            if (memberAssignments == null) {
                return new SynthonSetDownsamplingResult.ClusterInfo(representative, members, minSimilarityWithinCluster);
            }
            return new SynthonSetDownsamplingResult.ClusterInfo(representative,
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
