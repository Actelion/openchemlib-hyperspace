package com.idorsia.research.chem.hyperspace.downsampling;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.CachedStereoMoleculeProvider;
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

    private final CachedStereoMoleculeProvider moleculeProvider = new CachedStereoMoleculeProvider();
    private final RawSynthonDescriptorCache<int[]> descriptorCache;
    private final ThreadLocal<DescriptorHandlerBinarySkelSpheres> descriptorHandlers;

    public SkelSpheresKCentersRawDownsampler() {
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
        List<RawSynthon> synthons = synthonSet.getSynthons();
        if (synthons.isEmpty()) {
            return new SynthonSetDownsamplingResult(synthonSet.toFragType(),
                    Collections.emptyList(),
                    Collections.emptyList(),
                    0);
        }

        List<RawSynthon> working = new ArrayList<>(synthons);
        Collections.shuffle(working, new Random(request.getRandomSeed()));

        List<CenterRecord> centers = new ArrayList<>();
        for (RawSynthon candidate : working) {
            int[] descriptor = descriptorCache.getOrCompute(candidate);
            if (descriptor == null) {
                continue;
            }
            CenterWithSimilarity best = findBestCenter(candidate, descriptor, centers, request.isEnforceConnectorEquivalence());
            boolean belowThreshold = best == null || best.similarity < request.getMinSimilarity();
            boolean canCreateCenter = request.getMaxCenters() <= 0 || centers.size() < request.getMaxCenters();
            if ((best == null || belowThreshold) && canCreateCenter) {
                centers.add(new CenterRecord(candidate, descriptor));
                continue;
            }
            if (best != null) {
                best.center.registerAssignment(best.similarity);
            } else if (!centers.isEmpty()) {
                CenterWithSimilarity fallback = findBestCenter(candidate, descriptor, centers, false);
                if (fallback != null) {
                    fallback.center.registerAssignment(fallback.similarity);
                }
            } else {
                centers.add(new CenterRecord(candidate, descriptor));
            }
        }

        List<SynthonSpace.FragId> representatives = new ArrayList<>();
        List<SynthonSetDownsamplingResult.ClusterInfo> clusters = new ArrayList<>();
        for (CenterRecord center : centers) {
            SynthonSpace.FragId frag = center.representative.toFragId();
            representatives.add(frag);
            clusters.add(new SynthonSetDownsamplingResult.ClusterInfo(frag,
                    center.members,
                    center.minSimilarityWithinCluster));
        }

        return new SynthonSetDownsamplingResult(synthonSet.toFragType(),
                representatives,
                clusters,
                synthons.size());
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
            StereoMolecule mol = moleculeProvider.parseIDCode(synthon.getIdcode());
            mol.ensureHelperArrays(Molecule.cHelperCIP);
            return descriptorHandlers.get().createDescriptor(mol);
        } catch (Exception e) {
            return null;
        }
    }

    private static final class CenterRecord {
        private final RawSynthon representative;
        private final int[] descriptor;
        private int members = 1;
        private double minSimilarityWithinCluster = 1.0;

        private CenterRecord(RawSynthon representative, int[] descriptor) {
            this.representative = representative;
            this.descriptor = descriptor;
        }

        private void registerAssignment(double similarity) {
            members++;
            minSimilarityWithinCluster = Math.min(minSimilarityWithinCluster, similarity);
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
