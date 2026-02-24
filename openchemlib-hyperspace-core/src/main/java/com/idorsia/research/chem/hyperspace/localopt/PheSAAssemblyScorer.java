package com.idorsia.research.chem.hyperspace.localopt;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

/**
 * Default AssemblyScorer implementation that aligns assemblies with a query PheSA descriptor.
 */
public class PheSAAssemblyScorer implements AssemblyScorer {

    private final PheSAMolecule queryDescriptor;
    private final double minSimilarity;
    private final LongAdder comparisonCounter;
    private final ThreadLocal<DescriptorHandlerShape> descriptorHandlers = ThreadLocal.withInitial(DescriptorHandlerShape::new);
    private final ThreadLocal<ConformerSetGenerator> conformerGenerators = ThreadLocal.withInitial(() -> new ConformerSetGenerator(1));

    public PheSAAssemblyScorer(PheSAMolecule queryDescriptor, double minSimilarity) {
        this(queryDescriptor, minSimilarity, null);
    }

    public PheSAAssemblyScorer(PheSAMolecule queryDescriptor,
                               double minSimilarity,
                               LongAdder comparisonCounter) {
        this.queryDescriptor = queryDescriptor;
        this.minSimilarity = minSimilarity;
        this.comparisonCounter = comparisonCounter;
    }

    @Override
    public LocalOptimizationResult.BeamEntry score(String reactionId,
                                                   List<SynthonSpace.FragId> fragments,
                                                   int originatingRound,
                                                   LocalOptimizationLogger logger) {
        try {
            List<StereoMolecule> parts = new ArrayList<>(fragments.size());
            for (SynthonSpace.FragId frag : fragments) {
                StereoMolecule mol = new StereoMolecule();
                SynchronizedIDCodeParser.parse(mol, frag.idcode);
                mol.ensureHelperArrays(Molecule.cHelperCIP);
                parts.add(mol);
            }
            StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts);
            assembled.ensureHelperArrays(Molecule.cHelperCIP);
            int atomCount = assembled.getAtoms();
            int rotatable = countRotatableBonds(assembled);
            ConformerSet conformers = conformerGenerators.get().generateConformerSet(assembled);
            if (conformers.isEmpty()) {
                return null;
            }
            DescriptorHandlerShape handler = descriptorHandlers.get();
            PheSAMolecule descriptor = handler.createDescriptor(conformers);
            if (handler.calculationFailed(descriptor)) {
                return null;
            }
            double similarity = handler.getSimilarity(queryDescriptor, descriptor);
            if (comparisonCounter != null) {
                comparisonCounter.increment();
            }
            if (similarity < minSimilarity) {
                return null;
            }
            List<String> fragmentIds = fragments.stream().map(f -> f.fragment_id).collect(Collectors.toList());
            logger.logCandidate(fragmentIds, similarity, assembled.getIDCode(), originatingRound);
            return new LocalOptimizationResult.BeamEntry(fragments,
                    similarity,
                    atomCount,
                    rotatable,
                    assembled.getIDCode(),
                    originatingRound);
        } catch (Exception e) {
            return null;
        }
    }

    private int countRotatableBonds(StereoMolecule molecule) {
        molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
        boolean[] rotatable = new boolean[molecule.getBonds()];
        TorsionDB.findRotatableBonds(molecule, true, rotatable);
        int count = 0;
        for (boolean b : rotatable) {
            if (b) {
                count++;
            }
        }
        return count;
    }
}
