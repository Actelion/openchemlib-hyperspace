package com.idorsia.research.chem.hyperspace.tools.chembl;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/**
 * One atom-mapped, overlapping decomposition used to train a latent molecular assembler.
 * Product atom indices refer to the atom order obtained by parsing {@link #productIdcode()}.
 */
public record LatentAssemblyDecomposition(
        String schemaVersion,
        String sourceMoleculeId,
        String productIdcode,
        int productAtomCount,
        int productNonHydrogenAtomCount,
        List<Integer> scaffoldProductAtoms,
        List<Arm> arms,
        View scaffoldView,
        List<View> assemblyViews,
        List<CutBond> cutBonds,
        SelectionMetrics selectionMetrics
) {
    public static final String SCHEMA_VERSION = "latent-assembly-decomposition-v1";

    public LatentAssemblyDecomposition {
        schemaVersion = Objects.requireNonNull(schemaVersion, "schemaVersion");
        sourceMoleculeId = Objects.requireNonNull(sourceMoleculeId, "sourceMoleculeId");
        productIdcode = Objects.requireNonNull(productIdcode, "productIdcode");
        scaffoldProductAtoms = List.copyOf(scaffoldProductAtoms);
        arms = List.copyOf(arms);
        assemblyViews = List.copyOf(assemblyViews);
        cutBonds = List.copyOf(cutBonds);
        scaffoldView = Objects.requireNonNull(scaffoldView, "scaffoldView");
        selectionMetrics = Objects.requireNonNull(selectionMetrics, "selectionMetrics");

        if (!SCHEMA_VERSION.equals(schemaVersion)) {
            throw new IllegalArgumentException("Unsupported schema version: " + schemaVersion);
        }
        if (arms.size() < 2 || arms.size() > 3) {
            throw new IllegalArgumentException("A decomposition must contain two or three arms");
        }
        if (assemblyViews.size() != arms.size() || cutBonds.size() != arms.size()) {
            throw new IllegalArgumentException("Every arm requires one assembly view and one cut bond");
        }
        if (productAtomCount <= 0 || productNonHydrogenAtomCount <= 0
                || productNonHydrogenAtomCount > productAtomCount) {
            throw new IllegalArgumentException("Invalid product atom counts");
        }

        Set<Integer> partition = new HashSet<>();
        addPartitionAtoms(partition, scaffoldProductAtoms, productAtomCount, "scaffold");
        for (int i = 0; i < arms.size(); i++) {
            Arm arm = arms.get(i);
            if (arm.index() != i) {
                throw new IllegalArgumentException("Arm indices must be consecutive and ordered");
            }
            addPartitionAtoms(partition, arm.productAtoms(), productAtomCount, "arm " + i);
        }
        if (partition.size() != productAtomCount) {
            throw new IllegalArgumentException("Scaffold and arms do not partition all product atoms");
        }
    }

    private static void addPartitionAtoms(Set<Integer> partition,
                                          List<Integer> atoms,
                                          int productAtomCount,
                                          String label) {
        if (atoms.isEmpty()) {
            throw new IllegalArgumentException(label + " is empty");
        }
        for (int atom : atoms) {
            if (atom < 0 || atom >= productAtomCount) {
                throw new IllegalArgumentException(label + " contains an invalid product atom: " + atom);
            }
            if (!partition.add(atom)) {
                throw new IllegalArgumentException("Product atom occurs in multiple components: " + atom);
            }
        }
    }

    public record Arm(int index, String label, List<Integer> productAtoms) {
        public Arm {
            label = Objects.requireNonNull(label, "label");
            productAtoms = List.copyOf(productAtoms);
        }
    }

    /**
     * A canonical, hydrogen-capped molecular view. The mapping list is indexed by
     * the atom order obtained when reparsing {@code idcode}.
     */
    public record View(
            String label,
            String idcode,
            List<Integer> productAtomByViewAtom,
            List<Integer> openScaffoldAttachmentViewAtoms
    ) {
        public View {
            label = Objects.requireNonNull(label, "label");
            idcode = Objects.requireNonNull(idcode, "idcode");
            productAtomByViewAtom = List.copyOf(productAtomByViewAtom);
            openScaffoldAttachmentViewAtoms = List.copyOf(openScaffoldAttachmentViewAtoms);
            Set<Integer> mapped = new HashSet<>(productAtomByViewAtom);
            if (mapped.size() != productAtomByViewAtom.size()) {
                throw new IllegalArgumentException("A view contains duplicate product atom mappings");
            }
            for (int viewAtom : openScaffoldAttachmentViewAtoms) {
                if (viewAtom < 0 || viewAtom >= productAtomByViewAtom.size()) {
                    throw new IllegalArgumentException("Invalid open attachment view atom: " + viewAtom);
                }
            }
        }
    }

    public record CutBond(
            int productBondIndex,
            int scaffoldProductAtom,
            int armProductAtom,
            int armIndex,
            int bondType,
            int bondOrder
    ) {}

    public record SelectionMetrics(
            int scaffoldNonHydrogenAtoms,
            int scaffoldRingAtoms,
            List<Integer> armNonHydrogenAtoms,
            int minimumAttachmentDistance,
            double scaffoldFraction
    ) {
        public SelectionMetrics {
            armNonHydrogenAtoms = List.copyOf(armNonHydrogenAtoms);
        }
    }

    static List<Integer> atomList(java.util.BitSet atoms) {
        List<Integer> result = new ArrayList<>(atoms.cardinality());
        for (int atom = atoms.nextSetBit(0); atom >= 0; atom = atoms.nextSetBit(atom + 1)) {
            result.add(atom);
        }
        return List.copyOf(result);
    }
}
