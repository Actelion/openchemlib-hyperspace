package com.idorsia.research.chem.hyperspace3d.feature;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import java.util.ArrayDeque;
import java.util.Arrays;

/** OCL-native implementation of the Deepspace7 V2 50D graph contract. */
public final class OCLDeepSpaceFeaturizer {
    private static final int N = DeepSpaceFeatureSchema.MAX_ATOMS;
    private static final int AD = DeepSpaceFeatureSchema.ATOM_FEATURE_DIM;
    private static final int PD = DeepSpaceFeatureSchema.PAIR_FEATURE_DIM;

    public FeaturizationResult featurize(StereoMolecule input) {
        if (input == null) return FeaturizationResult.rejected("NULL_MOLECULE");
        StereoMolecule mol = removeExplicitHydrogens(input);
        mol.ensureHelperArrays(Molecule.cHelperCIP);
        String invalid = validate(mol);
        if (invalid != null) return FeaturizationResult.rejected(invalid);

        int atomCount = mol.getAtoms();
        float[] atom = new float[N * AD];
        float[] pair = new float[N * N * PD];
        boolean[] mask = new boolean[N];
        Arrays.fill(mask, 0, atomCount, true);

        RingInfo rings = RingInfo.of(mol);
        for (int a = 0; a < atomCount; a++) {
            encodeAtom(mol, rings, a, atom);
        }
        int[][] distances = graphDistances(mol, atomCount);
        for (int a = 0; a < atomCount; a++) {
            for (int b = 0; b < atomCount; b++) {
                encodePair(mol, rings, distances, a, b, pair);
            }
        }
        return new FeaturizationResult(atom, pair, mask, null);
    }

    private static StereoMolecule removeExplicitHydrogens(StereoMolecule input) {
        StereoMolecule mol = new StereoMolecule(input);
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            if (mol.getAtomicNo(atom) == 1) mol.markAtomForDeletion(atom);
        }
        mol.deleteMarkedAtomsAndBonds();
        return mol;
    }

    private static String validate(StereoMolecule mol) {
        int atoms = mol.getAtoms();
        if (atoms < 6) return "TOO_FEW_HEAVY_ATOMS";
        if (atoms > N) return "TOO_MANY_HEAVY_ATOMS";
        int charge = 0;
        boolean carbon = false;
        for (int a = 0; a < atoms; a++) {
            int atomicNo = mol.getAtomicNo(a);
            if (!DeepSpaceFeatureSchema.ELEMENT_SET.contains(atomicNo)) {
                return "UNSUPPORTED_ELEMENT:" + mol.getAtomLabel(a);
            }
            carbon |= atomicNo == 6;
            charge += mol.getAtomCharge(a);
            if (mol.getAtomCharge(a) < -3 || mol.getAtomCharge(a) > 3) {
                return "UNSUPPORTED_ATOM_CHARGE:" + mol.getAtomCharge(a);
            }
            if (mol.getAtomRadical(a) != Molecule.cAtomRadicalStateNone) return "RADICAL";
            if (mol.getAtomMass(a) != 0 && !mol.isNaturalAbundance(a)) return "ISOTOPE";
            if (mol.getAtomMapNo(a) != 0) return "ATOM_MAP";
        }
        if (!carbon) return "NO_CARBON";
        if (charge < -2 || charge > 2) return "UNSUPPORTED_TOTAL_CHARGE:" + charge;
        if (!connected(mol, atoms)) return "MULTIPLE_FRAGMENTS";
        return null;
    }

    private static boolean connected(StereoMolecule mol, int atoms) {
        boolean[] seen = new boolean[atoms];
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        seen[0] = true;
        queue.add(0);
        int count = 0;
        while (!queue.isEmpty()) {
            int a = queue.removeFirst();
            count++;
            for (int i = 0; i < mol.getConnAtoms(a); i++) {
                int b = mol.getConnAtom(a, i);
                if (!seen[b]) { seen[b] = true; queue.addLast(b); }
            }
        }
        return count == atoms;
    }

    private static void encodeAtom(StereoMolecule mol, RingInfo rings,
                                   int a, float[] out) {
        int base = a * AD;
        int element = DeepSpaceFeatureSchema.ELEMENTS.indexOf(mol.getAtomicNo(a));
        out[base + 2 + element] = 1f;

        int charge = mol.getAtomCharge(a);
        out[base + 13 + (charge >= -3 && charge <= 3 ? charge + 3 : 7)] = 1f;
        out[base + 21 + hybridization(mol, a)] = 1f;
        out[base + 27 + chirality(mol, a)] = 1f;
        out[base + 31 + ringClass(rings.smallest[a])] = 1f;
        out[base + 39 + Math.min(rings.count[a], 3)] = 1f;
        out[base + 43] = mol.getConnAtoms(a) / 4f;
        out[base + 44] = (mol.getOccupiedValence(a) + mol.getAllHydrogens(a)) / 6f;
        out[base + 45] = mol.getAllHydrogens(a) / 4f;
        out[base + 46] = mol.isAromaticAtom(a) ? 1f : 0f;
        out[base + 47] = rings.count[a] > 0 ? 1f : 0f;
        out[base + 48] = rings.small[a] ? 1f : 0f;
        out[base + 49] = rings.large[a] ? 1f : 0f;
    }

    private static int hybridization(StereoMolecule mol, int atom) {
        int pi = mol.getAtomPi(atom);
        int z = mol.getAtomicNo(atom);
        if (mol.isAromaticAtom(atom)) return 1;
        int coordination = mol.getConnAtoms(atom) + mol.getAllHydrogens(atom);
        if (z == 15 || z == 16) return coordination >= 6 ? 4 : coordination == 5 ? 3 : 2;
        if (pi >= 2) return 0;
        if (pi == 1 || ((z == 7 || z == 8)
                && hasPiNeighbor(mol, atom))) return 1;
        if (coordination >= 6) return 4;
        if (coordination == 5) return 3;
        if (coordination <= 4) return 2;
        return 5;
    }

    private static boolean hasPiNeighbor(StereoMolecule mol, int atom) {
        for (int i = 0; i < mol.getConnAtoms(atom); i++) {
            int neighbor = mol.getConnAtom(atom, i);
            if (mol.getAtomicNo(neighbor) == 6
                    && (mol.getAtomPi(neighbor) > 0 || mol.isAromaticAtom(neighbor))) return true;
        }
        return false;
    }

    private static int chirality(StereoMolecule mol, int atom) {
        return switch (mol.getAtomCIPParity(atom)) {
            case Molecule.cAtomCIPParityNone -> 0;
            case Molecule.cAtomCIPParityRorM -> 1;
            case Molecule.cAtomCIPParitySorP -> 2;
            default -> 3;
        };
    }

    private static int ringClass(int size) {
        if (size == 0) return 0;
        if (size >= 3 && size <= 8) return size - 2;
        return 7;
    }

    private static void encodePair(StereoMolecule mol, RingInfo rings, int[][] distances,
                                   int a, int b, float[] out) {
        int base = ((a * N) + b) * PD;
        out[base] = a == b ? 1f : 0f;
        int bond = a == b ? -1 : mol.getBond(a, b);
        out[base + 1] = bond >= 0 ? 1f : 0f;
        out[base + 2 + bondType(mol, bond)] = 1f;
        out[base + 7 + bondStereo(mol, bond)] = 1f;
        out[base + 13 + Math.min(distances[a][b], 9)] = 1f;
        int shared = rings.shared[a][b];
        out[base + 23 + ringClass(shared)] = 1f;
        if (bond >= 0) {
            out[base + 31] = isConjugated(mol, bond) ? 1f : 0f;
            out[base + 32] = mol.isAromaticBond(bond) ? 1f : 0f;
        }
        out[base + 33] = shared > 0 ? 1f : 0f;
        out[base + 34] = shared >= 3 && shared <= 8 ? 1f : 0f;
        out[base + 35] = shared > 8 ? 1f : 0f;
    }

    private static int bondType(StereoMolecule mol, int bond) {
        if (bond < 0) return 0;
        if (mol.isAromaticBond(bond)) return 4;
        return switch (mol.getBondOrder(bond)) {
            case 1 -> 1;
            case 2 -> 2;
            case 3 -> 3;
            default -> 0;
        };
    }

    private static int bondStereo(StereoMolecule mol, int bond) {
        if (bond < 0) return 0;
        return switch (mol.getBondCIPParity(bond)) {
            case Molecule.cBondCIPParityEorP -> 1;
            case Molecule.cBondCIPParityZorM -> 2;
            case Molecule.cBondCIPParityProblem -> 5;
            default -> 0;
        };
    }

    private static boolean isConjugated(StereoMolecule mol, int bond) {
        if (mol.isAromaticBond(bond) || mol.isDelocalizedBond(bond)) return true;
        int a = mol.getBondAtom(0, bond);
        int b = mol.getBondAtom(1, bond);
        if (mol.getAtomicNo(a) == 15 || mol.getAtomicNo(a) == 16
                || mol.getAtomicNo(b) == 15 || mol.getAtomicNo(b) == 16) return false;
        if (mol.getBondOrder(bond) > 1) return hasAdjacentPi(mol, a, bond)
                || hasAdjacentPi(mol, b, bond) || hasAdjacentHetero(mol, a, bond)
                || hasAdjacentHetero(mol, b, bond) || mol.isStabilizedAtom(a)
                || mol.isStabilizedAtom(b);
        return (mol.getAtomPi(a) > 0 && mol.getAtomPi(b) > 0
                && mol.getAtomicNo(a) != 15 && mol.getAtomicNo(a) != 16 && mol.getAtomicNo(b) != 15 && mol.getAtomicNo(b) != 16)
                || (mol.getAtomPi(a) > 0 && mol.getAtomicNo(a) != 15 && mol.getAtomicNo(a) != 16
                    && mol.getAtomicNo(b) != 6 && mol.isStabilizedAtom(b))
                || (mol.getAtomPi(b) > 0 && mol.getAtomicNo(b) != 15 && mol.getAtomicNo(b) != 16
                    && mol.getAtomicNo(a) != 6 && (mol.isStabilizedAtom(a) || mol.getAtomCharge(a) != 0))
                || (mol.isAromaticAtom(a) && (mol.getAtomicNo(b) == 7 || mol.getAtomicNo(b) == 8))
                || (mol.isAromaticAtom(b) && (mol.getAtomicNo(a) == 7 || mol.getAtomicNo(a) == 8));
    }

    private static boolean hasAdjacentPi(StereoMolecule mol, int atom, int excludedBond) {
        for (int i = 0; i < mol.getConnAtoms(atom); i++) {
            int bond = mol.getConnBond(atom, i);
            if (bond != excludedBond
                    && (mol.getBondOrder(bond) > 1 || mol.isAromaticBond(bond))) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasAdjacentHetero(StereoMolecule mol, int atom,
                                              int excludedBond) {
        for (int i = 0; i < mol.getConnAtoms(atom); i++) {
            int bond = mol.getConnBond(atom, i);
            int neighbor = mol.getConnAtom(atom, i);
            if (bond != excludedBond && mol.getBondOrder(bond) == 1
                    && mol.getAtomicNo(neighbor) != 6) return true;
        }
        return false;
    }

    private static int[][] graphDistances(StereoMolecule mol, int atoms) {
        int[][] result = new int[atoms][atoms];
        for (int start = 0; start < atoms; start++) {
            Arrays.fill(result[start], 99);
            result[start][start] = 0;
            ArrayDeque<Integer> queue = new ArrayDeque<>();
            queue.add(start);
            while (!queue.isEmpty()) {
                int a = queue.removeFirst();
                for (int i = 0; i < mol.getConnAtoms(a); i++) {
                    int b = mol.getConnAtom(a, i);
                    if (result[start][b] == 99) {
                        result[start][b] = result[start][a] + 1;
                        queue.addLast(b);
                    }
                }
            }
        }
        return result;
    }

    private static final class RingInfo {
        final int[] smallest;
        final int[] count;
        final boolean[] small;
        final boolean[] large;
        final int[][] shared;

        private RingInfo(int atoms) {
            smallest = new int[atoms];
            count = new int[atoms];
            small = new boolean[atoms];
            large = new boolean[atoms];
            shared = new int[atoms][atoms];
        }

        static RingInfo of(StereoMolecule mol) {
            int atoms = mol.getAtoms();
            RingInfo info = new RingInfo(atoms);
            RingCollection rings = mol.getRingSet();
            for (int r = 0; r < rings.getSize(); r++) {
                int size = rings.getRingSize(r);
                int[] members = rings.getRingAtoms(r);
                for (int atom : members) {
                    info.count[atom]++;
                    if (info.smallest[atom] == 0 || size < info.smallest[atom]) info.smallest[atom] = size;
                    info.small[atom] |= size <= 8;
                    info.large[atom] |= size > 8;
                }
                for (int a : members) for (int b : members) {
                    if (info.shared[a][b] == 0 || size < info.shared[a][b]) info.shared[a][b] = size;
                }
            }
            return info;
        }
    }
}
