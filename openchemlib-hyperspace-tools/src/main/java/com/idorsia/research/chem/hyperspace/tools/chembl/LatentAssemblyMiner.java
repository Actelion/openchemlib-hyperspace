package com.idorsia.research.chem.hyperspace.tools.chembl;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * Mines conservative A+B+scaffold and A+B+C+scaffold training decompositions.
 * Candidate cuts are non-ring heavy-atom single bonds. Retained cutsets must form
 * a component star whose center is a complete, ring-containing scaffold.
 */
public final class LatentAssemblyMiner {
    public static final int HARD_MAX_PRODUCT_NON_HYDROGEN_ATOMS = 32;

    public Result mine(Options options) throws IOException {
        List<LatentAssemblyDecomposition> decompositions = new ArrayList<>();
        Stats stats = mine(options, decompositions::add);
        return new Result(List.copyOf(decompositions), stats);
    }

    public Stats mine(Options options, DecompositionSink sink) throws IOException {
        Objects.requireNonNull(options, "options");
        Objects.requireNonNull(sink, "sink");
        MutableStats stats = new MutableStats();
        IDCodeParser parser = new IDCodeParser();

        try (BufferedReader reader = newReader(options.input())) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Input file is empty: " + options.input());
            }
            HeaderIndex columns = HeaderIndex.parse(header);
            int idColumn = columns.required(options.idColumn());
            int structureColumn = columns.required(options.structureColumn());

            String line;
            while ((line = reader.readLine()) != null) {
                if (options.maxMolecules() > 0 && stats.parsedMolecules >= options.maxMolecules()) {
                    break;
                }
                stats.inputRows++;
                String[] fields = line.split("\t", -1);
                if (idColumn >= fields.length || structureColumn >= fields.length) {
                    stats.skippedRows++;
                    continue;
                }
                String sourceMoleculeId = fields[idColumn].trim();
                String sourceIdcode = fields[structureColumn].trim();
                if (sourceMoleculeId.isEmpty() || sourceIdcode.isEmpty()) {
                    stats.skippedRows++;
                    continue;
                }

                StereoMolecule molecule;
                try {
                    molecule = canonicalMolecule(parser, sourceIdcode);
                } catch (RuntimeException ex) {
                    stats.parseFailures++;
                    continue;
                }
                stats.parsedMolecules++;

                int nonHydrogenAtoms = nonHydrogenAtomCount(molecule, allAtoms(molecule));
                if (nonHydrogenAtoms < options.minProductNonHydrogenAtoms()) {
                    stats.rejectedProductsTooSmall++;
                    continue;
                }
                if (nonHydrogenAtoms > options.maxProductNonHydrogenAtoms()) {
                    stats.rejectedProductsTooLarge++;
                    continue;
                }
                if (connectedComponents(molecule, new BitSet()).size() != 1) {
                    stats.rejectedDisconnectedProducts++;
                    continue;
                }
                if (options.onlyCommonDrugElements() && !containsOnlyCommonDrugElements(molecule)) {
                    stats.rejectedUncommonElements++;
                    continue;
                }
                stats.eligibleProducts++;

                List<Candidate> candidates = mineMolecule(molecule, options, stats);
                if (!candidates.isEmpty()) {
                    stats.productsWithDecompositions++;
                }
                int count = Math.min(options.maxDecompositionsPerMolecule(), candidates.size());
                for (int i = 0; i < count; i++) {
                    sink.accept(toDecomposition(sourceMoleculeId, molecule, nonHydrogenAtoms, candidates.get(i)));
                    stats.emittedDecompositions++;
                }
            }
        }
        return stats.freeze();
    }

    private static StereoMolecule canonicalMolecule(IDCodeParser parser, String idcode) {
        StereoMolecule parsed = new StereoMolecule();
        parser.parse(parsed, idcode);
        parsed.ensureHelperArrays(Molecule.cHelperCIP);
        String canonicalIdcode = parsed.getIDCode();
        StereoMolecule canonical = new StereoMolecule();
        parser.parse(canonical, canonicalIdcode);
        canonical.ensureHelperArrays(Molecule.cHelperCIP);
        return canonical;
    }

    private static List<Candidate> mineMolecule(StereoMolecule molecule,
                                                 Options options,
                                                 MutableStats stats) {
        List<Integer> eligibleBonds = candidateCutBonds(molecule);
        Map<String, Candidate> unique = new LinkedHashMap<>();
        if (options.armMode().includes(2)) {
            for (int i = 0; i < eligibleBonds.size(); i++) {
                for (int j = i + 1; j < eligibleBonds.size(); j++) {
                    Candidate candidate = analyzeCutset(molecule,
                            new int[]{eligibleBonds.get(i), eligibleBonds.get(j)}, options);
                    if (candidate != null) {
                        stats.acceptedTwoArmCandidates++;
                        unique.putIfAbsent(candidate.key(), candidate);
                    }
                }
            }
        }
        if (options.armMode().includes(3)) {
            for (int i = 0; i < eligibleBonds.size(); i++) {
                for (int j = i + 1; j < eligibleBonds.size(); j++) {
                    for (int k = j + 1; k < eligibleBonds.size(); k++) {
                        Candidate candidate = analyzeCutset(molecule,
                                new int[]{eligibleBonds.get(i), eligibleBonds.get(j), eligibleBonds.get(k)}, options);
                        if (candidate != null) {
                            stats.acceptedThreeArmCandidates++;
                            unique.putIfAbsent(candidate.key(), candidate);
                        }
                    }
                }
            }
        }
        List<Candidate> candidates = new ArrayList<>(unique.values());
        candidates.sort(Candidate.ORDER);
        return candidates;
    }

    private static Candidate analyzeCutset(StereoMolecule molecule, int[] cutBonds, Options options) {
        BitSet cutMask = new BitSet(molecule.getBonds());
        for (int bond : cutBonds) {
            cutMask.set(bond);
        }
        List<BitSet> components = connectedComponents(molecule, cutMask);
        if (components.size() != cutBonds.length + 1) {
            return null;
        }

        int[] componentByAtom = new int[molecule.getAtoms()];
        Arrays.fill(componentByAtom, -1);
        for (int component = 0; component < components.size(); component++) {
            for (int atom = components.get(component).nextSetBit(0);
                 atom >= 0;
                 atom = components.get(component).nextSetBit(atom + 1)) {
                componentByAtom[atom] = component;
            }
        }

        int[] componentDegree = new int[components.size()];
        for (int bond : cutBonds) {
            int component1 = componentByAtom[molecule.getBondAtom(0, bond)];
            int component2 = componentByAtom[molecule.getBondAtom(1, bond)];
            if (component1 == component2) {
                return null;
            }
            componentDegree[component1]++;
            componentDegree[component2]++;
        }

        int scaffoldComponent = -1;
        for (int component = 0; component < componentDegree.length; component++) {
            if (componentDegree[component] == cutBonds.length) {
                if (scaffoldComponent >= 0) {
                    return null;
                }
                scaffoldComponent = component;
            }
        }
        if (scaffoldComponent < 0) {
            return null;
        }
        for (int component = 0; component < componentDegree.length; component++) {
            if (component != scaffoldComponent && componentDegree[component] != 1) {
                return null;
            }
        }

        BitSet scaffold = components.get(scaffoldComponent);
        List<CandidateArm> arms = new ArrayList<>();
        for (int bond : cutBonds) {
            int atom1 = molecule.getBondAtom(0, bond);
            int atom2 = molecule.getBondAtom(1, bond);
            int scaffoldAtom = componentByAtom[atom1] == scaffoldComponent ? atom1 : atom2;
            int armAtom = scaffoldAtom == atom1 ? atom2 : atom1;
            int armComponent = componentByAtom[armAtom];
            arms.add(new CandidateArm(components.get(armComponent), bond, scaffoldAtom, armAtom));
        }
        arms.sort(Comparator.comparingInt(CandidateArm::scaffoldAtom)
                .thenComparingInt(CandidateArm::armAtom)
                .thenComparingInt(CandidateArm::cutBond));

        int productNonHydrogenAtoms = nonHydrogenAtomCount(molecule, allAtoms(molecule));
        int scaffoldNonHydrogenAtoms = nonHydrogenAtomCount(molecule, scaffold);
        int scaffoldRingAtoms = ringAtomCount(molecule, scaffold);
        double scaffoldFraction = scaffoldNonHydrogenAtoms / (double) productNonHydrogenAtoms;
        if (scaffoldNonHydrogenAtoms < options.minScaffoldNonHydrogenAtoms()
                || scaffoldNonHydrogenAtoms > options.maxScaffoldNonHydrogenAtoms()
                || scaffoldRingAtoms < options.minScaffoldRingAtoms()
                || scaffoldFraction < options.minScaffoldFraction()
                || scaffoldFraction > options.maxScaffoldFraction()) {
            return null;
        }

        Set<Integer> scaffoldAttachmentAtoms = new HashSet<>();
        List<Integer> armSizes = new ArrayList<>();
        for (CandidateArm arm : arms) {
            int armSize = nonHydrogenAtomCount(molecule, arm.atoms());
            if (armSize < options.minArmNonHydrogenAtoms()
                    || armSize > options.maxArmNonHydrogenAtoms()
                    || !scaffoldAttachmentAtoms.add(arm.scaffoldAtom())) {
                return null;
            }
            armSizes.add(armSize);
        }

        int minimumAttachmentDistance = Integer.MAX_VALUE;
        for (int i = 0; i < arms.size(); i++) {
            for (int j = i + 1; j < arms.size(); j++) {
                minimumAttachmentDistance = Math.min(minimumAttachmentDistance,
                        molecule.getPathLength(arms.get(i).scaffoldAtom(), arms.get(j).scaffoldAtom()));
            }
        }
        if (minimumAttachmentDistance < options.minAttachmentDistance()) {
            return null;
        }
        return new Candidate(scaffold, arms, scaffoldNonHydrogenAtoms, scaffoldRingAtoms,
                armSizes, minimumAttachmentDistance, scaffoldFraction);
    }

    private static LatentAssemblyDecomposition toDecomposition(String sourceMoleculeId,
                                                               StereoMolecule molecule,
                                                               int productNonHydrogenAtoms,
                                                               Candidate candidate) {
        List<LatentAssemblyDecomposition.Arm> arms = new ArrayList<>();
        List<LatentAssemblyDecomposition.CutBond> cuts = new ArrayList<>();
        List<Integer> allScaffoldAttachments = candidate.arms().stream()
                .map(CandidateArm::scaffoldAtom).distinct().sorted().toList();

        LatentAssemblyDecomposition.View scaffoldView = buildView(
                molecule, "S", candidate.scaffold(), allScaffoldAttachments);
        List<LatentAssemblyDecomposition.View> assemblyViews = new ArrayList<>();
        for (int i = 0; i < candidate.arms().size(); i++) {
            CandidateArm arm = candidate.arms().get(i);
            String armLabel = Character.toString((char) ('A' + i));
            arms.add(new LatentAssemblyDecomposition.Arm(i, armLabel,
                    LatentAssemblyDecomposition.atomList(arm.atoms())));

            BitSet viewAtoms = (BitSet) candidate.scaffold().clone();
            viewAtoms.or(arm.atoms());
            List<Integer> openAttachments = new ArrayList<>();
            for (int j = 0; j < candidate.arms().size(); j++) {
                if (j != i) {
                    openAttachments.add(candidate.arms().get(j).scaffoldAtom());
                }
            }
            assemblyViews.add(buildView(molecule, "S+" + armLabel, viewAtoms, openAttachments));
            cuts.add(new LatentAssemblyDecomposition.CutBond(
                    arm.cutBond(), arm.scaffoldAtom(), arm.armAtom(), i,
                    molecule.getBondType(arm.cutBond()), molecule.getBondOrder(arm.cutBond())));
        }

        return new LatentAssemblyDecomposition(
                LatentAssemblyDecomposition.SCHEMA_VERSION,
                sourceMoleculeId,
                molecule.getIDCode(),
                molecule.getAtoms(),
                productNonHydrogenAtoms,
                LatentAssemblyDecomposition.atomList(candidate.scaffold()),
                arms,
                scaffoldView,
                assemblyViews,
                cuts,
                new LatentAssemblyDecomposition.SelectionMetrics(
                        candidate.scaffoldNonHydrogenAtoms(),
                        candidate.scaffoldRingAtoms(),
                        candidate.armNonHydrogenAtoms(),
                        candidate.minimumAttachmentDistance(),
                        candidate.scaffoldFraction())
        );
    }

    private static LatentAssemblyDecomposition.View buildView(StereoMolecule product,
                                                              String label,
                                                              BitSet includedProductAtoms,
                                                              List<Integer> openScaffoldProductAtoms) {
        boolean[] include = new boolean[product.getAtoms()];
        for (int atom = includedProductAtoms.nextSetBit(0);
             atom >= 0;
             atom = includedProductAtoms.nextSetBit(atom + 1)) {
            include[atom] = true;
        }
        int[] productToOldViewAtom = new int[product.getAtoms()];
        StereoMolecule view = new StereoMolecule();
        product.copyMoleculeByAtoms(view, include, true, productToOldViewAtom);
        view.setFragment(false);
        view.ensureHelperArrays(Molecule.cHelperCIP);

        int[] productByOldViewAtom = new int[view.getAtoms()];
        Arrays.fill(productByOldViewAtom, -1);
        for (int productAtom = 0; productAtom < productToOldViewAtom.length; productAtom++) {
            int oldViewAtom = productToOldViewAtom[productAtom];
            if (oldViewAtom >= 0) {
                productByOldViewAtom[oldViewAtom] = productAtom;
            }
        }

        Canonizer canonizer = new Canonizer(view);
        String idcode = canonizer.getIDCode();
        int[] oldViewAtomByCanonicalAtom = canonizer.getGraphAtoms();
        List<Integer> productAtomByViewAtom = new ArrayList<>(oldViewAtomByCanonicalAtom.length);
        Map<Integer, Integer> canonicalViewAtomByProductAtom = new HashMap<>();
        for (int canonicalAtom = 0; canonicalAtom < oldViewAtomByCanonicalAtom.length; canonicalAtom++) {
            int productAtom = productByOldViewAtom[oldViewAtomByCanonicalAtom[canonicalAtom]];
            productAtomByViewAtom.add(productAtom);
            canonicalViewAtomByProductAtom.put(productAtom, canonicalAtom);
        }

        List<Integer> openViewAtoms = openScaffoldProductAtoms.stream()
                .map(productAtom -> {
                    Integer viewAtom = canonicalViewAtomByProductAtom.get(productAtom);
                    if (viewAtom == null) {
                        throw new IllegalStateException("Open scaffold atom is absent from view: " + productAtom);
                    }
                    return viewAtom;
                })
                .sorted()
                .toList();
        return new LatentAssemblyDecomposition.View(label, idcode, productAtomByViewAtom, openViewAtoms);
    }

    private static List<Integer> candidateCutBonds(StereoMolecule molecule) {
        List<Integer> bonds = new ArrayList<>();
        molecule.ensureHelperArrays(Molecule.cHelperRings);
        for (int bond = 0; bond < molecule.getBonds(); bond++) {
            int atom1 = molecule.getBondAtom(0, bond);
            int atom2 = molecule.getBondAtom(1, bond);
            if (molecule.getAtomicNo(atom1) <= 1 || molecule.getAtomicNo(atom2) <= 1) {
                continue;
            }
            if (molecule.isRingBond(bond) || molecule.isAromaticBond(bond)
                    || molecule.getBondOrder(bond) != 1) {
                continue;
            }
            bonds.add(bond);
        }
        return bonds;
    }

    private static List<BitSet> connectedComponents(StereoMolecule molecule, BitSet cutBonds) {
        List<BitSet> components = new ArrayList<>();
        BitSet unvisited = allAtoms(molecule);
        while (!unvisited.isEmpty()) {
            int start = unvisited.nextSetBit(0);
            unvisited.clear(start);
            BitSet component = new BitSet(molecule.getAtoms());
            component.set(start);
            ArrayDeque<Integer> queue = new ArrayDeque<>();
            queue.add(start);
            while (!queue.isEmpty()) {
                int atom = queue.removeFirst();
                for (int i = 0; i < molecule.getConnAtoms(atom); i++) {
                    if (cutBonds.get(molecule.getConnBond(atom, i))) {
                        continue;
                    }
                    int neighbor = molecule.getConnAtom(atom, i);
                    if (neighbor >= molecule.getAtoms() || !unvisited.get(neighbor)) {
                        continue;
                    }
                    unvisited.clear(neighbor);
                    component.set(neighbor);
                    queue.add(neighbor);
                }
            }
            components.add(component);
        }
        return components;
    }

    private static BitSet allAtoms(StereoMolecule molecule) {
        BitSet atoms = new BitSet(molecule.getAtoms());
        atoms.set(0, molecule.getAtoms());
        return atoms;
    }

    private static int nonHydrogenAtomCount(StereoMolecule molecule, BitSet atoms) {
        int count = 0;
        for (int atom = atoms.nextSetBit(0); atom >= 0; atom = atoms.nextSetBit(atom + 1)) {
            if (molecule.getAtomicNo(atom) > 1) {
                count++;
            }
        }
        return count;
    }

    private static int ringAtomCount(StereoMolecule molecule, BitSet atoms) {
        int count = 0;
        for (int atom = atoms.nextSetBit(0); atom >= 0; atom = atoms.nextSetBit(atom + 1)) {
            if (molecule.getAtomicNo(atom) > 1 && molecule.isRingAtom(atom)) {
                count++;
            }
        }
        return count;
    }

    private static boolean containsOnlyCommonDrugElements(StereoMolecule molecule) {
        for (int atom = 0; atom < molecule.getAtoms(); atom++) {
            int atomicNo = molecule.getAtomicNo(atom);
            if (atomicNo != 1 && atomicNo != 5 && atomicNo != 6 && atomicNo != 7
                    && atomicNo != 8 && atomicNo != 9 && atomicNo != 14 && atomicNo != 15
                    && atomicNo != 16 && atomicNo != 17 && atomicNo != 34
                    && atomicNo != 35 && atomicNo != 53) {
                return false;
            }
        }
        return true;
    }

    private static BufferedReader newReader(Path path) throws IOException {
        InputStream input = Files.newInputStream(path);
        if (path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz")) {
            input = new GZIPInputStream(input);
        }
        return new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
    }

    @FunctionalInterface
    public interface DecompositionSink {
        void accept(LatentAssemblyDecomposition decomposition) throws IOException;
    }

    public record Result(List<LatentAssemblyDecomposition> decompositions, Stats stats) {
        public Result {
            decompositions = List.copyOf(decompositions);
            stats = Objects.requireNonNull(stats, "stats");
        }
    }

    public record Stats(
            long inputRows,
            long skippedRows,
            long parseFailures,
            long parsedMolecules,
            long rejectedProductsTooSmall,
            long rejectedProductsTooLarge,
            long rejectedDisconnectedProducts,
            long rejectedUncommonElements,
            long eligibleProducts,
            long productsWithDecompositions,
            long acceptedTwoArmCandidates,
            long acceptedThreeArmCandidates,
            long emittedDecompositions
    ) {}

    private static final class MutableStats {
        long inputRows;
        long skippedRows;
        long parseFailures;
        long parsedMolecules;
        long rejectedProductsTooSmall;
        long rejectedProductsTooLarge;
        long rejectedDisconnectedProducts;
        long rejectedUncommonElements;
        long eligibleProducts;
        long productsWithDecompositions;
        long acceptedTwoArmCandidates;
        long acceptedThreeArmCandidates;
        long emittedDecompositions;

        Stats freeze() {
            return new Stats(inputRows, skippedRows, parseFailures, parsedMolecules,
                    rejectedProductsTooSmall, rejectedProductsTooLarge,
                    rejectedDisconnectedProducts, rejectedUncommonElements,
                    eligibleProducts, productsWithDecompositions,
                    acceptedTwoArmCandidates, acceptedThreeArmCandidates,
                    emittedDecompositions);
        }
    }

    private record CandidateArm(BitSet atoms, int cutBond, int scaffoldAtom, int armAtom) {
        CandidateArm {
            atoms = (BitSet) atoms.clone();
        }
    }

    private record Candidate(
            BitSet scaffold,
            List<CandidateArm> arms,
            int scaffoldNonHydrogenAtoms,
            int scaffoldRingAtoms,
            List<Integer> armNonHydrogenAtoms,
            int minimumAttachmentDistance,
            double scaffoldFraction
    ) {
        static final Comparator<Candidate> ORDER = (left, right) -> {
            int comparison = Integer.compare(right.minimumArmSize(), left.minimumArmSize());
            if (comparison != 0) return comparison;
            comparison = Integer.compare(right.minimumAttachmentDistance, left.minimumAttachmentDistance);
            if (comparison != 0) return comparison;
            comparison = Integer.compare(left.armSizeImbalance(), right.armSizeImbalance());
            if (comparison != 0) return comparison;
            comparison = Double.compare(Math.abs(left.scaffoldFraction - 0.5),
                    Math.abs(right.scaffoldFraction - 0.5));
            if (comparison != 0) return comparison;
            return left.key().compareTo(right.key());
        };

        Candidate {
            scaffold = (BitSet) scaffold.clone();
            arms = List.copyOf(arms);
            armNonHydrogenAtoms = List.copyOf(armNonHydrogenAtoms);
        }

        int minimumArmSize() {
            return armNonHydrogenAtoms.stream().mapToInt(Integer::intValue).min().orElse(0);
        }

        int armSizeImbalance() {
            int min = armNonHydrogenAtoms.stream().mapToInt(Integer::intValue).min().orElse(0);
            int max = armNonHydrogenAtoms.stream().mapToInt(Integer::intValue).max().orElse(0);
            return max - min;
        }

        String key() {
            StringBuilder key = new StringBuilder(scaffold.toString());
            arms.forEach(arm -> key.append('|').append(arm.atoms()).append('@').append(arm.scaffoldAtom()));
            return key.toString();
        }
    }

    private static final class HeaderIndex {
        private final Map<String, Integer> columns;

        private HeaderIndex(Map<String, Integer> columns) {
            this.columns = columns;
        }

        static HeaderIndex parse(String header) {
            String[] names = header.split("\t", -1);
            Map<String, Integer> columns = new HashMap<>();
            for (int i = 0; i < names.length; i++) {
                columns.put(names[i], i);
            }
            return new HeaderIndex(columns);
        }

        int required(String name) {
            Integer index = columns.get(name);
            if (index == null) {
                throw new IllegalArgumentException("Missing required input column: " + name);
            }
            return index;
        }
    }

    public enum ArmMode {
        TWO,
        THREE,
        BOTH;

        boolean includes(int armCount) {
            return this == BOTH || (this == TWO && armCount == 2) || (this == THREE && armCount == 3);
        }

        public static ArmMode parse(String raw) {
            return switch (raw.toLowerCase(Locale.ROOT)) {
                case "2", "two", "ab" -> TWO;
                case "3", "three", "abc" -> THREE;
                case "both", "2,3", "ab,abc" -> BOTH;
                default -> throw new IllegalArgumentException("Unsupported arm mode: " + raw);
            };
        }
    }

    public static final class Options {
        private final Path input;
        private final String idColumn;
        private final String structureColumn;
        private final int maxMolecules;
        private final int minProductNonHydrogenAtoms;
        private final int maxProductNonHydrogenAtoms;
        private final int minScaffoldNonHydrogenAtoms;
        private final int maxScaffoldNonHydrogenAtoms;
        private final int minScaffoldRingAtoms;
        private final double minScaffoldFraction;
        private final double maxScaffoldFraction;
        private final int minArmNonHydrogenAtoms;
        private final int maxArmNonHydrogenAtoms;
        private final int minAttachmentDistance;
        private final int maxDecompositionsPerMolecule;
        private final boolean onlyCommonDrugElements;
        private final ArmMode armMode;

        private Options(Builder builder) {
            input = Objects.requireNonNull(builder.input, "input");
            idColumn = requireText(builder.idColumn, "idColumn");
            structureColumn = requireText(builder.structureColumn, "structureColumn");
            maxMolecules = requireNonNegative(builder.maxMolecules, "maxMolecules");
            minProductNonHydrogenAtoms = requirePositive(builder.minProductNonHydrogenAtoms,
                    "minProductNonHydrogenAtoms");
            maxProductNonHydrogenAtoms = requirePositive(builder.maxProductNonHydrogenAtoms,
                    "maxProductNonHydrogenAtoms");
            if (maxProductNonHydrogenAtoms > HARD_MAX_PRODUCT_NON_HYDROGEN_ATOMS) {
                throw new IllegalArgumentException("maxProductNonHydrogenAtoms must not exceed the hard limit of "
                        + HARD_MAX_PRODUCT_NON_HYDROGEN_ATOMS);
            }
            if (minProductNonHydrogenAtoms > maxProductNonHydrogenAtoms) {
                throw new IllegalArgumentException("Product atom bounds are inverted");
            }
            minScaffoldNonHydrogenAtoms = requirePositive(builder.minScaffoldNonHydrogenAtoms,
                    "minScaffoldNonHydrogenAtoms");
            maxScaffoldNonHydrogenAtoms = requirePositive(builder.maxScaffoldNonHydrogenAtoms,
                    "maxScaffoldNonHydrogenAtoms");
            if (minScaffoldNonHydrogenAtoms > maxScaffoldNonHydrogenAtoms) {
                throw new IllegalArgumentException("Scaffold atom bounds are inverted");
            }
            minScaffoldRingAtoms = requireNonNegative(builder.minScaffoldRingAtoms, "minScaffoldRingAtoms");
            minScaffoldFraction = requireFraction(builder.minScaffoldFraction, "minScaffoldFraction");
            maxScaffoldFraction = requireFraction(builder.maxScaffoldFraction, "maxScaffoldFraction");
            if (minScaffoldFraction > maxScaffoldFraction) {
                throw new IllegalArgumentException("Scaffold fractions are inverted");
            }
            minArmNonHydrogenAtoms = requirePositive(builder.minArmNonHydrogenAtoms,
                    "minArmNonHydrogenAtoms");
            maxArmNonHydrogenAtoms = requirePositive(builder.maxArmNonHydrogenAtoms,
                    "maxArmNonHydrogenAtoms");
            if (minArmNonHydrogenAtoms > maxArmNonHydrogenAtoms) {
                throw new IllegalArgumentException("Arm atom bounds are inverted");
            }
            minAttachmentDistance = requireNonNegative(builder.minAttachmentDistance, "minAttachmentDistance");
            maxDecompositionsPerMolecule = requirePositive(builder.maxDecompositionsPerMolecule,
                    "maxDecompositionsPerMolecule");
            onlyCommonDrugElements = builder.onlyCommonDrugElements;
            armMode = Objects.requireNonNull(builder.armMode, "armMode");
        }

        public static Builder builder() { return new Builder(); }
        public Path input() { return input; }
        public String idColumn() { return idColumn; }
        public String structureColumn() { return structureColumn; }
        public int maxMolecules() { return maxMolecules; }
        public int minProductNonHydrogenAtoms() { return minProductNonHydrogenAtoms; }
        public int maxProductNonHydrogenAtoms() { return maxProductNonHydrogenAtoms; }
        public int minScaffoldNonHydrogenAtoms() { return minScaffoldNonHydrogenAtoms; }
        public int maxScaffoldNonHydrogenAtoms() { return maxScaffoldNonHydrogenAtoms; }
        public int minScaffoldRingAtoms() { return minScaffoldRingAtoms; }
        public double minScaffoldFraction() { return minScaffoldFraction; }
        public double maxScaffoldFraction() { return maxScaffoldFraction; }
        public int minArmNonHydrogenAtoms() { return minArmNonHydrogenAtoms; }
        public int maxArmNonHydrogenAtoms() { return maxArmNonHydrogenAtoms; }
        public int minAttachmentDistance() { return minAttachmentDistance; }
        public int maxDecompositionsPerMolecule() { return maxDecompositionsPerMolecule; }
        public boolean onlyCommonDrugElements() { return onlyCommonDrugElements; }
        public ArmMode armMode() { return armMode; }

        private static String requireText(String value, String name) {
            if (value == null || value.isBlank()) throw new IllegalArgumentException(name + " is required");
            return value;
        }
        private static int requirePositive(int value, String name) {
            if (value <= 0) throw new IllegalArgumentException(name + " must be positive");
            return value;
        }
        private static int requireNonNegative(int value, String name) {
            if (value < 0) throw new IllegalArgumentException(name + " must not be negative");
            return value;
        }
        private static double requireFraction(double value, String name) {
            if (!Double.isFinite(value) || value < 0.0 || value > 1.0) {
                throw new IllegalArgumentException(name + " must be between zero and one");
            }
            return value;
        }
    }

    public static final class Builder {
        private Path input;
        private String idColumn = "chembl_id";
        private String structureColumn = "idcode";
        private int maxMolecules;
        private int minProductNonHydrogenAtoms = 12;
        private int maxProductNonHydrogenAtoms = HARD_MAX_PRODUCT_NON_HYDROGEN_ATOMS;
        private int minScaffoldNonHydrogenAtoms = 6;
        private int maxScaffoldNonHydrogenAtoms = 24;
        private int minScaffoldRingAtoms = 5;
        private double minScaffoldFraction = 0.30;
        private double maxScaffoldFraction = 0.75;
        private int minArmNonHydrogenAtoms = 3;
        private int maxArmNonHydrogenAtoms = 12;
        private int minAttachmentDistance = 2;
        private int maxDecompositionsPerMolecule = 1;
        private boolean onlyCommonDrugElements = true;
        private ArmMode armMode = ArmMode.TWO;

        public Builder input(Path value) { input = value; return this; }
        public Builder idColumn(String value) { idColumn = value; return this; }
        public Builder structureColumn(String value) { structureColumn = value; return this; }
        public Builder maxMolecules(int value) { maxMolecules = value; return this; }
        public Builder minProductNonHydrogenAtoms(int value) { minProductNonHydrogenAtoms = value; return this; }
        public Builder maxProductNonHydrogenAtoms(int value) { maxProductNonHydrogenAtoms = value; return this; }
        public Builder minScaffoldNonHydrogenAtoms(int value) { minScaffoldNonHydrogenAtoms = value; return this; }
        public Builder maxScaffoldNonHydrogenAtoms(int value) { maxScaffoldNonHydrogenAtoms = value; return this; }
        public Builder minScaffoldRingAtoms(int value) { minScaffoldRingAtoms = value; return this; }
        public Builder minScaffoldFraction(double value) { minScaffoldFraction = value; return this; }
        public Builder maxScaffoldFraction(double value) { maxScaffoldFraction = value; return this; }
        public Builder minArmNonHydrogenAtoms(int value) { minArmNonHydrogenAtoms = value; return this; }
        public Builder maxArmNonHydrogenAtoms(int value) { maxArmNonHydrogenAtoms = value; return this; }
        public Builder minAttachmentDistance(int value) { minAttachmentDistance = value; return this; }
        public Builder maxDecompositionsPerMolecule(int value) { maxDecompositionsPerMolecule = value; return this; }
        public Builder onlyCommonDrugElements(boolean value) { onlyCommonDrugElements = value; return this; }
        public Builder armMode(ArmMode value) { armMode = value; return this; }
        public Options build() { return new Options(this); }
    }
}
