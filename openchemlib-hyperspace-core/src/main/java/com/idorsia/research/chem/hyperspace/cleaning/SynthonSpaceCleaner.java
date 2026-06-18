package com.idorsia.research.chem.hyperspace.cleaning;

import com.actelion.research.chem.ExtendedMoleculeFunctions;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * High-level workflow that cleans every synthon within a {@link RawSynthonSpace}
 * by assembling temporary example molecules, delegating to a {@link StructureCleaner},
 * and projecting the normalized atom state back onto the individual synthons.
 */
public final class SynthonSpaceCleaner {

    private final StructureCleaner.Provider cleanerProvider;
    private final SynthonCleaningOptions options;

    public SynthonSpaceCleaner(StructureCleaner structureCleaner) {
        this(StructureCleaner.Provider.singleton(structureCleaner), SynthonCleaningOptions.builder().build());
    }

    public SynthonSpaceCleaner(StructureCleaner structureCleaner, SynthonCleaningOptions options) {
        this(StructureCleaner.Provider.singleton(structureCleaner), options);
    }

    public SynthonSpaceCleaner(StructureCleaner.Provider cleanerProvider) {
        this(cleanerProvider, SynthonCleaningOptions.builder().build());
    }

    public SynthonSpaceCleaner(StructureCleaner.Provider cleanerProvider, SynthonCleaningOptions options) {
        this.cleanerProvider = Objects.requireNonNull(cleanerProvider, "structureCleanerProvider");
        this.options = options == null ? SynthonCleaningOptions.builder().build() : options;
    }

    public SynthonCleaningOptions getOptions() {
        return options;
    }

    public RawSynthonSpace clean(RawSynthonSpace source) {
        Objects.requireNonNull(source, "source");
        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(source.getName())
                .version(source.getVersion());
        source.getMetadata().forEach(builder::putMetadata);
        if (options.getMaxSynthonAtoms() > 0) {
            builder.putMetadata("cleaning.maxSynthonAtoms", Integer.toString(options.getMaxSynthonAtoms()));
        }

        ExecutorService executor = Executors.newFixedThreadPool(options.getMaxParallelism());
        try {
            Random seedGenerator = options.newRandom();
            for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : source.getReactions().entrySet()) {
                processReaction(builder, entry.getKey(), entry.getValue(), seedGenerator, executor);
            }
            return builder.build();
        } finally {
            executor.shutdown();
        }
    }

    private void processReaction(RawSynthonSpace.Builder builder,
                                 String reactionId,
                                 RawSynthonSpace.ReactionData data,
                                 Random seedGenerator,
                                 ExecutorService executor) {
        Map<Integer, List<RawSynthon>> cleanedSets = new LinkedHashMap<>();
        Map<String, Map<String, String>> cleanedAttributes = new LinkedHashMap<>();

        Map<Integer, List<SynthonFragment>> fragmentsBySet = filterByMaxSynthonAtoms(reactionId,
                indexFragments(data.getRawFragmentSets()),
                executor);
        Map<String, Map<String, String>> originalAttributes = data.getFragmentAttributes();

        for (Map.Entry<Integer, List<SynthonFragment>> entry : fragmentsBySet.entrySet()) {
            if (entry.getValue().isEmpty()) {
                System.out.printf("Skipping reaction %s because set %d is empty after max-synthon-atoms filtering%n",
                        reactionId, entry.getKey());
                return;
            }
            List<RawSynthon> cleanedSet = cleanSynthonSet(reactionId,
                    entry.getKey(),
                    entry.getValue(),
                    fragmentsBySet,
                    originalAttributes,
                    cleanedAttributes,
                    seedGenerator,
                    executor);
            if (cleanedSet.isEmpty()) {
                System.out.printf("Skipping reaction %s because set %d has no cleaned synthons%n",
                        reactionId, entry.getKey());
                return;
            }
            cleanedSets.put(entry.getKey(), cleanedSet);
            System.out.printf("Cleaned reaction %s set %d: %d synthons -> %d cleaned variants%n",
                    reactionId, entry.getKey(), entry.getValue().size(), cleanedSet.size());
        }

        cleanedSets.forEach((idx, list) -> builder.addRawFragments(reactionId, idx, list));
        if (!data.getExampleScaffolds().isEmpty()) {
            builder.addExampleScaffolds(reactionId, data.getExampleScaffolds());
        }
        data.getPartialAssemblies().forEach((missingIdx, assemblies) -> builder.addPartialAssemblies(reactionId, missingIdx, assemblies));
        if (!data.getRepresentativeCompounds().isEmpty()) {
            builder.addRepresentativeCompounds(reactionId, data.getRepresentativeCompounds());
        }
        data.getDescriptors().forEach((key, value) -> builder.addReactionDescriptor(reactionId, key, value));
        data.getReactionMetadata().forEach((key, value) -> builder.addReactionMetadata(reactionId, key, value));
        cleanedAttributes.forEach((fragmentId, attrs) ->
                attrs.forEach((key, value) -> builder.addFragmentAttribute(reactionId, fragmentId, key, value)));
    }

    private Map<Integer, List<SynthonFragment>> indexFragments(Map<Integer, List<RawSynthon>> rawSets) {
        Map<Integer, List<SynthonFragment>> fragmentsBySet = new LinkedHashMap<>();
        rawSets.forEach((idx, list) -> {
            List<SynthonFragment> fragments = new ArrayList<>(list.size());
            for (RawSynthon raw : list) {
                fragments.add(new SynthonFragment(raw));
            }
            fragmentsBySet.put(idx, fragments);
        });
        return fragmentsBySet;
    }

    private Map<Integer, List<SynthonFragment>> filterByMaxSynthonAtoms(String reactionId,
                                                                       Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                                                       ExecutorService executor) {
        int maxAtoms = options.getMaxSynthonAtoms();
        if (maxAtoms <= 0) {
            return fragmentsBySet;
        }
        Map<Integer, List<SynthonFragment>> filtered = new LinkedHashMap<>();
        fragmentsBySet.forEach((setIndex, fragments) -> {
            FilterTaskResult result = filterSynthonSetByMaxAtoms(reactionId, setIndex, fragments, maxAtoms, executor);
            int skipped = result.skipped;
            if (skipped > 0) {
                System.out.printf("Filtered reaction %s set %d: skipped %d synthons above %d atoms%n",
                        reactionId, setIndex, skipped, maxAtoms);
            }
            filtered.put(setIndex, result.retained);
        });
        return filtered;
    }

    private FilterTaskResult filterSynthonSetByMaxAtoms(String reactionId,
                                                        int setIndex,
                                                        List<SynthonFragment> fragments,
                                                        int maxAtoms,
                                                        ExecutorService executor) {
        if (fragments.isEmpty()) {
            return new FilterTaskResult(Collections.emptyList(), 0);
        }

        int taskCount = Math.min(fragments.size(), options.getMaxParallelism());
        List<Future<FilterTaskResult>> futures = new ArrayList<>(taskCount);
        for (int task = 0; task < taskCount; task++) {
            int from = task * fragments.size() / taskCount;
            int to = (task + 1) * fragments.size() / taskCount;
            futures.add(executor.submit(() -> filterSynthonRangeByMaxAtoms(fragments, from, to, maxAtoms)));
        }

        List<SynthonFragment> retained = new ArrayList<>(fragments.size());
        int skipped = 0;
        for (Future<FilterTaskResult> future : futures) {
            FilterTaskResult result;
            try {
                result = future.get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException("Max-synthon-atoms filtering interrupted", e);
            } catch (ExecutionException e) {
                throw new IllegalStateException("Max-synthon-atoms filtering failed for reaction "
                        + reactionId + " set " + setIndex, e.getCause());
            }
            retained.addAll(result.retained);
            skipped += result.skipped;
        }
        return new FilterTaskResult(retained, skipped);
    }

    private FilterTaskResult filterSynthonRangeByMaxAtoms(List<SynthonFragment> fragments,
                                                          int from,
                                                          int to,
                                                          int maxAtoms) {
        IDCodeParser parser = new IDCodeParser();
        List<SynthonFragment> retained = new ArrayList<>(to - from);
        int skipped = 0;
        for (int idx = from; idx < to; idx++) {
            SynthonFragment fragment = fragments.get(idx);
            StereoMolecule molecule = decodeMolecule(parser, fragment.rawSynthon);
            if (molecule.getAtoms() <= maxAtoms) {
                retained.add(fragment);
            } else {
                skipped++;
            }
        }
        return new FilterTaskResult(retained, skipped);
    }

    private List<RawSynthon> cleanSynthonSet(String reactionId,
                                             int fragmentIndex,
                                             List<SynthonFragment> fragments,
                                             Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                             Map<String, Map<String, String>> originalAttributes,
                                             Map<String, Map<String, String>> cleanedAttributes,
                                             Random seedGenerator,
                                             ExecutorService executor) {
        List<Future<CleanTaskResult>> futures = new ArrayList<>(fragments.size());
        for (SynthonFragment fragment : fragments) {
            long fragmentSeed = seedGenerator.nextLong();
            futures.add(executor.submit(() -> cleanFragment(reactionId,
                    fragmentIndex,
                    fragment,
                    fragmentsBySet,
                    fragmentSeed)));
        }

        List<RawSynthon> cleaned = new ArrayList<>();
        for (int idx = 0; idx < futures.size(); idx++) {
            CleanTaskResult result;
            try {
                result = futures.get(idx).get();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new IllegalStateException("Cleaning interrupted", e);
            } catch (ExecutionException e) {
                throw new IllegalStateException("Cleaning failed for reaction " + reactionId, e.getCause());
            }
            List<RawSynthon> replacements = result.replacements;
            if (replacements.isEmpty()) {
                replacements = Collections.singletonList(result.fragment.rawSynthon);
            }
            cleaned.addAll(replacements);
            Map<String, String> attrs = originalAttributes.get(result.fragment.rawSynthon.getFragmentId());
            for (RawSynthon repl : replacements) {
                if (attrs != null && !attrs.isEmpty()) {
                    cleanedAttributes.put(repl.getFragmentId(), new LinkedHashMap<>(attrs));
                }
            }
        }

        return cleaned;
    }

    private CleanTaskResult cleanFragment(String reactionId,
                                          int fragmentIndex,
                                          SynthonFragment fragment,
                                          Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                          long randomSeed) {
        Random random = new Random(randomSeed);
        IDCodeParser parser = new IDCodeParser();
        StereoMolecule targetMolecule = decodeMolecule(parser, fragment.rawSynthon);
        StructureCleaner cleaner = cleanerProvider.acquire();
        try {
            List<StereoMolecule> cleanedSynthonMolecules = new ArrayList<>();
            int assemblies = options.getAssembliesPerSynthon();
            for (int run = 0; run < assemblies; run++) {
                SynthonAssemblyResult assembly = assembleExample(parser, targetMolecule, fragmentsBySet, fragmentIndex, random);
                int[] atomMap = assembly.atomMaps.get(0);
                StereoMolecule cleanerInput = new StereoMolecule(assembly.assembled);
                cleanerInput.ensureHelperArrays(Molecule.cHelperCIP);
                List<StereoMolecule> cleanedAssemblies = cleaner.cleanStructure(cleanerInput);
                if (cleanedAssemblies == null || cleanedAssemblies.isEmpty()) {
                    cleanedAssemblies = Collections.singletonList(assembly.assembled);
                }
                for (StereoMolecule cleanedAssembly : cleanedAssemblies) {
                    cleanedAssembly.ensureHelperArrays(Molecule.cHelperCIP);
                    cleanedSynthonMolecules.add(projectSynthon(targetMolecule, cleanedAssembly, atomMap));
                }
            }
            if (cleanedSynthonMolecules.isEmpty()) {
                cleanedSynthonMolecules.add(new StereoMolecule(targetMolecule));
            }
            List<RawSynthon> replacements = deduplicateByIdcode(
                    toRawSynthons(reactionId, fragmentIndex, fragment.rawSynthon.getFragmentId(), cleanedSynthonMolecules));
            return new CleanTaskResult(fragment, replacements);
        } catch (RuntimeException ex) {
            throw new IllegalStateException("Cleaning failed for fragment " + fragment.rawSynthon.getFragmentId(), ex);
        } finally {
            cleanerProvider.release(cleaner);
        }
    }

    private SynthonAssemblyResult assembleExample(IDCodeParser parser,
                                                  StereoMolecule targetMolecule,
                                                  Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                                  int targetIndex,
                                                  Random random) {
        List<StereoMolecule> parts = new ArrayList<>();
        List<int[]> atomMaps = new ArrayList<>();
        parts.add(new StereoMolecule(targetMolecule));
        fragmentsBySet.forEach((idx, list) -> {
            if (idx == targetIndex || list.isEmpty()) {
                return;
            }
            SynthonFragment partner = list.size() == 1
                    ? list.get(0)
                    : list.get(random.nextInt(list.size()));
            parts.add(decodeMolecule(parser, partner.rawSynthon));
        });
        StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts, atomMaps);
        assembled.ensureHelperArrays(Molecule.cHelperCIP);
        return new SynthonAssemblyResult(assembled, atomMaps);
    }

    private StereoMolecule decodeMolecule(IDCodeParser parser, RawSynthon rawSynthon) {
        StereoMolecule molecule = parser.getCompactMolecule(rawSynthon.getIdcode());
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule;
    }

    private StereoMolecule projectSynthon(StereoMolecule template,
                                          StereoMolecule cleanedAssembly,
                                          int[] atomMap) {
        StereoMolecule result = new StereoMolecule(template);
        result.ensureHelperArrays(Molecule.cHelperNeighbours);
        for (int atom = 0; atom < result.getAtoms() && atom < atomMap.length; atom++) {
            int assemblyAtom = atomMap[atom];
            if (assemblyAtom < 0 || assemblyAtom >= cleanedAssembly.getAtoms()) {
                continue;
            }
            result.setAtomCharge(atom, cleanedAssembly.getAtomCharge(assemblyAtom));
            result.setAtomParity(atom, cleanedAssembly.getAtomParity(assemblyAtom),
                    cleanedAssembly.isAtomParityPseudo(assemblyAtom));
            result.setAtomConfigurationUnknown(atom, cleanedAssembly.isAtomConfigurationUnknown(assemblyAtom));
            result.setAtomESR(atom, cleanedAssembly.getAtomESRType(assemblyAtom),
                    cleanedAssembly.getAtomESRGroup(assemblyAtom));
            result.setAtomRadical(atom, cleanedAssembly.getAtomRadical(assemblyAtom));
        }

        for (int bond = 0; bond < result.getBonds(); bond++) {
            int atom1 = result.getBondAtom(0, bond);
            int atom2 = result.getBondAtom(1, bond);
            if (atom1 >= atomMap.length || atom2 >= atomMap.length) {
                continue;
            }
            int assemblyAtom1 = atomMap[atom1];
            int assemblyAtom2 = atomMap[atom2];
            if (assemblyAtom1 < 0 || assemblyAtom2 < 0) {
                continue;
            }
            int assemblyBond = ExtendedMoleculeFunctions.getBondNo(cleanedAssembly, assemblyAtom1, assemblyAtom2);
            if (assemblyBond < 0) {
                continue;
            }
            result.setBondType(bond, cleanedAssembly.getBondType(assemblyBond));
            result.setBondParity(bond, cleanedAssembly.getBondParity(assemblyBond),
                    cleanedAssembly.isBondParityPseudo(assemblyBond));
            result.setBondESR(bond, cleanedAssembly.getBondESRType(assemblyBond),
                    cleanedAssembly.getBondESRGroup(assemblyBond));
        }

        result.ensureHelperArrays(Molecule.cHelperCIP);
        return result;
    }

    private List<RawSynthon> toRawSynthons(String reactionId,
                                           int fragmentIndex,
                                           String baseFragmentId,
                                           List<StereoMolecule> molecules) {
        List<RawSynthon> result = new ArrayList<>(molecules.size());
        boolean reuseOriginalId = options.reuseFragmentIdWhenSingleResult() && molecules.size() == 1;
        int duplicateIndex = 0;
        for (StereoMolecule molecule : molecules) {
            String fragmentId;
            if (reuseOriginalId && duplicateIndex == 0) {
                fragmentId = baseFragmentId;
            } else {
                fragmentId = options.formatFragmentId(baseFragmentId, duplicateIndex);
            }
            result.add(RawSynthon.fromMolecule(reactionId, fragmentIndex, fragmentId, molecule));
            duplicateIndex++;
        }
        return result;
    }

    private List<RawSynthon> deduplicateByIdcode(List<RawSynthon> synthons) {
        Map<String, RawSynthon> unique = new LinkedHashMap<>();
        for (RawSynthon synthon : synthons) {
            unique.putIfAbsent(synthon.getIdcode(), synthon);
        }
        return new ArrayList<>(unique.values());
    }

    private static final class SynthonFragment {
        private final RawSynthon rawSynthon;

        private SynthonFragment(RawSynthon rawSynthon) {
            this.rawSynthon = rawSynthon;
        }
    }

    private static final class SynthonAssemblyResult {
        private final StereoMolecule assembled;
        private final List<int[]> atomMaps;

        private SynthonAssemblyResult(StereoMolecule assembled, List<int[]> atomMaps) {
            this.assembled = assembled;
            this.atomMaps = atomMaps;
        }
    }

    private static final class CleanTaskResult {
        private final SynthonFragment fragment;
        private final List<RawSynthon> replacements;

        private CleanTaskResult(SynthonFragment fragment, List<RawSynthon> replacements) {
            this.fragment = fragment;
            this.replacements = replacements;
        }
    }

    private static final class FilterTaskResult {
        private final List<SynthonFragment> retained;
        private final int skipped;

        private FilterTaskResult(List<SynthonFragment> retained, int skipped) {
            this.retained = retained;
            this.skipped = skipped;
        }
    }
}
