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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Random;
import java.util.concurrent.Callable;
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
        if (source.getDownsamplingAlgorithm().isPresent() || source.getDownsamplingRequest().isPresent()) {
            builder.withDownsamplingMetadata(source.getDownsamplingAlgorithm().orElse(null),
                    source.getDownsamplingRequest().orElse(null));
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
        Map<String, RawSynthon> cleanedById = new HashMap<>();
        Map<String, Map<String, String>> cleanedAttributes = new LinkedHashMap<>();

        Map<Integer, List<SynthonFragment>> fragmentsBySet = decodeFragments(data.getRawFragmentSets());
        Map<String, Map<String, String>> originalAttributes = data.getFragmentAttributes();

        for (Map.Entry<Integer, List<SynthonFragment>> entry : fragmentsBySet.entrySet()) {
            List<RawSynthon> cleanedSet = cleanSynthonSet(reactionId,
                    entry.getKey(),
                    entry.getValue(),
                    fragmentsBySet,
                    cleanedById,
                    originalAttributes,
                    cleanedAttributes,
                    seedGenerator,
                    executor);
            cleanedSets.put(entry.getKey(), cleanedSet);
            System.out.printf("Cleaned reaction %s set %d: %d synthons -> %d cleaned variants%n",
                    reactionId, entry.getKey(), entry.getValue().size(), cleanedSet.size());
        }

        cleanedSets.forEach((idx, list) -> builder.addRawFragments(reactionId, idx, list));
        rebuildDownsampledSets(builder, reactionId, data.getRawDownsampledSets(), cleanedById);
        if (!data.getExampleScaffolds().isEmpty()) {
            builder.addExampleScaffolds(reactionId, data.getExampleScaffolds());
        }
        data.getPartialAssemblies().forEach((missingIdx, assemblies) -> builder.addPartialAssemblies(reactionId, missingIdx, assemblies));
        if (!data.getRepresentativeCompounds().isEmpty()) {
            builder.addRepresentativeCompounds(reactionId, data.getRepresentativeCompounds());
        }
        data.getDescriptors().forEach((key, value) -> builder.addReactionDescriptor(reactionId, key, value));
        cleanedAttributes.forEach((fragmentId, attrs) ->
                attrs.forEach((key, value) -> builder.addFragmentAttribute(reactionId, fragmentId, key, value)));
    }

    private void rebuildDownsampledSets(RawSynthonSpace.Builder builder,
                                        String reactionId,
                                        Map<Integer, List<RawSynthon>> originalDownsampled,
                                        Map<String, RawSynthon> cleanedById) {
        originalDownsampled.forEach((idx, fragments) -> {
            List<RawSynthon> rebuilt = new ArrayList<>();
            for (RawSynthon fragment : fragments) {
                RawSynthon replacement = cleanedById.get(fragment.getFragmentId());
                if (replacement != null) {
                    rebuilt.add(replacement);
                }
            }
            if (!rebuilt.isEmpty()) {
                builder.addRawDownsampledFragments(reactionId, idx, rebuilt);
            }
        });
    }

    private Map<Integer, List<SynthonFragment>> decodeFragments(Map<Integer, List<RawSynthon>> rawSets) {
        IDCodeParser parser = new IDCodeParser();
        Map<Integer, List<SynthonFragment>> fragmentsBySet = new LinkedHashMap<>();
        rawSets.forEach((idx, list) -> {
            List<SynthonFragment> fragments = new ArrayList<>(list.size());
            for (RawSynthon raw : list) {
                StereoMolecule mol = parser.getCompactMolecule(raw.getIdcode());
                mol.ensureHelperArrays(Molecule.cHelperCIP);
                fragments.add(new SynthonFragment(raw, mol));
            }
            fragmentsBySet.put(idx, fragments);
        });
        return fragmentsBySet;
    }

    private List<RawSynthon> cleanSynthonSet(String reactionId,
                                             int fragmentIndex,
                                             List<SynthonFragment> fragments,
                                             Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                             Map<String, RawSynthon> cleanedById,
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
                cleanedById.put(repl.getFragmentId(), repl);
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
        StructureCleaner cleaner = cleanerProvider.acquire();
        try {
            List<StereoMolecule> cleanedSynthonMolecules = new ArrayList<>();
            int assemblies = options.getAssembliesPerSynthon();
            for (int run = 0; run < assemblies; run++) {
                SynthonAssemblyResult assembly = assembleExample(fragment, fragmentsBySet, fragmentIndex, random);
                int[] atomMap = assembly.atomMaps.get(0);
                StereoMolecule cleanerInput = new StereoMolecule(assembly.assembled);
                cleanerInput.ensureHelperArrays(Molecule.cHelperCIP);
                List<StereoMolecule> cleanedAssemblies = cleaner.cleanStructure(cleanerInput);
                if (cleanedAssemblies == null || cleanedAssemblies.isEmpty()) {
                    cleanedAssemblies = Collections.singletonList(assembly.assembled);
                }
                for (StereoMolecule cleanedAssembly : cleanedAssemblies) {
                    cleanedAssembly.ensureHelperArrays(Molecule.cHelperCIP);
                    cleanedSynthonMolecules.add(projectSynthon(fragment.molecule, cleanedAssembly, atomMap));
                }
            }
            if (cleanedSynthonMolecules.isEmpty()) {
                cleanedSynthonMolecules.add(new StereoMolecule(fragment.molecule));
            }
            List<RawSynthon> replacements = toRawSynthons(reactionId, fragmentIndex, fragment.rawSynthon.getFragmentId(), cleanedSynthonMolecules);
            return new CleanTaskResult(fragment, replacements);
        } catch (RuntimeException ex) {
            throw new IllegalStateException("Cleaning failed for fragment " + fragment.rawSynthon.getFragmentId(), ex);
        } finally {
            cleanerProvider.release(cleaner);
        }
    }

    private SynthonAssemblyResult assembleExample(SynthonFragment target,
                                                  Map<Integer, List<SynthonFragment>> fragmentsBySet,
                                                  int targetIndex,
                                                  Random random) {
        List<StereoMolecule> parts = new ArrayList<>();
        List<int[]> atomMaps = new ArrayList<>();
        parts.add(new StereoMolecule(target.molecule));
        fragmentsBySet.forEach((idx, list) -> {
            if (idx == targetIndex || list.isEmpty()) {
                return;
            }
            SynthonFragment partner = list.size() == 1
                    ? list.get(0)
                    : list.get(random.nextInt(list.size()));
            parts.add(new StereoMolecule(partner.molecule));
        });
        StereoMolecule assembled = SynthonAssembler.assembleSynthons_faster(parts, atomMaps);
        assembled.ensureHelperArrays(Molecule.cHelperCIP);
        return new SynthonAssemblyResult(assembled, atomMaps);
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

    private static final class SynthonFragment {
        private final RawSynthon rawSynthon;
        private final StereoMolecule molecule;

        private SynthonFragment(RawSynthon rawSynthon, StereoMolecule molecule) {
            this.rawSynthon = rawSynthon;
            this.molecule = molecule;
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
}
