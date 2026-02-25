package com.idorsia.research.chem.hyperspace.cleaning;

import com.actelion.research.chem.StereoMolecule;

import java.util.List;

/**
 * Strategy interface that wraps an external normalization workflow.
 * Implementations are expected to return at least one molecule and to preserve
 * atom ordering (i.e. the n-th atom in the returned molecules corresponds to
 * the n-th atom in the input {@link StereoMolecule}).
 */
@FunctionalInterface
public interface StructureCleaner {

    /**
     * Cleans the provided structure and returns one or more enumerated variants.
     *
     * @param structure molecule to clean (callers will provide a defensive copy)
     * @return list of normalized structures; implementors must keep the atom
     * ordering identical to the input. Returning {@code null} or an empty list
     * signals that the input should be used as-is.
     */
    List<StereoMolecule> cleanStructure(StereoMolecule structure);
}
