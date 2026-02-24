package com.idorsia.research.chem.hyperspace.rawspace;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Objects;

/**
 * Minimal synthon representation used in {@link RawSynthonSpace}.
 */
public final class RawSynthon implements Serializable {

    private final String reactionId;
    private final int fragmentIndex;
    private final String fragmentId;
    private final String idcode;
    private final BitSet connectors;

    public RawSynthon(String reactionId,
                      int fragmentIndex,
                      String fragmentId,
                      String idcode,
                      BitSet connectors) {
        this.reactionId = Objects.requireNonNull(reactionId, "reactionId");
        this.fragmentIndex = fragmentIndex;
        this.fragmentId = Objects.requireNonNull(fragmentId, "fragmentId");
        this.idcode = Objects.requireNonNull(idcode, "idcode");
        this.connectors = connectors == null ? new BitSet() : (BitSet) connectors.clone();
    }

    public String getReactionId() {
        return reactionId;
    }

    public int getFragmentIndex() {
        return fragmentIndex;
    }

    public String getFragmentId() {
        return fragmentId;
    }

    public String getIdcode() {
        return idcode;
    }

    public BitSet getConnectors() {
        return (BitSet) connectors.clone();
    }

    public SynthonSpace.FragId toFragId() {
        return new SynthonSpace.FragId(reactionId,
                fragmentIndex,
                idcode,
                fragmentId,
                (BitSet) connectors.clone(),
                null,
                null);
    }

    public static RawSynthon fromFragId(SynthonSpace.FragId fragId) {
        return new RawSynthon(fragId.rxn_id,
                fragId.frag,
                fragId.fragment_id,
                fragId.idcode,
                fragId.getConnectors());
    }

    public static RawSynthon fromMolecule(String reactionId,
                                          int fragmentIndex,
                                          String fragmentId,
                                          StereoMolecule molecule) {
        StereoMolecule working = new StereoMolecule(molecule);
        working.ensureHelperArrays(Molecule.cHelperNeighbours);
        BitSet connectors = SynthonSpace.computeConnectorBitSet(working);
        return new RawSynthon(reactionId, fragmentIndex, fragmentId, working.getIDCode(), connectors);
    }
}
