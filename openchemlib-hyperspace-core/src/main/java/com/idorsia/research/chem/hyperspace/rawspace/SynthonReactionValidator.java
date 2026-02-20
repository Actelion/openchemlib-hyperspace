package com.idorsia.research.chem.hyperspace.rawspace;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Validates individual synthon reactions to ensure connector consistency.
 */
public final class SynthonReactionValidator {

    private SynthonReactionValidator() {
    }

    public static void validate(Map<Integer, List<Object>> molecules) throws Exception {
        int[] connectorCounts = new int[200];
        Map<Integer, Integer> connectorBondTypes = new HashMap<>();

        IDCodeParser parser = new IDCodeParser();

        for (Integer fragmentIdx : molecules.keySet()) {
            Set<Integer> connectorSet = null;

            Object first = molecules.get(fragmentIdx).iterator().next();
            StereoMolecule template;
            if (first instanceof StereoMolecule) {
                template = (StereoMolecule) first;
            } else if (first instanceof String) {
                template = new StereoMolecule();
                parser.parse(template, (String) first);
            } else {
                throw new IllegalArgumentException("Unsupported template object: " + first);
            }
            connectorSet = SynthonSpace.computeConnectorSet(template);
            for (Integer connector : connectorSet) {
                connectorCounts[connector] += 1;
            }

            for (Object entry : molecules.get(fragmentIdx)) {
                StereoMolecule molecule = null;
                if (entry == null) {
                    System.out.println("Null molecule object supplied -> skip");
                    continue;
                } else if (entry instanceof StereoMolecule) {
                    molecule = (StereoMolecule) entry;
                } else if (entry instanceof String) {
                    molecule = new StereoMolecule();
                    parser.parse(molecule, (String) entry);
                } else {
                    System.out.println("Incompatible molecule object supplied -> skip");
                    continue;
                }

                Set<Integer> connectors = SynthonSpace.computeConnectorSet(molecule);
                if (!connectorSet.equals(connectors)) {
                    throw new Exception("Connector Set " + fragmentIdx + " has inconsistent connectors");
                }
                Map<Integer, Integer> connectorPositions = SynthonSpace.computeConnectorPositions(molecule);
                for (Integer connector : connectorPositions.keySet()) {
                    int atomIndex = connectorPositions.get(connector);
                    if (molecule.getConnAtoms(atomIndex) != 1) {
                        throw new Exception("Connector with degree != 1");
                    }
                    int bondOrder = molecule.getBondOrder(molecule.getBond(atomIndex, molecule.getConnAtom(atomIndex, 0)));

                    if (!connectorBondTypes.containsKey(connector)) {
                        connectorBondTypes.put(connector, bondOrder);
                    } else if (bondOrder != connectorBondTypes.get(connector)) {
                        throw new Exception("Inconsistent bond order detected! " + bondOrder
                                + " != " + connectorBondTypes.get(connector));
                    }
                }
            }
        }

        if (connectorBondTypes.isEmpty()) {
            throw new Exception("No connectors detected");
        }

        List<Integer> connectors = connectorBondTypes.keySet().stream().sorted().collect(Collectors.toList());
        if (connectors.get(0) != 92) {
            throw new Exception("Smallest connector is not 92");
        }
        for (int idx = 0; idx < connectors.size(); idx++) {
            if (connectors.get(idx) != 92 + idx) {
                throw new Exception("Connectors are not ascending from 92: " + connectors);
            }
        }

        for (int connector : connectors) {
            if (connectorCounts[connector] != 2) {
                throw new Exception("Found " + connectorCounts[connector] + " connectors of type "
                        + connector + " instead of 2");
            }
        }
    }
}
