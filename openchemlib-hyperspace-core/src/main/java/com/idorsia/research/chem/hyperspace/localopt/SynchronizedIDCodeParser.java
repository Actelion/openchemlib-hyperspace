package com.idorsia.research.chem.hyperspace.localopt;

import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.StereoMolecule;

/**
 * Serialize IDCode parsing through a single parser to avoid non-thread-safe coordinate generation.
 */
final class SynchronizedIDCodeParser {

    private static final IDCodeParserWithoutCoordinateInvention PARSER = new IDCodeParserWithoutCoordinateInvention();

    private SynchronizedIDCodeParser() {}

    static synchronized void parse(StereoMolecule target, String idcode) {
        PARSER.parse(target, idcode);
    }
}
