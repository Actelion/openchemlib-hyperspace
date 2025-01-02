package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.DWARFileParser;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class InputStructureParser {

    /**
     * The ending of the input file determines the type of file that we are trying to load:
     * .txt or .idc will be parsed as line-separated idcodess.
     * .dwar will be parsed as dwar file
     *
     * @return
     */
    public static List<StereoMolecule> defaultParseInput(String inputFile) {

        if(inputFile.toLowerCase().endsWith(".dwar")) {
            System.out.println("[INFO] input= "+inputFile+" -> parse dwar");
            return parseDWARStructures(inputFile);
        }
        if( inputFile.toLowerCase().endsWith(".idc") || inputFile.toLowerCase().endsWith(".txt") ) {
            System.out.println("[INFO] input= "+inputFile+" -> parse idcodes");
            return parseIDC(inputFile);
        }
        if( inputFile.toLowerCase().endsWith(".smi") ) {
            System.out.println("[INFO] input= "+inputFile+" -> parse smiles");
            return parseSMILES(inputFile);
        }
        // last resort, also try to parse idcodes..
        return parseIDC(inputFile);
    }

    public static List<StereoMolecule> parseDWARStructures(String inputFile) {
        return parseDWARStructures(inputFile,"Structure");
    }

    /**
     * Tries to find structure column with colName, if no such column can be found it searches for
     * any structure column.
     *
     * @param inputFile
     * @param colName
     * @return
     */
    public static List<StereoMolecule> parseDWARStructures(String inputFile, String colName) {
        System.out.println("[INFO] Parse dwar: "+inputFile);
        List<StereoMolecule> structures = new ArrayList<>();
        DWARFileParser in = new DWARFileParser(inputFile);
        int sidx = in.getSpecialFieldIndex(colName);

        while(in.next()) {
            try {
                if (sidx >= 0) {
                    structures.add(HyperspaceUtils.parseIDCode(in.getSpecialFieldData(sidx)));
                } else {
                    structures.add(in.getMolecule());
                }
            }
            catch(Exception ex) {
                System.out.println("[WARN] error parsing structure");
            }
        }

        System.out.println("[INFO] Done parsing dwar file, parsed "+structures.size()+" molecules");
        return structures;
    }

    public static List<StereoMolecule> parseSMILES(String inputFile) {
        System.out.println("[INFO] Parse smiles: "+inputFile);
        List<StereoMolecule> structures = new ArrayList<>();
        SmilesParser parser = new SmilesParser();
        try(FileReader f_in = new FileReader(inputFile)) {
            BufferedReader in = new BufferedReader(f_in);
            String line = null;
            while( (line=in.readLine()) != null ) {
                try {
                    StereoMolecule mi = new StereoMolecule();
                    parser.parse(mi, line.trim());
                    mi.ensureHelperArrays(Molecule.cHelperCIP);
                    mi.setFragment(true);
                    structures.add(mi);
                }
                catch(Exception ex) {
                    System.out.println("[WARN] Could not parse input line: "+line+" -> skip");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("[INFO] Parse smiles: done, parsed structures: "+structures.size());

        return structures;

    }

    public static List<StereoMolecule> parseIDC(String inputFile) {
        System.out.println("[INFO] Parse idcodes: "+inputFile);
        List<StereoMolecule> structures = new ArrayList<>();
        IDCodeParser icp  = new IDCodeParser();
        try(FileReader f_in = new FileReader(inputFile)) {
            BufferedReader in = new BufferedReader(f_in);
            String line = null;
            while( (line=in.readLine()) != null ) {
                try {
                    StereoMolecule mi = new StereoMolecule();
                    icp.parse(mi, line.trim());
                    mi.ensureHelperArrays(Molecule.cHelperCIP);
                    structures.add(mi);
                }
                catch(Exception ex) {
                    System.out.println("[WARN] Could not parse input line: "+line+" -> skip");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("[INFO] Parse idcodes: done, parsed structures: "+structures.size());

        return structures;
    }
}
