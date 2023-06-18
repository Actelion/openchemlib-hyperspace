package com.idorsia.research.chem.hyperspace.reparametrization;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.ScaffoldHelper;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import java.util.ArrayList;
import java.util.List;

public class SynthonSpaceScaffold {

    final String assembled_idcode;
    final String assembled_murcko_scaffold_idcode;
    final String assembled_no_exit_vectors_idcode;
    final String s1_idcode;
    final String s2_idcode;
    final String s3_idcode;
    final List<String> s1_variants;
    final List<String> s2_variants;
    final List<String> s3_variants;

    final String rxn;

    public SynthonSpaceScaffold(String assembled_idcode, String assembled_murcko_scaffold_idcode, String assembled_no_exit_vectors_idcode, String s1_idcode, String s2_idcode, String s3_idcode, List<String> s1_variants, List<String> s2_variants, List<String> s3_variants, String rxn) {
        this.assembled_idcode = assembled_idcode;
        this.assembled_murcko_scaffold_idcode = assembled_murcko_scaffold_idcode;
        this.assembled_no_exit_vectors_idcode = assembled_no_exit_vectors_idcode;
        this.s1_idcode = s1_idcode;
        this.s2_idcode = s2_idcode;
        this.s3_idcode = s3_idcode;
        this.s1_variants = s1_variants;
        this.s2_variants = s2_variants;
        this.s3_variants = s3_variants;
        this.rxn = rxn;
    }

    public SynthonSpaceScaffold(String rxn, String assembled_ic, String s1_ic, String s2_ic, String s3_ic,
                          List<String> s1_v, List<String> s2_v, List<String> s3_v) {
        this.rxn = rxn;
        this.assembled_idcode = assembled_ic;
        this.s1_idcode=s1_ic;
        this.s2_idcode=s2_ic;
        this.s3_idcode=s3_ic;
        this.s1_variants = s1_v;
        this.s2_variants = s2_v;
        this.s3_variants = s3_v;

        // compute the remaining values..
        StereoMolecule mi = new StereoMolecule();
        mi = HyperspaceUtils.parseIDCode(assembled_idcode);
        StereoMolecule murcko = new StereoMolecule(mi);
        murcko.ensureHelperArrays(Molecule.cHelperCIP);
        ScaffoldHelper.createMurckoScaffold(murcko,false);
        //StereoMolecule no_exit_vectors = new StereoMolecule(mi);
        //no_exit_vectors.ensureHelperArrays(Molecule.cHelperCIP);
        StereoMolecule no_exit_vectors = HyperspaceUtils.removeConnectorsFromStructure( mi , false);
        this.assembled_murcko_scaffold_idcode = murcko.getIDCode();
        this.assembled_no_exit_vectors_idcode = no_exit_vectors.getIDCode();
    }

    public int countConnectors() {
        IDCodeParser icp = new IDCodeParser();
        StereoMolecule m = new StereoMolecule();
        icp.parse(m,assembled_idcode);
        m.ensureHelperArrays(Molecule.cHelperNeighbours);
        int count_th = 0;
        for(int zi=0;zi<m.getAtoms();zi++) {
            if(m.getAtomicNo(zi)==90){ count_th++;}
        }
        return count_th;
    }

    public String toCSV_short() {
        List<String> parts = new ArrayList<>();
        parts.add(this.assembled_idcode);
        parts.add(this.assembled_murcko_scaffold_idcode);
        parts.add(this.assembled_no_exit_vectors_idcode);
        parts.add(this.s1_idcode);
        parts.add(this.s2_idcode);
        parts.add(this.s3_idcode);
        parts.add(""+this.s1_variants.size());
        parts.add(""+this.s2_variants.size());
        parts.add(""+this.s3_variants.size());
        long ns1 = Math.max(1,this.s1_variants.size());
        long ns2 = Math.max(1,this.s2_variants.size());
        long ns3 = Math.max(1,this.s3_variants.size());
        parts.add(""+ (ns1*ns2*ns3) );
        parts.add(""+this.rxn);
        //parts.add(""+countConnectors());
        return String.join(",",parts);
    }

}