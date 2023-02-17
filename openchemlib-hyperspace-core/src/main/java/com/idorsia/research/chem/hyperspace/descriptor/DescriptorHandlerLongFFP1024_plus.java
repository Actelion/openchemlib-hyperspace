package com.idorsia.research.chem.hyperspace.descriptor;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.*;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * This class augments either the PFP or FFP descriptor by additional bits that encode the
 * part of the molecule around the connector points. This results in much faster runtime for
 * substructure search, as it helps to eliminate non-matching targets efficiently.
 *
 */

public class DescriptorHandlerLongFFP1024_plus extends AbstractDescriptorHandlerLongFP<StereoMolecule> implements Serializable {

    public static final String VERSION = SSSearcherWithIndex.cIndexVersion;
    private static DescriptorHandlerLongFFP1024_plus sDefaultInstance;
    private static final int sLongCount = (SSSearcherWithIndex.getNoOfKeys() + 63) / 64;

    public static DescriptorHandlerLongFFP1024_plus getDefaultInstance() {
        synchronized(DescriptorHandlerLongFFP1024_plus.class) {
            if (sDefaultInstance == null)
                sDefaultInstance = new DescriptorHandlerLongFFP1024_plus("ffp");
        }
        return sDefaultInstance;
    }

    /**
     * mode can be: "ffp" or "pfp"
     *
     * @param mode
     */
    public DescriptorHandlerLongFFP1024_plus(String mode) {
        this.setMode(mode);
    }

    public DescriptorInfo getInfo() {
        return new DescriptorInfo("FFP1024_plus_"+this.mMode,"FFP1024_plus_"+this.mMode,
                "-1",0,true,false,false,true);
        //throw new Error("not implemented");
        //return DescriptorConstants.DESCRIPTOR_FFP512;
    }

    public String getVersion() {
        return VERSION;
    }

    public long[] decode(String s) {
        throw new Error("not implemented");
//        long[] descriptor = (s != null && s.length() == 128) ?
//                SSSearcherWithIndex.getLongIndexFromHexString(s) : super.decode(s);
//        return (descriptor != null && descriptor.length == sLongCount) ? descriptor : null;
    }

    public long[] decode(byte[] bytes) {
        throw new Error("not implemented");
//        long[] descriptor = (bytes != null && bytes.length == 128) ?
//                SSSearcherWithIndex.getLongIndexFromHexString(bytes) : super.decode(bytes);
//        return (descriptor != null && descriptor.length == sLongCount) ? descriptor : null;
    }

    private FingerPrintGenerator fpg = new FingerPrintGenerator();//new FingerPrintGenerator(false,true);


    private String mMode = "ffp";

    public void setMode(String mode) {
        if(! (mode.equals("ffp") || mode.equals("pfp") ) ) {
            throw new Error("Invalid mode: "+mode);
        }
        this.mMode = mode;
    }

    public String getMode() {
        return mMode;
    }

    public long[] createDescriptor(StereoMolecule mol) {

        long[] descriptor_part_standard = null;
        if(mMode.equals("ffp")) {
            descriptor_part_standard = new SSSearcherWithIndex().createLongIndex(mol);
        }
        else if(mMode.equals("pfp")) {
            descriptor_part_standard = fpg.getFingerprint(mol, 512).toLongArray();
        }


        // create path fp from all 92/93/94:
        List<String> path_features = new ArrayList<>();
        for(int zi=0;zi<mol.getAtoms();zi++) {
            if(mol.getAtomicNo(zi)==92){
                boolean visited[] = new boolean[mol.getAtoms()];
                path_features.add("xB92_");
                createPathFingerprintFeature_df_recursive( mol ,  -1, zi,10, "xB92_",path_features , visited);
            }
            if(mol.getAtomicNo(zi)==93){
                boolean visited[] = new boolean[mol.getAtoms()];
                path_features.add("yR93_");
                createPathFingerprintFeature_df_recursive( mol ,  -1, zi,10, "yR93_",path_features , visited );
            }
            if(mol.getAtomicNo(zi)==94){
                boolean visited[] = new boolean[mol.getAtoms()];
                path_features.add("zO94_");
                createPathFingerprintFeature_df_recursive( mol ,  -1, zi,10, "zO94_",path_features , visited);
            }
        }
        // hash the feature strings..
        BitSet bsi = new BitSet(512);
        for(String pi : path_features) {
            int index = new java.util.Random(pi.hashCode()).nextInt(512);
            bsi.set(index);
        }
        long descriptor_part2[] = bsi.toLongArray();

        long descriptor_complete[] = new long[16];
        for(int zi=0;zi<8;zi++){ if(descriptor_part_standard.length>zi){descriptor_complete[zi]   = descriptor_part_standard[zi];} }
        for(int zi=0;zi<descriptor_part2.length;zi++){ descriptor_complete[8+zi]             = descriptor_part2[zi]; }

        return (descriptor_part_standard == null) ? FAILED_OBJECT : descriptor_complete;
    }


    public static List<String> createPathFingerprintFeatures(StereoMolecule m, int pos, int length) {
        boolean visited[] = new boolean[m.getAtoms()];
        ArrayList<String> features = new ArrayList<>();
        createPathFingerprintFeature_df_recursive(m,-1,pos,length,"<"+pos+">",features,visited);
        return features;
    }

    private static void createPathFingerprintFeature_df_recursive(StereoMolecule m, int last_pos, int pos,
                                                                  int remaining_steps, String current_path,
                                                                  List<String> features, boolean visited[]) {
        if(remaining_steps==0){return;};

        String next_path = current_path;
        if(last_pos>=0) {
            int bond_type = m.getBondType( m.getBond(last_pos,pos));
            next_path += (bond_type==3)?"=3=":"--";
        }

        next_path += "A"+m.getAtomicNo(pos);
        features.add(next_path);
        visited[pos] = true;


        for( int zi=0;zi<m.getConnAtoms(pos);zi++) {
            boolean visited_next[] = Arrays.copyOf(visited,visited.length);
            int next_pos = m.getConnAtom(pos,zi);
            if(visited[next_pos]){continue;}
            createPathFingerprintFeature_df_recursive(m,pos,next_pos,remaining_steps-1,next_path,features,visited_next);
        }
    }


    public DescriptorHandler<long[], StereoMolecule> getThreadSafeCopy() {
        return this;
    }


    public static void main(String args[]) {
        test_01();
    }

    public static void test_01() {
        StereoMolecule ma = new StereoMolecule();
        StereoMolecule mb = new StereoMolecule();

        SmilesParser sp = new SmilesParser();

//        String smi_a = "NC(C(Cc1ccccc1)N(N=N[U])[Np])=O";
//        String smi_b = "[U]/N=N\\[Np]";

//        String smi_a = "[Np]CCC1CCCCC1";
//        String smi_b = "CCC(CC(CC[Np])C1)CC1N";

        String smi_a = "[Np]CCC(CCCC1)C1[Pu]";
        String smi_b = "CNCCC(CC(CO)CC1CC[Np])C1[Pu]";



        try {
            sp.parse(ma,smi_a);
            sp.parse(mb,smi_b);
        } catch (Exception e) {
            e.printStackTrace();
        }

        DescriptorHandlerLongFFP1024_plus dh = DescriptorHandlerLongFFP1024_plus.getDefaultInstance();

        long[] bits_a = dh.createDescriptor(ma);
        long[] bits_b = dh.createDescriptor(mb);

        BitSet bsa = BitSet.valueOf(bits_a);
        BitSet bsb = BitSet.valueOf(bits_a);

        //check that a is subset of b:
        BitSet b2 = (BitSet) bsb.clone();
        b2.or(bsa);
        if(b2.equals(bsb)) {
            System.out.println("ok, is subset!");
        }
        else {
            System.out.println("hmm.. is NO subset!");
        }

        System.out.println(bsa);
        System.out.println(bsb);

        System.out.println("ok");
    }

}
