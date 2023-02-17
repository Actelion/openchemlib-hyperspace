package com.idorsia.research.chem.hyperspace.fragment;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class DefaultFragmentGenerator implements FragmentGenerator {
    private static final double CORRECTION_FACTOR = 0.7;

    private static final byte[] FAILED_OBJECT = new byte[0];
    //private static final int MAX_SPHERE_COUNT = 5;
    //private static final int EXACT_SPHERE_COUNT = 4;
    //private static final int SKELETON_SPHERE_COUNT = 5;
    //private static final int HASH_BITS = 10;
    //private static final int HASH_INIT = 13;
    //private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);

    int sphere_count = 6;

    public DefaultFragmentGenerator() {
    }

    public DefaultFragmentGenerator(int sphere_count) {
        this.sphere_count = sphere_count;
    }

    public List<Frag> createFragments(StereoMolecule mol) {
        return createFragments(mol,this.sphere_count);
    }

    /**
     * This descriptor requires proper up/down bonds, because it encodes stereo parities.
     * If a passed molecule is generated from idcode parsing, make sure that coordinates
     * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with
     * the respective option.
     */
    public static List<Frag> createFragments(StereoMolecule mol, int sphere_count) {

        int MAX_SPHERE_COUNT   = sphere_count;
        int EXACT_SPHERE_COUNT = sphere_count;

        if (mol == null)
            return null;

        //List<StereoMolecule> skeleton_frags = new ArrayList<>();
        List<Frag> skeleton_frags = new ArrayList<>();

        mol.ensureHelperArrays(Molecule.cHelperRings);
        StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

        //byte[] descriptor = new byte[DESCRIPTOR_SIZE];

//System.out.println("descriptor skeleton spheres:");
        int[] atomList = new int[mol.getAtoms()];
        boolean[] atomMask = new boolean[mol.getAtoms()];
        for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
            if (rootAtom != 0)
                Arrays.fill(atomMask, false);

            int min = 0;
            int max = 0;

            for (int sphere=0; sphere<MAX_SPHERE_COUNT && max<mol.getAtoms(); sphere++) {
                if (max == 0) {
                    atomList[0] = rootAtom;
                    atomMask[rootAtom] = true;
                    max = 1;
                }
                else {
                    int newMax = max;
                    for (int i=min; i<max; i++) {
                        int atom = atomList[i];
                        for (int j=0; j<mol.getConnAtoms(atom); j++) {
                            int connAtom = mol.getConnAtom(atom, j);
                            if (!atomMask[connAtom]) {
                                atomMask[connAtom] = true;
                                atomList[newMax++] = connAtom;
                            }
                        }
                    }
                    min = max;
                    max = newMax;
                }

                mol.copyMoleculeByAtoms(fragment, atomMask, true, null);

                // take fragment as it is
                if (sphere < EXACT_SPHERE_COUNT) {
                    //String idcode = new Canonizer(fragment).getIDCode();
                    //skeleton_frags.add( new StereoMolecule(fragment) );
                    skeleton_frags.add( new Frag(mol, new StereoMolecule(fragment) , atomMask) );

//                    int h = BurtleHasher.hashlittle(idcode, HASH_INIT);
//                    h = (h & BurtleHasher.hashmask(HASH_BITS));
//                    if (descriptor[h] < DescriptorEncoder.MAX_COUNT_VALUE)
//                        descriptor[h]++;
//System.out.println("atom:"+rootAtom+"\tfragment\tradius:"+sphere+"\thash:"+h+"\t"+idcode);
                }

//                // take atomic no reduced fragment skeleton also
//                if (sphere < SKELETON_SPHERE_COUNT) {
//                    for (int atom=0; atom<fragment.getAllAtoms(); atom++) {
//                        fragment.setAtomicNo(atom, 6);
//                    }
//                    skeleton_frags.add(fragment);
////                    String idcode = new Canonizer(fragment).getIDCode();
////                    int h = BurtleHasher.hashlittle(idcode, HASH_INIT);
////                    h = (h & BurtleHasher.hashmask(HASH_BITS));
////                    if (descriptor[h] < DescriptorEncoder.MAX_COUNT_VALUE)
////                        descriptor[h]++;
////System.out.println("atom:"+rootAtom+"\tskeleton\tradius:"+sphere+"\thash:"+h+"\t"+idcode);
//                }
            }
        }

        return skeleton_frags;
        //return descriptor;
    }


    public static void main(String args[]) {
        String idc_a = "fegi`@DXUkpAYrJJJJXiYIQPqILehpymUSTuUT@C@@@"; // [from chembl]
        IDCodeParser icp = new IDCodeParser();
        StereoMolecule m = new StereoMolecule();
        icp.parse(m,idc_a);

        for(int zi=1;zi<40;zi++) {
            //List<String> skelspheres = DefaultFragmentGenerator.createFragments(m,zi).stream().map(mi->mi.frag.getIDCode()).distinct().collect(Collectors.toList());
            List<String> skelspheres = DefaultFragmentGenerator.createFragments(m,zi).stream().map(mi->mi.frag_idcode).distinct().collect(Collectors.toList());
            System.out.println("Size="+zi+" -> Frags: "+skelspheres.size());
            System.out.println("Frag[idcode]");
            for(String si : skelspheres) {
                System.out.println(si);
            }
        }

    }

}