package com.idorsia.research.chem.hyperspace.fragment;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public interface FragmentGenerator {

    //public static boolean keep_fragment_stereo_molecules_in_memory = true;
    public static Map<String,StereoMolecule> cached_fragments = new ConcurrentHashMap<>();

    public List<Frag> createFragments(StereoMolecule mol);

    public static class Frag {
        //public StereoMolecule parent = null;
        //public StereoMolecule frag = null;

        public final boolean atom_map[];

        public final String parent_idcode;
        public final String frag_idcode;

        public final int frag_num_atoms;
        public final int frag_num_bonds;

        public Frag( StereoMolecule parent, StereoMolecule frag, boolean atom_map[] ) {
            //this.parent = parent;
            //this.frag = frag;
            this.parent_idcode = parent.getIDCode();
            this.frag_idcode   = frag.getIDCode();
            this.atom_map = Arrays.copyOf( atom_map , atom_map.length );
            this.frag_num_atoms = frag.getAtoms();
            this.frag_num_bonds = frag.getBonds();

            if(!cached_fragments.containsKey(this.frag_idcode)) {
                cached_fragments.put(this.frag_idcode,frag);
            }


        }

        public boolean equals(Object o) {
            if(!(o instanceof Frag)) {
                return false;
            }
            Frag fo = (Frag) o;
            return  parent_idcode.equals( fo.parent_idcode ) && frag_idcode.equals( fo.frag_idcode );
        }

        public int hashCode() {
            return parent_idcode.hashCode() ^ frag_idcode.hashCode();
        }

        public StereoMolecule getParent() {
            IDCodeParser icp = new IDCodeParser();
            StereoMolecule mi = new StereoMolecule();
            icp.parse(mi,this.parent_idcode);
            return mi;
        }

        public StereoMolecule getFrag() {
            if(cached_fragments.containsKey(this.frag_idcode)) {
                StereoMolecule mi = new StereoMolecule( cached_fragments.get(this.frag_idcode) );
                mi.ensureHelperArrays(Molecule.cHelperCIP);
                return mi;
            }
            else {
                IDCodeParser icp = new IDCodeParser();
                StereoMolecule mi = new StereoMolecule();
                icp.parse(mi, this.frag_idcode);
                cached_fragments.put(this.frag_idcode,mi);
                return mi;
            }
        }

    }

}
