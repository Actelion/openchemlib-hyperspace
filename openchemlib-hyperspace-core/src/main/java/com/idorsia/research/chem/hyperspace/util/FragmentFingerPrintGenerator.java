/**
 *
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 */


package com.idorsia.research.chem.hyperspace.util;



import com.actelion.research.chem.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 *
 * NOT THREAD-SAFE!
 *
 * Generates fragment fingerprints.
 *
 * The type of considered fragments as well as the information included
 * in the fingerprints can be configured:
 *
 * 1. considered fragments are selected by three values:
 * - min_size
 * - max_size
 * - only_paths: if true, then only fragments without branching are included
 *
 * 2. information included in the fingerprints:
 *
 * - atom label
 * - bond type (single/double/triple/aromatic)
 *
 * - num molecule neighbors : number of heavy-weight atom neighbors in the molecule
 * - num fragment neighbors : number of heavy-weight atom neighbors in the fragment
 *
 * - molecule valence state (valence state considering the whole molecule, (sp,db,tp,arom) -> (2,4,6,3) )
 * - fragment valence state (valence state considering only fragment-internal bonds, (sp,db,tp,arom) -> (2,4,6,3) )
 *
 * - stereo (labels for atoms: "no sc" , "sc, not specified", "sc-r", "sc-s", "sc-e", "sc-z")
 * - ring property: true if the atom is contained in a ring (in the whole molecule)
 *
 *
 */
public class FragmentFingerPrintGenerator {


    private int mMinSize = 0;
    private int mMaxSize = 8;

    private boolean mIncludeAtomLabel     = true;
    private boolean mIncludeBondType_a      = true;
    private boolean mIncludeBondType_b      = true;

    private boolean mIncludeMoleculeValenceState  = true;
    private boolean mIncludeFragmentValenceState  = true;

    private boolean mIncludeNumMoleculeNeighbors  = true;
    private boolean mIncludeNumFragmentNeighbors  = true;

    private boolean mIncludeStereo = true;
    private boolean mIncludeRingProperty = true;

    private int mNumBits = 2048;

    // currently processed stereo molecule
    private StereoMolecule M;

    // currently processed fragment StereoMolecule
    private StereoMolecule F;

    // map of frag atoms to molecule atoms of currently processed fragment
    private int[] fragToMol;

    // atom labels of the current fragment
    private String[] customFragmentAtomLabels;

    public FragmentFingerPrintGenerator(int min_size, int max_size) {
        this.mMinSize = min_size;
        this.mMaxSize = max_size;
    }


    public void setNumBits(int bits) {
        if(bits%64!=0) {
            throw new Error("Number of bits must be multiple of 64");
        }
        this.mNumBits = bits;
    }

    public long[] getFingerPrint(StereoMolecule m) {
        if(mNumBits%64!=0) {
            throw new Error("Number of bits must be multiple of 64");
        }

        String[] labels = getFingerPrintStrings(m);

        BitSet bs = new BitSet(mNumBits);
        for(int zi=0;zi<labels.length;zi++) {
            int fancy_hash = new Random(labels[zi].hashCode()).nextInt(mNumBits);
            bs.set( fancy_hash );
            //int hash = labels[zi].hashCode();
            //bs.set( fancy_hash % mNumBits );
        }
        return bs.toLongArray();
    }

    public void setIncludeAtomLabel(boolean set) {this.mIncludeAtomLabel=set;}

    /**
     * distinguishes between triple bonds and all other bonds. If F is set, then we distinguish
     * between single/double/triple/delocalized
     * @param set
     */
    public void setIncludeBondType_a(boolean set) {this.mIncludeBondType_a=set;}
    public void setIncludeBondType_b(boolean set) {this.mIncludeBondType_b=set;}

    public void setIncludeNumFragmentNeighbors(boolean set) {this.mIncludeNumFragmentNeighbors=set;}
    public void setIncludeNumMoleculeNeighbors(boolean set) {this.mIncludeNumMoleculeNeighbors=set;}
    public void setIncludeRingProperty(boolean set) {this.mIncludeRingProperty=set;}
    public void setIncludeFragmentValenceState(boolean set) {this.mIncludeFragmentValenceState=set;}
    public void setIncludeMoleculeValenceState(boolean set) {this.mIncludeMoleculeValenceState=set;}



    /**
     * NOTE! Not thread-safe!
     *
     * @param m
     * @return
     */
    public String[] getFingerPrintStrings(StereoMolecule m) {

        Canonizer cm = new Canonizer(m);
        m = cm.getCanMolecule(false);

        this.M = m;
        this.M.ensureHelperArrays(StereoMolecule.cHelperCIP);

        // generate the "molecule" parts of the atom labels:
        this.atomLabelMoleculePart = new String[this.M.getAtoms()];
        for(int zi=0;zi<this.M.getAtoms();zi++) {
            this.atomLabelMoleculePart[zi] = generateAtomLabel_MoleculePart(zi);
        }

        // generate fragments
        FragmentGenerator fg = new FragmentGenerator(m,this.mMinSize,this.mMaxSize);

        String fp_strings[] = new String[fg.getFragments().size()];

        int zfi = 0;
        for( boolean bfi[] : fg.getFragmentsAsBooleanArrays() ) {
            StereoMolecule fm = new StereoMolecule();
            int atom_map[] = new int[bfi.length];
            m.copyMoleculeByAtoms(fm,bfi,true,atom_map);
            fm.ensureHelperArrays(ExtendedMolecule.cHelperCIP);

            this.F = fm;

            int map_frag_to_mol[] = new int[fm.getAtoms()];
            for(int zf=0;zf<atom_map.length;zf++) {
                if(atom_map[zf]>=0) { map_frag_to_mol[ atom_map[zf] ] = zf; }
            }
            this.fragToMol = map_frag_to_mol;


            // 1. compute atom labels
            // 2. attach them as custom atom labels, then canonize
            this.customFragmentAtomLabels = new String[fm.getAtoms()];
            for(int zi=0;zi<this.F.getAtoms();zi++) {
                customFragmentAtomLabels[zi] = this.customFragmentAtomLabels[zi] = generateAtomLabel(zi);
                fm.setAtomCustomLabel(zi,customFragmentAtomLabels[zi]);
            }

            fm.setFragment(false);
            Canonizer c = new Canonizer(fm,Canonizer.ENCODE_ATOM_CUSTOM_LABELS);
            //int graph[] = c.getGraphAtoms();
            int graph[] = c.getGraphAtoms();
            //int ranks[] = c.getSymmetryRank(1);
            fm.setFragment(true);


            // dfs over fragment and generate fingerprint:
            // note: we can visit a vertex a second time
            StringBuilder frag_code = new StringBuilder();

            Map<Integer,Integer> priority = new HashMap<>();
            for(int zi=0;zi<graph.length;zi++) { priority.put(graph[zi],(graph.length+1)-zi );}

            BitSet visited = new BitSet(graph.length);
            int v      = graph[0];
            int v_last = -1;

            dfs_fragment(fm,v,v_last,visited,priority,frag_code);

            fp_strings[zfi] = frag_code.toString();
            zfi++;
        }

        return fp_strings;
    }


    private String[] atomLabelMoleculePart;

    /**
     * Generates the atom label part of the molecule properties.
     *
     * @param molecule_atom : atom in molecule
     * @return
     */
    private String generateAtomLabel_MoleculePart(int molecule_atom) {
        StringBuilder s = new StringBuilder();

        if(mIncludeNumMoleculeNeighbors) {
            s.append("mn");
            s.append( this.M.getNonHydrogenNeighbourCount(molecule_atom) );
        }

        if(mIncludeRingProperty) {
            if(this.M.isRingAtom(molecule_atom)) {
                s.append("R+");
            }
            else {
                s.append("R-");
            }
        }

        if(mIncludeMoleculeValenceState){
            s.append("mv");
            int val = 0;
            for(int zi=0;zi<this.M.getConnAtoms(molecule_atom);zi++) {
                int b = this.M.getConnBond(molecule_atom,zi);
                if(this.M.isDelocalizedBond(b) || this.M.isAromaticBond(b)) {
                    val+=3;
                }
                else {
                    val+= 2*this.M.getBondOrder(b);
                }
            }
            s.append(val);
        }

        return s.toString();
    }

    /**
     * Generates the atom labels used in the Fingerprint, according to the
     * configuration of the FragmentFingerPrintGenerator.
     *
     * NOTE! You must call the function generateAtomLabel_MoleculePart before calling this!
     *
     * Argument is atom number in fragment
     * @param a
     */
    private String generateAtomLabel(int a) {

        StereoMolecule fragment = this.F;

        StringBuilder s = new StringBuilder();
        int molecule_atom = this.fragToMol[a];

        s.append("a");
        if(mIncludeAtomLabel) {
            s.append(fragment.getAtomicNo(a));
        }
        if(mIncludeNumFragmentNeighbors) {
            s.append("fn");
            s.append( fragment.getNonHydrogenNeighbourCount(a) );
        }

        if(mIncludeFragmentValenceState){
            s.append("fv");
            int val = 0;
            for(int zi=0;zi<this.F.getConnAtoms(a);zi++) {
                int b = this.F.getConnBond(a,zi);
                if(this.F.isDelocalizedBond(b) || this.F.isAromaticBond(b)) {
                    val+=3;
                }
                else {
                    val+= 2*this.F.getBondOrder(b);
                }
            }
            s.append(val);

            s.append(atomLabelMoleculePart[molecule_atom]);
        }

        if(mIncludeStereo) {
            // TODO, or omit stereo stuff..
        }

        return s.toString();
    }

    private void append_atom_properties(StereoMolecule fragment, int vertex, StringBuilder s) {

        s.append(this.customFragmentAtomLabels[vertex]);
        // debug:
        //s.append("v("+vertex+")");
    }

    private void append_bond_properties(StereoMolecule fragment, int vertex_last, int vertex, StringBuilder s) {

        int bond = fragment.getBond(vertex_last,vertex);
        if(!mIncludeBondType_a && !mIncludeBondType_b) {
            s.append("-");
            return;
        }
        if(mIncludeBondType_a && !mIncludeBondType_b) {
            if(fragment.getBondOrder(bond)==3){s.append("b3"); }
            else{s.append("b-");}
            return;
        }


        if( fragment.isAromaticBond(bond) || fragment.isDelocalizedBond(bond) ) {
            s.append(":");
            return;
        }
        int bo = fragment.getBondOrder(bond);
        s.append("b"+bo);
    }

    /**
     * Does four things:
     * 1. encode bond (if not first call)
     * 2. checks if we close a cycle (in this case we put a ring-closure identifier and return)
     * 3. encodes the given atom
     * 4. determines order of next vertices to visit
     * 5. for each next vertice:
     *     5.1 continue recursion
     *
     * @param fragment
     * @param vertex
     * @param last_vertex
     * @param visited
     * @param priority
     * @param frag_code
     */
    private void dfs_fragment(StereoMolecule fragment, int vertex, int last_vertex, BitSet visited, Map<Integer,Integer> priority, StringBuilder frag_code) {

        // encode bond (if it is not the first call)
        if(last_vertex>=0) {
            append_bond_properties(fragment,vertex,last_vertex,frag_code);
        }

        // check if we closed a cycle
        // in this case we add the vertex priority to indicate which ring closure happened
        if(visited.get(vertex)) {
            frag_code.append("rc");
            frag_code.append(priority.get(vertex));
            return;
        }

        // ENCODE THE ATOM PROPERTIES:
        append_atom_properties(fragment,vertex,frag_code);

        // update visited
        //BitSet visited_p = (BitSet) visited.clone();
        BitSet visited_p = visited;
        visited_p.set(vertex);

        // now determine next vertices to visit
        PriorityQueue<FVertex> vn = new PriorityQueue<>();
        int neighbors = fragment.getNonHydrogenNeighbourCount(vertex);
        for(int zi=0;zi<neighbors;zi++) {
            int nai = fragment.getConnAtom(vertex,zi);
            if(nai!=last_vertex) {
                vn.add(new FVertex(nai,priority.get(nai)));
            }
        }

        // RECURSION:
        while(!vn.isEmpty()) {
            FVertex vi = vn.poll();
            dfs_fragment(fragment,vi.a2,vertex,visited,priority,frag_code);
        }

    }

    // NOTE: Highest priority is the least element.
    private static final class FVertex implements Comparable<FVertex>{
        public final int a2;
        public final int priority;
        public FVertex(int a2, int priority) {
            this.a2 = a2;
            this.priority = priority;
        }
        @Override
        public int compareTo(FVertex o) {
            return -Integer.compare(this.priority,o.priority);
        }
    }


    public static void main(String args[]) {
        //test_benchmark_0();
        //test_benchmark_1();

        hyper_fp_test_0();
    }

    public static void hyper_fp_test_0() {
        String smi_a = "Cc(cc1)ccc1S(C(CCCC12)C2NC=CC1c(cc1)ccc1F)(=O)=O";
        String smi_b = "Cc(cc1)ccc1S(C(CCCC1C2c(cc3)ccc3F)C1Nc1c2c(CCC2)c2s1)(=O)=O";

        SmilesParser sp = new SmilesParser();
        StereoMolecule sa = new StereoMolecule();
        StereoMolecule sb = new StereoMolecule();

        sa.setFragment(true);
        sb.setFragment(true);


        try {
            sp.parse(sa, smi_a);
            sp.parse(sb, smi_b);
        } catch (Exception e) {
            e.printStackTrace();
        }

        FragmentFingerPrintGenerator fp = FragmentFingerPrintGenerator.getHyperFP(1024);

        BitSet fp_a = BitSet.valueOf( fp.getFingerPrint(sa) );
        BitSet fp_b = BitSet.valueOf( fp.getFingerPrint(sb) );
        BitSet b2 = (BitSet) fp_b.clone();
        b2.or(fp_a);
        if (!b2.equals(fp_b)) {
            System.out.println("wrong?");
        }
        else {
            System.out.println("ok");
        }
    }

    public static void test_benchmark_1() {
        String filename_a = "C:\\datasets\\idcodes\\rand_600_from_osiris.txt";

        IDCodeParser parser = new IDCodeParser();


        List<StereoMolecule> list = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filename_a))) {
            String line;
            while ((line = br.readLine()) != null) {
                StereoMolecule mol = new StereoMolecule();
                parser.parse(mol, line.trim());
                list.add(mol);
            }
        } catch (Exception e) {
            System.out.println(e);
        }
        list = list.subList(0,500);

        //FragmentFingerPrintGenerator ffpg = new FragmentFingerPrintGenerator(1,8);
        FragmentFingerPrintGenerator ffpg = getHyperFP(1024);

        long ts_a = System.currentTimeMillis();

        for(int zi=0;zi<list.size();zi++) {
            String si[] = ffpg.getFingerPrintStrings(list.get(zi));
            System.out.println(si.length);
        }

        long ts_b = System.currentTimeMillis();
        System.out.println("Time= "+(ts_b-ts_a) +" , d.h. per molecule: " + ( (ts_b-ts_a) /(1.0*list.size()) ) );
    }


    public static void test_benchmark_0() {
        String filename_a = "C:\\datasets\\idcodes\\rand_600_from_osiris.txt";

        IDCodeParser parser = new IDCodeParser();


        List<StereoMolecule> list = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filename_a))) {
            String line;
            while ((line = br.readLine()) != null) {
                StereoMolecule mol = new StereoMolecule();
                parser.parse(mol, line.trim());
                list.add(mol);
            }
        } catch (Exception e) {
            System.out.println(e);
        }
        list = list.subList(0,500);

        FragmentFingerPrintGenerator ffpg = new FragmentFingerPrintGenerator(1,8);

        long ts_a = System.currentTimeMillis();

        for(int zi=0;zi<list.size();zi++) {
            String si[] = ffpg.getFingerPrintStrings(list.get(zi));
            System.out.println(si.length);
        }

        long ts_b = System.currentTimeMillis();
        System.out.println("Time= "+(ts_b-ts_a) +" , d.h. per molecule: " + ( (ts_b-ts_a) /(1.0*list.size()) ) );
    }

    public static void main__(String args[]) {

        String smiles_A = "O1CCCCC1N1CCCCC1";
        String smiles_B = "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5";

        String smiles_cyclohexen = "C=1CCCCC=1";

        String smiles_C = "N1CC2CCCC2CC1";

        SmilesParser sp = new SmilesParser();
        StereoMolecule mA = new StereoMolecule();
        StereoMolecule mB = new StereoMolecule();
        StereoMolecule mC = new StereoMolecule();

        StereoMolecule mCyclohexen = new StereoMolecule();

        try {
            sp.parse(mA,smiles_A);
            sp.parse(mB,smiles_B);
            sp.parse(mCyclohexen,smiles_cyclohexen);
            sp.parse(mC,smiles_C);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //FragmentGenerator fg = new FragmentGenerator(mA,0,6);

        FragmentFingerPrintGenerator fpg = new FragmentFingerPrintGenerator(4,4);


        long ts_a = System.currentTimeMillis();
        //String[] sfp = fpg.getFingerPrintStrings(mCyclohexen);
        String[] sfp = fpg.getFingerPrintStrings(mC);
        long ts_b = System.currentTimeMillis();

        System.out.println("Num Fragments: "+sfp.length);
        System.out.println("Num Fragments (unique): "+Arrays.stream(sfp).distinct().count());
        System.out.println("Time: "+(ts_b-ts_a));

        for(Object fsi : Arrays.stream(sfp).distinct().toArray()) {
            System.out.println(fsi);
        }

    }

    public static void main_(String args[]) {

        SmilesParser sp = new SmilesParser();
        StereoMolecule mC = new StereoMolecule();
        String smiles_C = "N1CC2CCCC2CC1";
        try {
            sp.parse(mC,smiles_C);
        } catch (Exception e) {
            e.printStackTrace();
        }

        for(int zb=0;zb<mC.getBonds();zb++) {
            boolean frag[] = new boolean[mC.getAtoms()];
            frag[mC.getBondAtom(0,zb)]=true;
            frag[mC.getBondAtom(1,zb)]=true;
            StereoMolecule fc = new StereoMolecule();
            mC.copyMoleculeByAtoms(fc,frag,true,null);
            System.out.println( fc.getIDCode() );
        }
    }


    public static FragmentFingerPrintGenerator getHyperFP(int bits) {
        FragmentFingerPrintGenerator fp = new FragmentFingerPrintGenerator(0,5);
        fp.setNumBits(bits);
        fp.setIncludeAtomLabel(true);
        fp.setIncludeBondType_a(true);
        fp.setIncludeBondType_b(false);

        fp.setIncludeFragmentValenceState(false);
        fp.setIncludeMoleculeValenceState(false);

        fp.setIncludeRingProperty(false);
        fp.setIncludeNumFragmentNeighbors(false);
        fp.setIncludeNumMoleculeNeighbors(false);

        return fp;
    }

}
