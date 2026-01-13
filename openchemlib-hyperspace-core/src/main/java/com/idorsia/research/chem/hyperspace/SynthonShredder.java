package com.idorsia.research.chem.hyperspace;


import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

/**
 *
 * Note: generating all splits with all connector permutations is a two-step
 * procedure:
 * 1. use the function trySplit() to generate "preliminary" Splits, where
 *    all connector nodes are just represented as Uranium atoms.
 * 2. use getAllSplitVariations(..) to enumerate all splits with unique
 *    connectors. This works up to 4 connectors!
 *
 */

public class SynthonShredder {

    public static int logLevel_A = 0;

    /**
     * For a given SplitResult, it computes all connector assignments.
     * NOTE! In the SplitResult from the trySplit function, the connectors are just represented as Uranium atoms.
     * This function now generates all Uranium/Neptunium/Plutonium/Americanum connector variants.
     *
     * @param si
     * @return
     */
    public static List<SplitResult> getAllSplitVariations(SplitResult si, int num_connectors) {

        List<SplitResult> results = new ArrayList<>();

        // now we have to consider all possible connector assignments..
        //List<StereoMolecule> variants = getAllConnectorAssignmentsForFragment(frags[0], num_splits);
        List<Integer> connectors = new ArrayList<>();

        //for(int zi=0;zi<num_splits;zi++){ connectors.add(zi); }
        for(int zi=0;zi<num_connectors;zi++){ connectors.add(zi); }

        List<List<Integer>> connector_permutations = SynthonSpace.all_permutations(connectors);

        // Now: blue gets first, red second, orange third. 0 means U, 1 means Np, 2 means Pu
        int connector_atoms[] = new int[]{92,93,94,95};

        for ( List<Integer> connector_perm :  connector_permutations ) {

            String cutset_perm_hash    = connector_perm.stream().map( xi -> ""+xi ).collect(Collectors.joining(";"));

            StereoMolecule frag_copy[] = new StereoMolecule[si.fragments.length];
            // adjust the fragments accordingly:
            //for (StereoMolecule frag_i : frags) {
            List<Set<Integer>> connector_config = new ArrayList<>();

            for(int fi=0;fi<si.fragments.length;fi++){
                StereoMolecule frag_i = si.fragments[fi];
            //for (StereoMolecule frag_i : si.fragments) {
                Set<Integer> connectors_i = new HashSet<>();
                for (int zi = 0; zi < frag_i.getAtoms(); zi++) {
                    if (frag_i.getAtomicNo(zi) == 92 || frag_i.getAtomicNo(zi) == 93 || frag_i.getAtomicNo(zi) == 94 || frag_i.getAtomicNo(zi) == 95 ) {
                        if (frag_i.getAtomColor(zi) == StereoMolecule.cAtomColorBlue) {
                            frag_i.setAtomicNo(zi, connector_atoms[connector_perm.get(0)]);
                            connectors_i.add(connector_atoms[connector_perm.get(0)]);
                        } else if (frag_i.getAtomColor(zi) == StereoMolecule.cAtomColorRed) {
                            frag_i.setAtomicNo(zi, connector_atoms[connector_perm.get(1)]);
                            connectors_i.add(connector_atoms[connector_perm.get(1)]);
                        } else if (frag_i.getAtomColor(zi) == StereoMolecule.cAtomColorOrange) {
                            frag_i.setAtomicNo(zi, connector_atoms[connector_perm.get(2)]);
                            connectors_i.add(connector_atoms[connector_perm.get(2)]);
                        } else if (frag_i.getAtomColor(zi) == StereoMolecule.cAtomColorMagenta) {
                            frag_i.setAtomicNo(zi, connector_atoms[connector_perm.get(3)]);
                            connectors_i.add(connector_atoms[connector_perm.get(3)]);
                        } else {
                            System.out.println("something wrong?..");
                        }
                    }
                }
                connector_config.add(connectors_i);

                StereoMolecule sm_cfi = new StereoMolecule(frag_i);
                sm_cfi.ensureHelperArrays(Molecule.cHelperCIP);
                frag_copy[fi] = sm_cfi;
            }
            results.add(new SplitResult(frag_copy,connector_config,si.cutset_hash+"_p_"+cutset_perm_hash,si.connector_positions));
        }

        return results;
    }

    /**
     * Generates a split with a single connector configuration.
     * All connectors are denoted by U, the matching connectors are colored with the same atomic color.
     * A split result returned by this function can be expanded into the splits for all connector configs
     * via the getAllSplitVariations function.
     *
     * NOTE: we return null in case that the hit is not valid
     *
     *
     * @param mol
     * @param bond_cutset
     * @param max_fragments if the split results in more than max_fragments, we do not consider it.
     * @return
     */
    public static SplitResult trySplit( StereoMolecule mol, int bond_cutset[] , int max_fragments) {

        //int BOND_HANDLING_MODE = 1; // this was the original mode, that just added a bond of the right order.
        int BOND_HANDLING_MODE = 2; // this is the mode using the copyBond operation;

        int si[]       = bond_cutset;
        int num_splits = bond_cutset.length;

        StereoMolecule mi = new StereoMolecule(mol);
        //mi.ensureHelperArrays(StereoMolecule.cHelperCIP);
        mi.ensureHelperArrays(StereoMolecule.cHelperNeighbours);

//        // 1. check if bonds from cutset are adjacent..
//        boolean has_adjacent = false;
//        for(int zi=0;zi<bond_cutset.length-1;zi++) {
//            for(int zj=zi+1;zj<bond_cutset.length;zj++) {
//                int ba_a = mi.getBondAtom(0,bond_cutset[zi]); int ba_b = mi.getBondAtom(1,bond_cutset[zi]);
//                int bb_a = mi.getBondAtom(0,bond_cutset[zj]); int bb_b = mi.getBondAtom(1,bond_cutset[zj]);
//                if(ba_a==bb_a || ba_a==bb_b || ba_b == bb_a || ba_b == bb_b) {
//                    has_adjacent = true; break;
//                }
//            }
//        }
//        if(has_adjacent) {
//            return null;
//        }

        Map<Integer, int[]> bond_atom_indices = new HashMap<>();
        Map<Integer, Integer> bond_types = new HashMap<>();
        Map<Integer, Integer> bond_qfs = new HashMap<>();
        for (int zi = 0; zi < si.length; zi++) {
            int ai[] = new int[]{mi.getBondAtom(0, si[zi]), mi.getBondAtom(1, si[zi])};
            bond_atom_indices.put(zi, ai);
            bond_types.put(zi, mi.getBondType(si[zi]));
            bond_qfs.put(zi,mi.getBondQueryFeatures(zi));
        }
        // delete bonds:
        for (int zi = 0; zi < si.length; zi++) {
            mi.markBondForDeletion(si[zi]);
        }
        int new_atom_pos[] = mi.deleteMarkedAtomsAndBonds();

        for (int zi = 0; zi < si.length; zi++) {
            //int ba = mol.getBondAtom(0,si[zi]);
            //int bb = mol.getBondAtom(1,si[zi]);
            int ba = new_atom_pos[bond_atom_indices.get(zi)[0]];
            int bb = new_atom_pos[bond_atom_indices.get(zi)[1]];

            //int bond_order           = mol.getBondOrder(si[zi]);
            //boolean bond_delocalized = mol.isAromaticBond(si[zi]) || mol.isDelocalizedBond(si[zi]);

            int a_new = -1;
            int b_new = -1;

            if(BOND_HANDLING_MODE == 1 ) {
                int bond_type = bond_types.get(zi);

                a_new = mi.addAtom(92);
                b_new = mi.addAtom(92);

                mi.addBond(ba, a_new, bond_type);
                mi.addBond(bb, b_new, bond_type);
            }
            else if (BOND_HANDLING_MODE == 2) {
                // uses the copy bond option
                a_new = mi.addAtom(92);
                b_new = mi.addAtom(92);

                mol.copyBond(mi,si[zi],-1,-1,ba,a_new,true);
                mol.copyBond(mi,si[zi],-1,-1,bb,b_new,true);
            }
            else if (BOND_HANDLING_MODE == 3){
                // handle supported bond query features separately
//                mi.copyBond
//                int bond_queries_raw = bond_qfs.get(zi);
//                if( (bond_queries_raw | Molecule.cBondQF)
            }

//                // color the uranium atoms, to check if they go to different fragments:
            if (zi == 0) {
                mi.setAtomColor(a_new, StereoMolecule.cAtomColorBlue);
                mi.setAtomColor(b_new, StereoMolecule.cAtomColorBlue);
            }
            if (zi == 1) {
                mi.setAtomColor(a_new, StereoMolecule.cAtomColorRed);
                mi.setAtomColor(b_new, StereoMolecule.cAtomColorRed);
            }
            if (zi == 2) {
                mi.setAtomColor(a_new, StereoMolecule.cAtomColorOrange);
                mi.setAtomColor(b_new, StereoMolecule.cAtomColorOrange);
            }
            if (zi == 3) {
                mi.setAtomColor(a_new, StereoMolecule.cAtomColorMagenta);
                mi.setAtomColor(b_new, StereoMolecule.cAtomColorMagenta);
            }
        }
//            try {
//                mi.validate();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
        //mi.ensureHelperArrays(StereoMolecule.cHelperCIP);
        mi.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
        StereoMolecule frags[] = mi.getFragments();

        //System.out.println("Fragments: ");
        for (StereoMolecule fi : frags) {
            fi.ensureHelperArrays(StereoMolecule.cHelperCIP);
            // System.out.print(fi.getAtoms()+" ");
        }

        if (num_splits > 0 && frags.length == 1) {
            //continue;
            return null;
        }

        if (max_fragments < frags.length) {
            //continue;
            return null;
        }

        // we must ensure that all same color uranium atom pairs are on different fragments:
        //Set<Set<Integer>> connector_config = new HashSet<>();
        List<int[]> connector_positions = new ArrayList<>();
        boolean split_unsuccessful = false;
        for (StereoMolecule fi : frags) {

            int connector_positions_i[] = new int[num_splits];
            Arrays.fill(connector_positions_i,-1);

            int cb = 0;
            int cr = 0;
            int co = 0;
            int cm = 0;
            for (int zi = 0; zi < fi.getAtoms(); zi++) {
                if (fi.getAtomColor(zi) == StereoMolecule.cAtomColorBlue) {
                    connector_positions_i[0] = zi;
                    cb++;
                }
                if (fi.getAtomColor(zi) == StereoMolecule.cAtomColorRed) {
                    connector_positions_i[1] = zi;
                    cr++;
                }
                if (fi.getAtomColor(zi) == StereoMolecule.cAtomColorOrange) {
                    connector_positions_i[2] = zi;
                    co++;
                }
                if (fi.getAtomColor(zi) == StereoMolecule.cAtomColorMagenta) {
                    connector_positions_i[3] = zi;
                    cm++;
                }

//                    if(fi.getAtomicNo(zi)==92){ cb++; }
//                    if(fi.getAtomicNo(zi)==93){ cb++; }
//                    if(fi.getAtomicNo(zi)==94){ cb++; }
            }
            if (cb >= 2 || cr >= 2 || co >= 2 || cm >= 2) {
                split_unsuccessful = true;
                //if(this.logLevel>0){ System.out.println("cb="+cb+" cr="+cr+" co="+co);}
            }
            else {
                connector_positions.add(connector_positions_i);
                //System.out.println("cb="+cb+" cr="+cr+" co="+co);
//                Set<Integer> connectors_i = new HashSet<>();
//                if(cb>=1){connectors_i.add(92);}
//                if(cr>=1){connectors_i.add(93);}
//                if(co>=1){connectors_i.add(94);}
            }
        }
        if (split_unsuccessful) {
            //if(this.logLevel>0){ System.out.println("Split not successful!");}
            return null;
        } else {
            if (logLevel_A > 1) {
                System.out.println("Split OK!");
                if(logLevel_A>2) {
                    System.out.println(Arrays.stream(frags).map(mxi -> HyperspaceUtils.idcodeToSmiles(mxi.getIDCode())).reduce((x, y) -> x + "." + y));
                    System.out.println("ok");
                }
            }
        }


//            for(int zi=0;zi<si.length;zi++) {
//                int col = 0;
//                if(zi==0){col=StereoMolecule.cAtomColorBlue;    }
//                if(zi==1){col=StereoMolecule.cAtomColorRed;     }
//                if(zi==2){col=StereoMolecule.cAtomColorOrange); }
//                int pa=-1; int pb=-1;
//                for(int zj=0;zj<mi.getAtoms();zj++) {
//                    if(mi.getAtomicNo(zj)==92){
//                        if(mi.getAtomColor(zj)==col) {
//                            if(pa>=0){pb=zj;}
//                            else{pa=zj;}
//                        }
//                    }
//                }
//                // check distance:
//                mi.getFrag
//            }

        if (logLevel_A > 0) {
            System.out.println("split: " + num_splits + " -> frags: " + frags.length);
            if(logLevel_A > 2) {
                for (int zi = 0; zi < frags.length; zi++) {
                    System.out.println(HyperspaceUtils.idcodeToSmiles(frags[zi].getIDCode()));
                    if (logLevel_A>1) {
                        if (HyperspaceUtils.idcodeToSmiles(frags[zi].getIDCode()).equals("Cc1ccc(C(c(cccc2)c2C([U])=O)[U])cc1")) {
                            System.out.println("ok");
//                            boolean ok_1 = bst.testSubset(getFP(frags[0]), new BitSetTree.Node[1]);
//                            frags[0] = (new Canonizer(frags[0])).getCanMolecule();
//                            boolean ok_2 = bst.testSubset(getFP(frags[0]), new BitSetTree.Node[1]);
//                            frags[0].setFragment(true);
//                            boolean ok_3 = bst.testSubset(getFP(frags[0]), new BitSetTree.Node[1]);
//                            System.out.println("hmm..");
                        }
                    }
                }
            }
        }

        // We must ensure helper arrays for fragments! This enables the correct sort.
        Map<StereoMolecule,Integer> frag_order_before_sort = new HashMap<>();
        for(int zi=0;zi<frags.length;zi++){
            frags[zi].ensureHelperArrays(Molecule.cHelperNeighbours);
            frag_order_before_sort.put(frags[zi],zi);
        }

        // find the largest fragment: (AND WE USE THE IDCODE LEXICOGRAPHIC ORDER AS TIEBREAKER, TO ENSURE CANONICAL ORDER!)
        Arrays.sort(frags, (a, b) -> (-Integer.compare(a.getBonds(), b.getBonds()) != 0) ? (-Integer.compare(a.getBonds(), b.getBonds())) : a.getIDCode().compareTo(b.getIDCode()));

        // order the connector positions according to the sorting of the fragments:
        List<int[]> sorted_connector_pos = new ArrayList<>();
        for(StereoMolecule fi : frags) {
            sorted_connector_pos.add( connector_positions.get(frag_order_before_sort.get(fi)));
        }

        // does not yet have a connector config (this will be computed only in the second step..)
        String hash_bond_cutset =  "cs_" + Arrays.stream(bond_cutset).mapToObj( xi -> ""+xi ).collect(Collectors.joining(","));
        return new SplitResult(frags,null,hash_bond_cutset,sorted_connector_pos);
    }

    /**
     * NOTE! Objects of this class are used for different purpose in different stages of the algorithm!
     *       However, in all stages, we expect the fragments to be in their canonical order (size, biggest first,
     *       then lexicographic IDCode as tiebreaker)
     *
     *
     * The objects generated in the trySplit(..) function contain colored Urarnium atoms (blue/red/orange) to indicate
     * the connector pairs!
     *
     * The objects created in the getAllSplitVariations(..) then contain (all possible variations) of U/Np/Pu atoms
     * to indicate the connector pairs! These fragments can then be compared against the library synthon reaction fragments!
     *
     * connector_positions: the i'th array contains the connector positions for the i'th fragment. If the
     * i'th fragment contains connector i, then the i'th entry in the array contains the position of
     * the connector atom, otherwise -1.
     *
     */
    public final static class SplitResult implements Serializable {
        public final StereoMolecule fragments[];
        public final List<Integer> connector_counts; // counts of connectors, sorted DESCENDING.
        public final String connector_config;
        public final String cutset_hash;
        public final List<int[]> connector_positions;
        public SplitResult(StereoMolecule fragments[], List<Set<Integer>> connector_config_pre, String cutset_hash, List<int[]> connector_positions) {
            List<BitSet> connector_config = new ArrayList<>();
            if(connector_config_pre!=null) {
                for (Set<Integer> si : connector_config_pre) {
                    BitSet bsi = new BitSet();
                    si.stream().forEach(zi -> bsi.set( (zi>=92)?zi-92:zi ));
                    connector_config.add(bsi);
                }
            }
            this.connector_config = connector_config==null?null: SynthonSpace.encodeConnectorConfig(connector_config);
            this.connector_counts = new ArrayList<>();
            for( int cpi[] : connector_positions) {
                int cnt_a = 0; for(int zi=0;zi<cpi.length;zi++) { if(cpi[zi]>=0) { cnt_a++; } }
                this.connector_counts.add(cnt_a);
            }
            this.connector_counts.sort( (x,y) -> -Integer.compare(x,y) );

            this.fragments = fragments;
            this.cutset_hash = cutset_hash;
            this.connector_positions = connector_positions;
        }

        /**
         * This method computes all possible splits with unique connector atoms.
         * This is computed based on the connector_positions data.
         *
         * @return
         */
        public List<StereoMolecule[]> getAllSplitsWithUniqueConnectors() {
            int num_connectors = this.connector_positions.iterator().next().length;
            List<Integer> connectors_to_use = new ArrayList<>();
            for(int zi=92;zi<92+num_connectors;zi++) { connectors_to_use.add(zi); }
            return this.getAllSplitsWithUniqueConnectors(connectors_to_use);
        }

        public List<StereoMolecule[]> getAllSplitsWithUniqueConnectors(List<Integer> connectors_to_use) {
            int num_connectors = this.connector_positions.iterator().next().length;
            List<StereoMolecule[]> results = new ArrayList<>();
            List<List<Integer>> all_conni_perms = SynthonSpace.all_permutations(connectors_to_use);
            for(List<Integer> pi : all_conni_perms) {
                StereoMolecule[] si = new StereoMolecule[this.fragments.length];
                for(int zi=0;zi<this.fragments.length;zi++) {
                    StereoMolecule fi = new StereoMolecule(this.fragments[zi]);
                    fi.ensureHelperArrays(Molecule.cHelperNeighbours);
                    for( int zj=0;zj < num_connectors; zj++) {
                        int cpij = this.connector_positions.get(zi)[zj];
                        if(cpij>=0) {
                            fi.setAtomicNo( cpij , pi.get(zj) );
                        }
                    }
                    fi.ensureHelperArrays(Molecule.cHelperCIP); // hmm.. maybe omit? :)
                    si[zi] = fi;
                }
                results.add(si);
            }
            return results;
        }

        public String toString() {
            String smiles =  "Smiles= "+ Arrays.stream(this.fragments).map( fi -> HyperspaceUtils.idcodeToSmiles(fi.getIDCode()) ).collect(Collectors.joining(".")) ;
            String idc    =  "idcode= "+Arrays.stream(this.fragments).map( fi -> fi.getIDCode() ).collect(Collectors.joining(" ::: "));
            return smiles+" ;; "+idc;
        }
    }




    public static List<SplitResult> computeAllSplitResults(StereoMolecule mi, int num_splits, int max_fragments) {
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        int nb = mi.getBonds();

        List<int[]> combi_list = CombinationGenerator.getAllOutOf(nb, num_splits);
        // !! returns null if b > a..
        if(combi_list==null) { return new ArrayList<>(); }
        List<int[]> splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());

        List<SplitResult> split_results = new ArrayList<>();
        for(int[] split_pattern : splits) {
            SynthonShredder.SplitResult split_result = SynthonShredder.trySplit(mi,split_pattern,max_fragments);
            if(split_result!=null) {
                split_results.add(split_result);
            }
        }

        return split_results;
    }



    public static void main(String args[]) {
        run_test_01();
    }

    private static void run_test_01() {
        SmilesParser spa = new SmilesParser();
        StereoMolecule ma = new StereoMolecule();

        try {
            spa.parse(ma,"CCn1cc(C(=O)O)c(=O)c2cc(F)c(N3CC[N+]4(CN5C(C)=CC(C)=C(C(=O)O)C5=N4)C(C)C3)c(F)c21"); // from chembl
        } catch (Exception e) {
            e.printStackTrace();
        }

        int max_splits = 5;
        int max_fragments = 6;


//        long ts_2_a = System.currentTimeMillis();
//        List<Pair<SynthonShredder2.BinarySplitResult,SynthonShredder2.BinarySplitResult>> all_splits = new ArrayList<>();
//        List<SynthonShredder2.BinarySplitResult> splits_binary = SynthonShredder2.computeAllBinarySplits(ma,num_splits);
//        //expand to full, i.e. take all combinations that sum up to max num_splits splits and that do not
//        //share any edges
//        Set<Set<Integer>> reported_cutsets = new HashSet<>();
//        Map<Set<Integer>,SplitResult> split_results_via_binary_splits = new HashMap<>();
//
//        for(int za=0;za<splits_binary.size();za++) {
//            for(int zb=0;zb<splits_binary.size();zb++) {
//                // we do this asymmetrically, i.e. we check if sb fits into a partition of sa
//                SynthonShredder2.BinarySplitResult sa = splits_binary.get(za);
//                SynthonShredder2.BinarySplitResult sb = splits_binary.get(zb);
//
//                if(sa.cutset.size()+sb.cutset.size() > num_splits) {
//                    continue;
//                }
//
//                if(sa.bonds_in_splits[0].containsAll( sb.cutset )) {
//                    all_splits.add( Pair.of(sa , sb) );
//                    Set<Integer> cutset_i = new HashSet<>(); cutset_i.addAll(sa.cutset); cutset_i.addAll(sb.cutset);
//                    reported_cutsets.add(cutset_i);
//                }
//                else if(sa.bonds_in_splits[1].containsAll( sb.cutset )) {
//                    all_splits.add( Pair.of(sa , sb) );
//                    Set<Integer> cutset_i = new HashSet<>(); cutset_i.addAll(sa.cutset); cutset_i.addAll(sb.cutset);
//                    reported_cutsets.add(cutset_i);
//                }
//            }
//        }
//        // compute cutsets:
//        for(Set<Integer> cutset_i : reported_cutsets) {
//            SplitResult sri = SynthonShredder.trySplit(ma,cutset_i.stream().mapToInt(vi->vi.intValue()).toArray(),max_fragments);
//            if(sri!=null) {
//                split_results_via_binary_splits.put(cutset_i,sri);
//            }
//            else {
//                System.out.println("mkay..");
//            }
//        }
//
//
//        long ts_2_b = System.currentTimeMillis();

        long ts_2_a = System.currentTimeMillis();
        Map<Set<Integer>,SplitResult> split_results_via_binary_splits = SynthonShredder2.computeAllValidSplits(ma,max_splits,max_fragments);
        long ts_2_b = System.currentTimeMillis();

        System.out.println("Time= "+(ts_2_b-ts_2_a)+" ms");

        long ts_a = System.currentTimeMillis();

        List<int[]> combi_list = new ArrayList<>();
        for(int zs=1;zs<=max_splits;zs++) {
            combi_list.addAll( CombinationGenerator.getAllOutOf(ma.getBonds(), zs) );
        }

        List<SplitResult> good_splits = new ArrayList<>();
        Map<Set<Integer>,SplitResult> split_results_conventional = new HashMap<>();

        // !! returns null if b > a..
        //List<int[]> splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());
        List<int[]> splits = combi_list;

        for(int[] si : splits) {
            SynthonShredder.SplitResult split_result = SynthonShredder.trySplit(ma,si,max_fragments);
            Set<Integer> sis = Arrays.stream(si).mapToObj(ii->ii).collect(Collectors.toSet());
            if( split_result != null) {
                //System.out.println("ok: "+ split_result.cutset_hash);
                good_splits.add(split_result);
                split_results_conventional.put(sis,split_result);
            }
            else {
                //System.out.println("not_ok: ");
            }
        }
        long ts_b = System.currentTimeMillis();
        System.out.println("time= "+ (ts_b-ts_a) + " ms" );
        System.out.println("Bye!");


        // print results:
        Set<Set<Integer>> combined = new HashSet<>();
        combined.addAll(split_results_conventional.keySet());
        combined.addAll(split_results_via_binary_splits.keySet());

        for(Set<Integer> si : combined) {
            System.out.println();
            System.out.println("si= "+si.toString());
            boolean in_a = split_results_via_binary_splits.containsKey(si);
            boolean in_b = split_results_conventional.containsKey(si);
            System.out.println("A: "+in_a+"  B: "+in_b);
            if(!in_a && in_b) {
                System.out.println( split_results_conventional.get(si).toString() );
            }
            if( in_a && !in_b) {
                System.out.println( split_results_via_binary_splits.get(si).toString() );
            }
        }

        System.out.println("ok!");

    }


}
