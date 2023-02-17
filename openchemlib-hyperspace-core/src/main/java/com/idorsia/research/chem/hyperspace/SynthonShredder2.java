package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.util.*;

public class SynthonShredder2 {


    public static Map<Integer, List<SynthonShredder.SplitResult>> computeSplitResults(StereoMolecule m, int max_splits, int max_fragments) {

        m.ensureHelperArrays(Molecule.cHelperNeighbours);

        //m.cop

        // create all 2 fragment splits:
        for(int zs=0;zs<max_splits;zs++) {



        }
        return null;
    }

    public static class BinarySplitResult {
        Set<Integer> cutset;
        //SynthonShredder.SplitResult split_result;
        StereoMolecule fragments[];
        Set<Integer>[] bonds_in_splits; // NOTE: there is no particular order of sets a and b..
        public BinarySplitResult(Set<Integer> cutset,StereoMolecule fragments[], Set<Integer>[] bonds_in_splits) {
            this.cutset = cutset;
            this.fragments = fragments;
            this.bonds_in_splits = bonds_in_splits;
        }
        public Set<Integer> getBondCutSet() {
            return this.cutset;
        }

        public Set<Integer>[] getBondsInSplits() { return bonds_in_splits;}
        public Set<Integer> getBondsInA() { return bonds_in_splits[0];}
        public Set<Integer> getBondsInB() { return bonds_in_splits[1];}

        BitSet[] bonds_bitset = null;
        public BitSet[] getBondSetsAsBitSets(int bitset_length) {
            if(bonds_bitset==null) {
                bonds_bitset = new BitSet[2];;
                bonds_bitset[0] = new BitSet(bitset_length);
                bonds_bitset[1] = new BitSet(bitset_length);
                for(int bii : bonds_in_splits[0]){bonds_bitset[0].set(bii,true);}
                for(int bii : bonds_in_splits[1]){bonds_bitset[1].set(bii,true);}
            }
            return bonds_bitset;
        }
    }

    public static List<BinarySplitResult> computeAllBinarySplits(StereoMolecule m_in, int max_splits) {
        m_in.ensureHelperArrays(Molecule.cHelperCIP);

        List<BinarySplitResult> all_split_sets = new ArrayList<>();

        Map<Integer,Set<Integer>> n_connectedness = new HashMap<>();
        for(int zi=1;zi<=max_splits;zi++) { n_connectedness.put(zi,new HashSet<>()); }

        List<Triple<Set<Integer>,StereoMolecule,int[]>> current_partial_splits = new ArrayList<>();
        // always contains the map from original m_in indeces to current molecule indeces
        int bond_order[] = new int[m_in.getBonds()];
        for(int zi=0;zi<bond_order.length;zi++) { bond_order[zi] = zi; }
        current_partial_splits.add(Triple.of(new HashSet<>(),m_in,bond_order));

        Set<Integer> all_bonds = new HashSet<>();
        for(int zi=0;zi<m_in.getBonds();zi++) { all_bonds.add(zi); }

        // create all 2 fragment splits:
        for(int zs=1;zs<=max_splits;zs++) {

            // this is always in "m_in" bond indeces !!
            Set<Integer> taken_bonds = new HashSet<>();
            for(int zxs=1;zxs<zs;zxs++) { taken_bonds.addAll( n_connectedness.get(zxs) ); }

            List<Triple<Set<Integer>,StereoMolecule,int[]>> next_partial_splits = new ArrayList<>();
            for(Triple<Set<Integer>,StereoMolecule,int[]> partial_split : current_partial_splits) {
                StereoMolecule m = partial_split.getMiddle();
                int bond_map[] = partial_split.getRight();
                //for (int za = 0; za < m.getAtoms(); za++) {
                //Set<Integer> remaining_bonds = new HashSet<>();
                // NOTE!! do not consider by mistake the Uranium bonds.. ;)
                //for(int zb=0;zb<m_in.getBonds();zb++){ remaining_bonds.add(zb); }
                for(int bi_original : all_bonds) { // just loop over all_bonds, thats safest..
                    int bi = bond_map[bi_original];
                    //System.out.println("bond to cut now: "+bi_original + " in molecule coordinates: "+bond_map[bi_original]);
                    if(bi<0) {
                        continue;
                    }
                    if( taken_bonds.contains(bi_original) ) {
                        continue;
                    }

                    //NOTE!! These are always in original molecule coordinates !!
                    Set<Integer> splitset_a = new HashSet<>(partial_split.getLeft());
                    if(splitset_a.contains(bi_original)) {
                        continue;
                    }

//                    // check adjacency (we sort out splits with adjacent bonds in the cutset):
//                    boolean adjacent = false; int bi_atom_a = m_in.getBondAtom(0,bi); int bi_atom_b = m_in.getBondAtom(1,bi);
//                    for(int bond_i : partial_split.getLeft()) {
//                        int bond_i_a = m_in.getBondAtom(0,bond_i); int bond_i_b = m_in.getBondAtom(1,bond_i);
//                        adjacent =  bi_atom_a == bond_i_a || bi_atom_a == bond_i_b || bi_atom_b == bond_i_a || bi_atom_b == bond_i_b ;
//                        if(adjacent) {break;}
//                    }
//                    if(adjacent) {
//                        continue;
//                    }

                    splitset_a.add(bi_original);




                    StereoMolecule m_copy_a = new StereoMolecule(m);
                    //System.out.println("current m: smiles= "+ HyperspaceUtils.idcodeToSmiles(m_copy_a.getIDCode()));

                    if (!m.isRingBond(bi)) {
                        int[] splitset = splitset_a.stream().mapToInt(si -> si.intValue()).toArray();

                        //SynthonShredder.SplitResult sri = SynthonShredder.trySplit(m_in, splitset, 2);

                        // create split manually from m:
                        StereoMolecule m_c = new StereoMolecule(m);
                        m_c.ensureHelperArrays(Molecule.cHelperNeighbours);
                        int bond_cut_map[] = new int[m_c.getBonds()];
                        cutBondInStereoMolecule(m_c,bi,bond_cut_map);

                        // update bond map ( (original mol indeces) -> (partial_split_mol_indeces) :
                        int updated_bond_map[] = new int[m_in.getBonds()];
                        for(int zi=0;zi<m_in.getBonds();zi++) {
                            updated_bond_map[zi] = -1;
                            if( bond_map[ zi ] >= 0 ) {
                                updated_bond_map[zi] = bond_cut_map[bond_map[zi]];
                            }
                        }


                        // now it gets a bit complicated:
                        // compute fragments, then loop over atoms to collect the bonds of the two fragments..
                        // then copy a and b by bonds, AND rely on the "bond order is retained" feature
                        // to find out bond numbers in the copied fragments.

                        // NOTE!! On the uranium connector atoms:
                        //        They are considered in the bonds_in_a etc. array and sets. HOWEVER, we do not
                        //        consider them in the bond map arrays (obviously, as they have no bond in the original mol.)

                        int fragment_numbers[] = new int[m_c.getAtoms()];
                        m_c.getFragmentNumbers(fragment_numbers,false,true);

                        // NOTE!! These two sets bonds_in_a and bonds_in_b CONTAIN the connector Uranium atom bonds!!
                        // NOTE!! They are coordinates in m, NOT in original
                        Set<Integer> bonds_in_a = new HashSet<>();
                        Set<Integer> bonds_in_b = new HashSet<>();
                        boolean bonds_a[] = new boolean[m_c.getBonds()];
                        boolean bonds_b[] = new boolean[m_c.getBonds()];
                        for(int zi=0;zi<fragment_numbers.length;zi++) {
                            if(fragment_numbers[zi]==0) {
                                for(int bxi=0;bxi<m_c.getConnAtoms(zi);bxi++) {
                                    bonds_in_a.add( m_c.getConnBond(zi,bxi) );
                                    bonds_a[m_c.getConnBond(zi,bxi)] = true;
                                }
                            }
                            else if(fragment_numbers[zi]==1) {
                                for(int bxi=0;bxi<m_c.getConnAtoms(zi);bxi++) {
                                    bonds_in_b.add( m_c.getConnBond(zi,bxi) );
                                    bonds_b[m_c.getConnBond(zi,bxi)] = true;
                                }
                            }
                            else {
                                System.out.println("[ERROR] Wrong number of fragments..");
                            }
                        }
                        // compute bond map:
                        int bond_map_fragment_a_to_original[] = new int[bonds_in_a.size()];
                        int bond_map_fragment_b_to_original[] = new int[bonds_in_b.size()];
                        int bond_map_mol_to_original[] = new int[m_c.getBonds()];
                        // init to -1..
                        for(int zi=0;zi<bond_map_mol_to_original.length;zi++) { bond_map_mol_to_original[zi]=-1; }
                        for(int zi=0;zi<updated_bond_map.length;zi++) {

                            if(updated_bond_map[zi]<0) {continue;}
                            bond_map_mol_to_original[ updated_bond_map[zi] ] = zi;
                        }

                        int a_found_bonds = 0;
                        int b_found_bonds = 0;
                        for(int zi=0;zi<bond_map_mol_to_original.length;zi++) {
                            // 1. map to mol
                            if(bonds_in_a.contains(zi) ) {
                                bond_map_fragment_a_to_original[ a_found_bonds ] = bond_map_mol_to_original[zi];
                                a_found_bonds++;
                            }
                            else if(bonds_in_b.contains(zi)) {
                                bond_map_fragment_b_to_original[ b_found_bonds ] = bond_map_mol_to_original[zi];
                                b_found_bonds++;
                            }
                            else {
                                System.out.println("[ERROR] something very wrong..");
                            }
                        }

                        // now copy mols by bonds:
                        StereoMolecule fragment_a = new StereoMolecule();
                        StereoMolecule fragment_b = new StereoMolecule();

                        m_c.copyMoleculeByBonds(fragment_a,bonds_a,true,null);
                        m_c.copyMoleculeByBonds(fragment_b,bonds_b,true,null);

                        fragment_a.ensureHelperArrays(Molecule.cHelperCIP);
                        fragment_b.ensureHelperArrays(Molecule.cHelperCIP);

                        Set<Integer> bonds_in_a_original = new HashSet<>();
                        Set<Integer> bonds_in_b_original = new HashSet<>();

                        //for(int bxia : bonds_in_a) { bonds_in_a_original.add( bond_map_fragment_a_to_original[] ) }
                        //for(int bxia : bonds_in_b) { bonds_in_b_original[bxia] = true; }
                        for(int bxia : bond_map_fragment_a_to_original) {
                            if(bxia>=0) {
                                bonds_in_a_original.add(bxia);
                            }
                        }
                        for(int bxia : bond_map_fragment_b_to_original) {
                            if(bxia>=0) {
                                bonds_in_b_original.add(bxia);
                            }
                        }


                        // here we do the sanity check by also copying from the original molecule:
                        if(false) {
                            boolean bonds_from_original_a[] = new boolean[m_in.getBonds()];
                            boolean bonds_from_original_b[] = new boolean[m_in.getBonds()];
                            for(int bxia : bond_map_fragment_a_to_original) {
                                if(bxia>=0) {
                                    bonds_from_original_a[bxia] = true;
                                }
                            }
                            for(int bxia : bond_map_fragment_b_to_original) {
                                if(bxia>=0) {
                                    bonds_from_original_b[bxia] = true;
                                }
                            }

                            StereoMolecule sanity_check_molecule_a = new StereoMolecule();
                            StereoMolecule sanity_check_molecule_b = new StereoMolecule();

                            m_in.copyMoleculeByBonds(sanity_check_molecule_a,bonds_from_original_a,true,null);
                            m_in.copyMoleculeByBonds(sanity_check_molecule_b,bonds_from_original_b,true,null);

                            // print smiles:
                            System.out.println("sanity check:");
                            System.out.println( HyperspaceUtils.idcodeToSmiles( fragment_a.getIDCode() ) );
                            System.out.println( HyperspaceUtils.idcodeToSmiles( sanity_check_molecule_a.getIDCode() ) );
                            System.out.println( HyperspaceUtils.idcodeToSmiles( fragment_b.getIDCode() ));
                            System.out.println( HyperspaceUtils.idcodeToSmiles( sanity_check_molecule_b.getIDCode() ) );
                            System.out.println("mkay..");
                        }

                        // create split manually from m:
//                        StereoMolecule m_c = new StereoMolecule(m);
//                        m_c.ensureHelperArrays(Molecule.cHelperCIP);
//                        cutBondInStereoMolecule(m_c,bi);
//                        System.out.println("after cutting bond manually: smiles = "+HyperspaceUtils.idcodeToSmiles(m_c.getIDCode()));

//                        if(sri==null) {
//                            System.out.println("[ERROR] something very wrong.. sri==null");
//                            System.out.println("[ERROR] Data: cutsetsize="+zs+" splitset= "+splitset_a.toString());
//                            continue;
//                        }
//                        if(sri.fragments.length!=2) {
//                            System.out.println("[ERROR] something very wrong.. num_fragments!=2 -> "+sri.fragments.length);
//                            System.out.println("[ERROR] Data: cutsetsize="+zs+" splitset= "+splitset_a.toString());
//                            continue;
//                        }

                        all_split_sets.add(new BinarySplitResult(splitset_a,
                                new StereoMolecule[]{fragment_a,fragment_b},
                                new Set[]{ bonds_in_a_original , bonds_in_b_original }));
                        n_connectedness.get(zs).addAll(splitset_a);

                        //System.out.println("Ok Split, cutsetsize="+zs+" splitset= "+splitset_a.toString());
                        //System.out.println(sri.toString());
                    } else {
                        // perform the split (cut bi) and add the molecule:
                        StereoMolecule m_partial = new StereoMolecule(m);
                        m_partial.ensureHelperArrays(Molecule.cHelperNeighbours);
                        int cut_bond_map[] = new int[m_partial.getBonds()];
                        cutBondInStereoMolecule(m_partial,bi,cut_bond_map);
                        // update bond map ( (original mol indeces) -> (partial_split_moL_indeces) :
                        int updated_bond_map[] = new int[m_in.getBonds()];
                        for(int zi=0;zi<m_in.getBonds();zi++) {
                            updated_bond_map[zi] = -1;
                            if( bond_map[ zi ] >= 0 ) {
                                updated_bond_map[zi] = cut_bond_map[bond_map[zi]];
                            }
                        }
                        next_partial_splits.add(Triple.of(splitset_a,m_partial,updated_bond_map));
                        //System.out.println("partial Split, cutsetsize="+splitset_a.size()+" cutset= "+splitset_a.toString());
                        //StereoMolecule sm_copy = new StereoMolecule(m_partial);
                        //System.out.println( HyperspaceUtils.idcodeToSmiles(sm_copy.getIDCode()));
                    }
                }
            }
            current_partial_splits = next_partial_splits;
        }

        return all_split_sets;
    }

    /**
     *
     * @param mi
     * @param bi
     * @param cut_bond_map maps from original indeces to new indeces, -1 means deleted
     */
    public static void cutBondInStereoMolecule(StereoMolecule mi , int bi, int cut_bond_map[]) {
        mi.ensureHelperArrays(Molecule.cHelperNeighbours);

        for( int zi=0;zi<mi.getBonds();zi++) {
            if(zi==bi) { cut_bond_map[zi]=-1; }
            else {
                if(zi<bi) {
                    cut_bond_map[zi] = zi;
                }
                else {
                    cut_bond_map[zi] = zi-1;
                }
            }
        }

//        if(bi==1) {
//            System.out.println( mi.getBondAtom(0,bi) +" / "+mi.getBondAtom(1,bi) );
//            System.out.println( mi.getBondType(bi) );
//        }

        int bond_type = mi.getBondType(bi);
        int bia1 = mi.getBondAtom(0,bi);
        int bia2 = mi.getBondAtom(1,bi);
        int idx_bia1 = mi.addAtom(92);
        int idx_bia2 = mi.addAtom(92);
        //System.out.println("new bond indeces: "+idx_bia1+" , "+idx_bia2);
        mi.addBond(bia1,idx_bia1,bond_type);
        mi.addBond(bia2,idx_bia2,bond_type);
        mi.deleteBond(bi);
        mi.ensureHelperArrays(Molecule.cHelperCIP);

//        if(bi==1) {
//            System.out.println("check, after:");
//            System.out.println( mi.getBondAtom(0,bi) +" / "+mi.getBondAtom(1,bi) );
//            System.out.println( mi.getBondType(bi) );
//            System.out.println("mkay");
//        }

    }


    public static class NestedBinarySplits {

        int bitset_length = 0;
        public Set<Integer> cut_set = new HashSet<>();
        //public Map<String,Set<Integer>> bond_sets = new HashMap<>();
        public Map<String,BitSet> bond_sets = new HashMap<>();

        public NestedBinarySplits(NestedBinarySplits si) {
            this.cut_set = new HashSet<>(si.cut_set);
            this.bond_sets = new HashMap<>(si.bond_sets);
        }

        public NestedBinarySplits(BinarySplitResult initial, int bitset_length) {
//            this.bond_sets.put(initial.cutset.toString()+"_a",new HashSet<>(initial.bonds_in_splits[0]));
//            this.bond_sets.put(initial.cutset.toString()+"_b",new HashSet<>(initial.bonds_in_splits[1]));
            this.bitset_length = bitset_length;
            BitSet bs_a = new BitSet(bitset_length);
            BitSet bs_b = new BitSet(bitset_length);
            for(Integer bxi : initial.bonds_in_splits[0]) {bs_a.set(bxi,true);}
            for(Integer bxi : initial.bonds_in_splits[1]) {bs_b.set(bxi,true);}
            this.bond_sets.put(initial.cutset.toString()+"_a",bs_a);
            this.bond_sets.put(initial.cutset.toString()+"_b",bs_b);

            this.cut_set.addAll(initial.cutset);

        }

        /**
         * Returns the bond_set id of the bond_set that can be split with this next_split.
         * Or null if no split is possible.
         *
         * @param next_split
         * @return
         */
        public String checkIfSplitIsPossible(BinarySplitResult next_split) {
            // loop over bond sets to see if we can use one:
            for(String bid : bond_sets.keySet()) {
//                if( bond_sets.get(bid).containsAll( next_split.cutset ) ) {
//                    return bid;
//                }
                boolean contains_all = true;
                BitSet bsi = bond_sets.get(bid);
                for(Integer ci : next_split.cutset) {
                    contains_all &= bsi.get(ci);
                }
                if(contains_all) { return bid; }
            }
            return null;
        }

        /**
         * Performs the split.
         *
         * @param bond_set_id
         * @param next_split
         */
        public void applyBondSetSplit(String bond_set_id, BinarySplitResult next_split) {
            // 1. split all bonds from bond_set according to next_split:
            //Set<Integer> new_bonds_a = new HashSet<>();
            //Set<Integer> new_bonds_b = new HashSet<>();

            //Set<Integer> to_split = this.bond_sets.remove(bond_set_id);
            BitSet to_split = this.bond_sets.remove(bond_set_id);

            BitSet bond_sets[] = next_split.getBondSetsAsBitSets(this.bitset_length);

            BitSet new_bonds_a = (BitSet) to_split.clone();
            BitSet new_bonds_b = (BitSet) to_split.clone();

            new_bonds_a.and(bond_sets[0]);
            new_bonds_b.and(bond_sets[0]);


//            for(int si : to_split) {
//                if(next_split.bonds_in_splits[0].contains(si)) {
//                    new_bonds_a.add(si);
//                }
//                else if(next_split.bonds_in_splits[1].contains(si)) {
//                    new_bonds_b.add(si);
//                }
//                else {
//                    // now this must be the connector: (we do this as sanity check..)
//                    if(!next_split.cutset.contains(si)) {
//                        System.out.println("[ERROR] something very very wrong");
//                    }
//                }
//            }
            this.cut_set.addAll(next_split.cutset);
            this.bond_sets.put( bond_set_id+"_"+next_split.cutset.toString()+"_a" , new_bonds_a );
            this.bond_sets.put( bond_set_id+"_"+next_split.cutset.toString()+"_b" , new_bonds_b );
        }

    }


    /**
     * NOTE! This does compute a little less splits compared to brute-force. The reason for this is,
     * that it is not possible to split edges directly next to connector edges again.
     *
     * @param mol
     * @param max_cuts
     * @param max_fragments
     * @return
     */
    public static Map<Set<Integer>, SynthonShredder.SplitResult> computeAllValidSplits(StereoMolecule mol, int max_cuts , int max_fragments) {

        long ts_2_a = System.currentTimeMillis();

        List<Pair<SynthonShredder2.BinarySplitResult,SynthonShredder2.BinarySplitResult>> all_splits = new ArrayList<>();
        List<SynthonShredder2.BinarySplitResult> splits_binary = SynthonShredder2.computeAllBinarySplits(mol,max_cuts);
        //expand to full, i.e. take all combinations that sum up to max num_splits splits and that do not
        //share any edges
        Map<Set<Integer>,Integer> all_reported_cutsets = new HashMap<>(); // maps cutsets to number of fragments

        Map<Set<Integer>, SynthonShredder.SplitResult> split_results_via_binary_splits = new HashMap<>();

        // now combine splits, such that we find all possible combinations..
        // this contains all that we newly created in the last round of expansion, i.e. these we have to further expand
        Set<NestedBinarySplits> current_cutsets = new HashSet<>();

        for(int zf = 1 ; zf <= max_fragments - 1 ; zf++) {
            if(zf==1) {
                for(BinarySplitResult bsr : splits_binary) {
                    all_reported_cutsets.put(bsr.cutset,2);
                    current_cutsets.add( new NestedBinarySplits(bsr,mol.getBonds()) );
                }
            }
            else {
                Set<NestedBinarySplits> next_cutsets = new HashSet<>();
                // try to expand all current_cutset
                for(NestedBinarySplits si : current_cutsets) {
                    for( BinarySplitResult next_bsri : splits_binary ) {

                        if(si.cut_set.size()+next_bsri.cutset.size() > max_cuts) {
                            continue;
                        }

                        Set<Integer> combined_connectors = new HashSet<>();
                        combined_connectors.addAll(si.cut_set);combined_connectors.addAll(next_bsri.cutset);
                        if(combined_connectors.size()<si.cut_set.size()+next_bsri.cutset.size()) {
                            continue;
                        }

                        String partition_to_split = si.checkIfSplitIsPossible(next_bsri);
                        if(partition_to_split!=null) {
                            NestedBinarySplits si_c = new NestedBinarySplits(si);
                            si_c.applyBondSetSplit(partition_to_split,next_bsri);
                            all_reported_cutsets.put(si_c.cut_set,si_c.bond_sets.size());
                            next_cutsets.add(si_c);
                        }
                    }
                }
                current_cutsets = next_cutsets;
            }
        }

        // now create all the SplitResults:
        for(Set<Integer> ci : all_reported_cutsets.keySet()) {
            SynthonShredder.SplitResult sri = SynthonShredder.trySplit(mol,ci.stream().mapToInt(ii->ii).toArray(),max_fragments);
            //SynthonShredder.SplitResult sri = null;
            if(sri!=null) {
                split_results_via_binary_splits.put(ci, sri);
            }
            else {
                //System.out.println("mkay..");
            }
        }

        long ts_2_b = System.currentTimeMillis();

        return split_results_via_binary_splits;
    }

}
