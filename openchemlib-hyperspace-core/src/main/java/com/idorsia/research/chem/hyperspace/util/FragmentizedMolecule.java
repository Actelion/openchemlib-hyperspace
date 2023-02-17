package com.idorsia.research.chem.hyperspace.util;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import org.apache.commons.lang3.tuple.Pair;

import java.util.*;
import java.util.stream.Collectors;

public class FragmentizedMolecule {

    StereoMolecule M;
    int frag_size;

    List<BitSet> F;
    List<StereoMolecule> F_asM;

    Map<String,List<Integer>> F_asIdCodeSet;

    Map<BitSet,StereoMolecule> F2;
    Map<String,Set<BitSet>> F_by_idcode;

    BitSet[] fragment_adjacency_matrix; // entry i contains all adjacency information (adjacent meaning fragment overlpa distance 1)

    // size: num_atoms * num_fragments
    // atom i, fragment j is at position ( i +  )
    List<BitSet> Map_atom_to_frag;

    public FragmentizedMolecule(StereoMolecule m, int frag_size) {
        this.M = m;
        this.frag_size = frag_size;
    }

    public void init() {
        FragmentGenerator fg = new FragmentGenerator(this.M,frag_size,frag_size);
        this.M = fg.getCanonizedInputMolecule();
        fg.computeFragments();
        this.F = fg.getFragments();
        Map_atom_to_frag = new ArrayList<>(this.M.getAtoms());
        for(int zi=0;zi<this.M.getAtoms();zi++) {
            this.Map_atom_to_frag.add( new BitSet( F.size() ));
        }

        this.F2 = new HashMap<>();
        this.F_asM = new ArrayList<>(this.F.size());
        this.F_by_idcode = new HashMap<>();
        this.F_asIdCodeSet = new HashMap<>();

        for( int zi=0;zi<this.F.size();zi++) {
            BitSet fi = this.F.get(zi);
            StereoMolecule mi = new StereoMolecule();
            boolean frag_as_bool_array[] = new boolean[this.M.getAtoms()];
            for(int za=0;za<frag_as_bool_array.length;za++) {
                frag_as_bool_array[za] = fi.get(za);
            }
            this.M.copyMoleculeByAtoms( mi , frag_as_bool_array , true , null );
            mi.ensureHelperArrays(Molecule.cHelperCIP);
            String fmi_idc = mi.getIDCode();
            this.F_asM.add(mi);
            this.F2.put(fi,mi);
            for(int za=0;za<frag_as_bool_array.length;za++) {
                this.Map_atom_to_frag.get(za).set( zi , fi.get(za) );
            }
            if(!this.F_by_idcode.containsKey(fmi_idc)){this.F_by_idcode.put(fmi_idc,new HashSet<>());}
            this.F_by_idcode.get(fmi_idc).add( fi );

            if(!this.F_asIdCodeSet.containsKey( fmi_idc )){this.F_asIdCodeSet.put(fmi_idc,new ArrayList<>());}
            this.F_asIdCodeSet.get(fmi_idc).add(zi);
        }


        // init adjacency information:
        this.fragment_adjacency_matrix = new BitSet[this.F.size()];
        for( int zi=0;zi<this.F.size();zi++) {
            BitSet adjacencies_i = new BitSet(this.F.size());
            BitSet f_i = this.F.get(zi);
            for( int zj=0;zj<this.F.size();zj++) {
                BitSet f_j = (BitSet) this.F.get(zj).clone();
                int card_j = f_j.cardinality();
                f_j.and(f_i);
                if(f_j.cardinality()==card_j-1) {
                    adjacencies_i.set(zj,true);
                }
            }
            this.fragment_adjacency_matrix[zi] = adjacencies_i;
        }
    }

    public StereoMolecule getM() {
        return this.M;
    }

    public BitSet getFragment(int i) {
        return this.F.get(i);
    }

    public boolean containsFragment(String fidc) {
        return this.F_by_idcode.containsKey(fidc);
    }

    public boolean testAdjacency(String fa, String fb) {
        List<Integer> fia = this.F_asIdCodeSet.get(fa);
        List<Integer> fib = this.F_asIdCodeSet.get(fb);
        boolean adjacent = false;
        for(int za=0;za<fia.size();za++) {
            for(int zb=0;zb<fib.size();zb++) {
                adjacent |= this.fragment_adjacency_matrix[fia.get(za)].get(fib.get(zb));
            }
        }
        return adjacent;
    }

    /**
     * This will return a value for any mask, not only for the fragment atom masks.
     *
     */
    public StereoMolecule getFragmentByAtomMask(BitSet atom_mask) {
        if(this.F2.containsKey(atom_mask)) {
            return this.F2.get(atom_mask);
        }
        else {
            StereoMolecule mi = new StereoMolecule();
            boolean frag_as_bool_array[] = new boolean[this.M.getAtoms()];
            for(int za=0;za< Math.min(atom_mask.size(),frag_as_bool_array.length) ;za++) {
                frag_as_bool_array[za] = atom_mask.get(za);
            }
            this.M.copyMoleculeByAtoms( mi , frag_as_bool_array , true , null );
            mi.ensureHelperArrays(Molecule.cHelperCIP);
            return mi;
        }
    }

    public List<StereoMolecule> getFragmentsAsStereoMolecules() {
        return this.F_asM;
    }


    /**
     * We start by assembling the uniquely defined ones into connected components.
     *
     * @param frag_idcodes
     * @param max_enumerated
     * @return
     */
    public List<StereoMolecule> assembleConnectedSubstructuresFromFragments(List<String> frag_idcodes, boolean accept_missing_idcodese, int max_enumerated) {
        List<List<BitSet>> bitsets = new ArrayList<>();
        for(String fidc : frag_idcodes) {
            Set<BitSet> bitsets_i = this.F_by_idcode.get(fidc);
            if(bitsets==null) {
                if(!accept_missing_idcodese) {
                    return null;
                }
                else {
                    continue;
                }
            }
            else {
                bitsets.add(new ArrayList<>(bitsets_i));
            }
        }


        // now start assembling unique ones.
        List<BitSet> bitsets_unique = bitsets.stream().filter( fi -> fi.size()==1 ).map( fi2 -> fi2.get(0) ).collect(Collectors.toList());

        List<BitSet> connected_components = new ArrayList<>();

        List<BitSet> available = new ArrayList<>(bitsets_unique);

        while( ! available.isEmpty() ) {
            BitSet current_component = available.iterator().next();
            available.remove(current_component);
            boolean expanding = true;
            while(expanding) {
                expanding = false;
                // find overlapping:
                for (BitSet bi_ext : available) {
                    BitSet bi_ext_i = (BitSet) bi_ext.clone();
                    int cardinality_bi_ext = bi_ext.cardinality();
                    bi_ext_i.and(current_component);
                    if (bi_ext_i.cardinality() == cardinality_bi_ext - 1) {
                        expanding = true;
                        current_component.or( bi_ext );
                        available.remove(bi_ext);
                        break;
                    }
                }
            }
            connected_components.add(current_component);
        }


        // then we consider just the largest connected component:
        connected_components.sort( (x,y) -> - Integer.compare( x.size() , y.size() ) );
        BitSet bitsets_unique_combined = connected_components.get(0);


        //BitSet bitsets_unique_combined = bitsets_unique.stream().reduce( (x,y) -> bitset_or(x,y) ).get();

//        System.out.println("frags[idcode]");
//        for(String fci : frag_idcodes) {
//            System.out.println(fci);
//        }
//        System.out.println();

        // then extend this in all possible ways
        List<List<BitSet>> bitsets_nonunique = bitsets.stream().filter( fi -> fi.size()>1 ).collect(Collectors.toList());

        // remove all non-unique ones that have an option that is already covered:
        List<List<BitSet>> bitsets_nonunique_not_covered = new ArrayList<>();

        for(List<BitSet> non_unique_i : bitsets_nonunique) {
            boolean covered = false;
            for(BitSet bsi : non_unique_i) {
                BitSet bsi2 = (BitSet) bsi.clone();
                bsi2.and(bitsets_unique_combined);
                if(bsi2.equals(bsi)) {
                    covered = true;
                    break;
                }
                else {
                }
            }
            if(!covered) {
                bitsets_nonunique_not_covered.add(non_unique_i);
            }
        }


        StereoMolecule m_connected = this.getFragmentByAtomMask( bitsets_unique_combined );
        System.out.println("connected: "+ HyperspaceUtils.idcodeToSmiles(m_connected.getIDCode()));

        List<Pair<BitSet,List<List<BitSet>>>> enumerated_and_available_bitsets = new ArrayList<>();

        enumerated_and_available_bitsets.add( Pair.of( bitsets_unique_combined , bitsets_nonunique_not_covered ) );

        boolean extending = true;
        while( extending ) {
            extending = false;
            List<Pair<BitSet,List<List<BitSet>>>> next_enumerated_and_available_bitsets = new ArrayList<>();
            for( Pair<BitSet,List<List<BitSet>>> pi : enumerated_and_available_bitsets) {
                BitSet atom_mask = pi.getLeft();
                for( List<BitSet> extensions : pi.getRight()) {
                    List<List<BitSet>> extensions_without_i = new ArrayList<>( pi.getRight() );
                    extensions_without_i.remove(extensions);
                    for( BitSet bi_ext : extensions ) {
                        BitSet bi_ext_i = (BitSet) bi_ext.clone();
                        int cardinality_bi_ext = bi_ext.cardinality();
                        bi_ext_i.and(atom_mask);
                        if( bi_ext_i.cardinality() == cardinality_bi_ext - 1 ) {
                            // this means we extend:
                            extending = true;
                            BitSet mask_extended = (BitSet) bi_ext.clone();
                            mask_extended.or(atom_mask);
                            next_enumerated_and_available_bitsets.add( Pair.of( mask_extended , extensions_without_i ) );
                        }
                    }
                }
            }
            if(extending) {
                enumerated_and_available_bitsets = next_enumerated_and_available_bitsets;
            }
        }

        Map<BitSet,StereoMolecule> mi = new HashMap<>();
        for( Pair<BitSet,List<List<BitSet>>> variant_i : enumerated_and_available_bitsets ) {
            if(!mi.containsKey(variant_i.getLeft())) {
                StereoMolecule mxi = this.getFragmentByAtomMask(variant_i.getLeft());
                mi.put(variant_i.getLeft(),mxi);
            }
        }

        return new ArrayList<>( mi.values() );
    }

    /**
     * NOTE: this is not necessary uniquely defined, it enumerates all versions.
     *
     * NOTE: returns null if it encounters unknown idcodes and accept_missing_idcodes is false.
     *
     * @param frag_idcodes
     * @return
     */
    public List<StereoMolecule> getSubmoleculesFromFragments(List<String> frag_idcodes, boolean accept_missing_idcodese, int max_enumerated) {
        List<List<BitSet>> bitsets = new ArrayList<>();
        for(String fidc : frag_idcodes) {
            Set<BitSet> bitsets_i = this.F_by_idcode.get(fidc);
            if(bitsets==null) {
                if(!accept_missing_idcodese) {
                    return null;
                }
                else {
                    continue;
                }
            }
            else {
                bitsets.add(new ArrayList<>(bitsets_i));
            }
        }

        List<BitSet> all_combined_bitsets = new ArrayList<>();
        if(bitsets.stream().mapToInt( fi -> fi.size() ).allMatch( si -> si==1 )) {
            BitSet or_bitsets= all_combined_bitsets.stream().reduce(  (x,y) -> bitset_or(x,y) ).get();
            all_combined_bitsets.add(or_bitsets);
        }
        else {
            List<BitSet> all_enumerated = bitsets.stream().reduce( (x,y) -> compute_all_or_versions(x,y,max_enumerated) ).get();
            all_combined_bitsets.addAll(all_enumerated);
        }


        // now build the molecules:
        Set<BitSet> all_bitsets = new HashSet<>(all_combined_bitsets);
        List<StereoMolecule> all_molecules = new ArrayList<>();
        for(BitSet bsi : all_bitsets) {
            StereoMolecule mi = this.getFragmentByAtomMask(bsi);
            all_molecules.add(mi);
        }

        return all_molecules;
    }

    public static List<BitSet> compute_all_or_versions( List<BitSet> x, List<BitSet> y , int max_size) {

        if(x.size()==1 && y.size() == 1) {
            List<BitSet> all_or_versions = new ArrayList<>();
            all_or_versions.add(bitset_or(x.get(0),y.get(0)));
            return all_or_versions;
        }

        int cnt = 0;
        Set<BitSet> all_or_versions = new HashSet<>();
        for(int zx=0;zx<x.size();zx++) {
            for(int zy=0;zy<y.size();zy++) {
                BitSet bsi = bitset_or( x.get(zx) , y.get(zy) );
                all_or_versions.add(bsi);
                if(cnt++>max_size) {
                    System.out.println("[WARNING] exceeded size in combining bitsets..");
                    break;
                }
            }
        }

        return new ArrayList<>(all_or_versions);
    }

    public static BitSet bitset_or(BitSet x, BitSet y) {
        BitSet xc = (BitSet) x.clone();
        xc.or(y);
        return xc;
    }

    /**
     * Distance is specified in overlap of bitsets. I.e. distance one means,
     * that all except for one atom overlap.
     *
     *
     * @param min_distance
     * @param max_distance
     * @return
     */
    public List<Pair<Integer,BitSet>> getNeighborFragments(BitSet fragment, int min_distance, int max_distance) {
        List<Pair<Integer,BitSet>> fragments = new ArrayList<>();
        for(int zi=0;zi<this.F.size();zi++) {
            BitSet fi = (BitSet) this.F.get(zi).clone();
            int cardinality_a = fi.cardinality();
            fi.and( fragment );
            int overlap = fi.cardinality();
            int distance = cardinality_a - overlap;
            if(distance >= min_distance && distance <= max_distance) {
                fragments.add(Pair.of(overlap,(BitSet) this.F.get(zi).clone()));
            }
        }
        return fragments;
    }

    /**
     * Distance is specified in overlap of bitsets. I.e. distance one means,
     * that all except for one atom overlap.
     *
     *
     * @param min_distance
     * @param max_distance
     * @return
     */
    public List<Pair<Integer,Integer>> getNeighborFragmentIndices(BitSet fragment, BitSet forbidden, int min_distance, int max_distance) {
        List<Pair<Integer,Integer>> fragments = new ArrayList<>();
        for(int zi=0;zi<this.F.size();zi++) {
            if(forbidden.get(zi)) {continue;}
            BitSet fi = (BitSet) this.F.get(zi).clone();
            int cardinality_a = fi.cardinality();
            fi.and( fragment );
            int overlap = fi.cardinality();
            int distance = cardinality_a - overlap;
            if(distance >= min_distance && distance <= max_distance) {
                fragments.add(Pair.of(overlap,zi));
            }
        }
        return fragments;
    }


    /**
     * Loops over all fragments that are neighboring fragments to scaffold fragments, and
     * checks if they cover the expansion atoms
     *
     * @param scaffold_fragments
     * @param atom_i
     */
    public BitSet computeFragmentsRequiredByExpansion(BitSet scaffold_fragments, int atom_i) {

        // we fetch the atom containing fragments, and the we sort out the ones
        // that are not adjacent to any

        BitSet atom_containing = this.Map_atom_to_frag.get(atom_i);

        BitSet required_fragments = new BitSet(this.F.size());

        for( int fi : scaffold_fragments.stream().toArray() ) {
            BitSet adjacent_i = (BitSet) this.fragment_adjacency_matrix[fi].clone();
            adjacent_i.or(atom_containing);
            required_fragments.or(  adjacent_i );
        }

        return required_fragments;
    }
}
