package com.idorsia.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */


import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.*;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerStandard2DFactory;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class SynthonSpace implements Serializable {

    private static final long SerialVersionUID=228844663377L;
    //private int logLevel = 0;

    private static int logLevel_hitExpansion   = 1;
    private static int logLevel_findCandidates = 1;

    //private final int BITS = 512;
    protected int BITS = 1024;

    public static final int MAX_CONNECTORS = 8;

    // max considered distance for the connector proximal region
    public static final int CONNECTOR_REGION_SIZE = 3;

    //private final int BITTREE_BIN_SIZE = 1024;
    private int BITTREE_BIN_SIZE = 64;

    public int getBits() {return this.BITS;}

    public String getSpaceInfoString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Reactions:        "+fragment_map.keySet().size() + "   Fragment Types: "+fragment_map.values().stream().flatMap( vi -> vi.values().stream() ).count() + "\n");
        sb.append("Total fragments:  "+fragment_map.values().stream().flatMap( vi -> vi.values().stream() ).mapToInt( li -> li.size() ).sum() + "\n" );

        long num_structures = 0;
        for(String ri : fragment_map.keySet()) {
            num_structures += fragment_map.get(ri).values().stream().mapToLong( li -> li.size() ).reduce( (x,y) -> x*y ).getAsLong();
        }
        sb.append("Structures:       "+ String.format( "%.2f Mio."  , (num_structures/1000000.0) )  + "\n\n");
        return sb.toString();
    }


    public static void setLogLevel(int p_logLevel_hitExpansion, int p_logLevel_findCandidates) {
        logLevel_hitExpansion    = p_logLevel_hitExpansion;
        logLevel_findCandidates  = p_logLevel_findCandidates;
    }

    //Map<String,Map<Integer,Map<BitSet,String>>> descriptors_sorted = new HashMap<>();

    public static final class FragType implements Serializable , Comparable<FragType> {
        public final String rxn_id;
        public final int frag;
        public final int hash;
        public FragType(String rxn_id,int frag) {
            this.rxn_id = rxn_id;
            this.frag   = frag;
            this.hash   = this.computHashCode();
        }
        public String toString() {
            return rxn_id + ":" + this.frag;
        }
        public boolean equals(Object o) {
            if(! (o instanceof FragType) ) {
                return false;
            }
            return o.toString().equals(this.toString());
        }
        public int hashCode() {
            return this.hash;
        }
        private int computHashCode() {
            return this.toString().hashCode();
        }

        @Override
        public int compareTo(FragType fragType) {
            return this.toString().compareTo(fragType.toString());
        }
    }

//    public Set<Integer> getConnectorSetForFragType(FragType ft) {
//        return this.fragment_map.get(ft.rxn_id).get(ft.frag).get(0).getConnectors();
//    }

    public BitSet getConnectorSetForFragType(FragType ft) {
        return this.fragment_map.get(ft.rxn_id).get(ft.frag).get(0).getConnectors();
    }



    public static final class FragId implements Serializable {
        private static final long SerialVersionUID=1116644677324L;

        public final String rxn_id;
        public final int frag;
        public final String idcode;
        public final String fragment_id;

        public final BitSet fp;
        public final BitSet fp_non_unique_connectors;

        // TODO: remove these..
        //public final int uranium_pos_blue;
        //public final int uranium_pos_red;
        //public final int uranium_pos_orange;

//        public final boolean has_blue;
//        public final boolean has_red;
//        public final boolean has_orange;
//        public final boolean has_green;

//        public final Set<Integer> connectors;
        BitSet connectors;

        public final int hash;

        public FragId(String rxn_id, int frag, StereoMolecule sm, String fragment_id, BitSet fp, BitSet fp_non_unique_connectors) {
            this.rxn_id = rxn_id;
            this.frag = frag;
            this.idcode = ""+new String( sm.getIDCode() );
            this.fragment_id = fragment_id;
            this.fp = fp;
            this.fp_non_unique_connectors = fp_non_unique_connectors;

            // determine the connector indicators:
//            boolean b_blue   = false;
//            boolean b_red    =  false;
//            boolean b_orange = false;
//            boolean b_green = false;

            BitSet connis = determineConnectorsInMolecule(sm);
            this.connectors = connis;

//            this.has_blue   = b_blue;
//            this.has_red    = b_red;
//            this.has_orange = b_orange;
//            this.has_green  = b_green;

//            HashSet<Integer> connector_set = new HashSet<>();
//            if(this.has_blue){   connector_set.add(92); }
//            if(this.has_red){    connector_set.add(93); }
//            if(this.has_orange){ connector_set.add(94); }
//            if(this.has_green){  connector_set.add(95); }
//            this.connectors = connector_set;

            this.hash       = this.toString().hashCode();
            //uranium_pos_blue = -1;
            //uranium_pos_red = -1;
            //uranium_pos_orange = -1;
            // identify and store the colors of the uranium placeholdes:
//            int pos_placeholder_blue    = -1;
//            int pos_placeholder_red     = -1;
//            int pos_placeholder_orange  = -1;
//            for(int zi=0;zi<sm.getAtoms();zi++) {
//                if(sm.getAtomicNo(zi)==92) {
//                    if(sm.getAtomColor(zi)==StereoMolecule.cAtomColorBlue) {
//                        pos_placeholder_blue = zi;
//                    }
//                    else if(sm.getAtomColor(zi)==StereoMolecule.cAtomColorRed) {
//                        pos_placeholder_red = zi;
//                    }
//                    else if(sm.getAtomColor(zi)==StereoMolecule.cAtomColorOrange) {
//                        pos_placeholder_orange = zi;
//                    }
//                }
//            }
//            uranium_pos_blue    = pos_placeholder_blue;
//            uranium_pos_red     = pos_placeholder_red;
//            uranium_pos_orange  = pos_placeholder_orange;
        }
        public FragId(String rxn_id, int frag, String idcode, String fragment_id, BitSet connector_bitset, BitSet fp,  BitSet fp_non_unique_connectors) {
            this.rxn_id = rxn_id;
            this.frag = frag;
            this.idcode = idcode;
            this.fragment_id = fragment_id;


            this.fp = fp;
            this.fp_non_unique_connectors = fp_non_unique_connectors;

            HashSet<Integer> connector_set = new HashSet<>();
//            if(this.has_blue){   connector_set.add(92); }
//            if(this.has_red){    connector_set.add(93); }
//            if(this.has_orange){ connector_set.add(94); }
//            if(this.has_green){  connector_set.add(95); }

            this.connectors = connector_bitset;

            this.hash       = this.toString().hashCode();

//            this.uranium_pos_blue    = pos_blue;
//            this.uranium_pos_red     = pos_red;
//            this.uranium_pos_orange  = pos_orange;
        }

        private void initConnectorIndicators() {

        }


//        public Set<Integer> getConnectors() {
//            return connectors;
//        }

        public BitSet getConnectors() {
            return connectors;
        }

        public String toString() {
            //return rxn_id + ":" + this.frag+" -> "+this.idcode + " ["+fragment_id+"]"+" ["+uranium_pos_blue+","+uranium_pos_red+","+uranium_pos_orange+"]";
            //return rxn_id + ":" + this.frag+" -> "+this.idcode + " ["+fragment_id+"]"+" [" + (has_blue?"1":"0")+","+(has_red?"1":"0")+","+(has_orange?"1":"0")+","+(has_green?"1":"0")+ "]";
            return rxn_id + ":" + this.frag+" -> "+this.idcode + " ["+fragment_id+"]"+" [" + createBitSetString(this.connectors,MAX_CONNECTORS) + "]";

        }
        public boolean equals(Object o) {
            if(! (o instanceof FragId) ) {
                return false;
            }
            return o.toString().equals(this.toString());
        }

        public int hashCode() {
            return this.hash;
        }

    }


    public static class RxnId implements Serializable {
        private static final long SerialVersionUID=7773366551928L;

        public final String id;
        public final int num_connectors;
        public final List<Integer> connector_counts;
        public final List<BitSet>  connectors;
        public RxnId(String id, int num_connectors, List<Integer> connector_counts, List<BitSet> connectors) {
            this.id = id;
            this.num_connectors = num_connectors;
            this.connector_counts = connector_counts;
            this.connectors = connectors;
        }
    }


    public static String createBitSetString(BitSet bs, int length) {
        StringBuilder sb = new StringBuilder();
        for(int zi=0;zi<length;zi++){ sb.append( (bs.get(zi)?"1":"0")); }
        return sb.toString();
    }
    public static BitSet parseBitSetString(String si) {
        char[] seq = si.toCharArray();
        BitSet bsi = new BitSet(si.length());
        for(int zi=0;zi<si.length();zi++) {
            if(seq[zi]=='1') {
                bsi.set(zi);
            }
            else if(seq[zi]=='0') {

            }
            else {
                System.out.println("ERROR: problem parsing BitSetString, encountered character: "+seq[zi]);
            }
        }
        return bsi;
    }

    /**
     * Creates connector bitset, i.e. it finds all atoms with
     * atomic weight inbetween 92 and 99.
     *
     * @return
     */
    public static BitSet determineConnectorsInMolecule(StereoMolecule mi) {
        mi.ensureHelperArrays(Molecule.cHelperNeighbours);
        BitSet bi = new BitSet(MAX_CONNECTORS);
        for(int zi=0;zi<mi.getAtoms();zi++) {
            int an = mi.getAtomicNo(zi);
            if(an>=92 && an<(92+MAX_CONNECTORS) ) {
                bi.set(an-92);
            }
        }
        return bi;
    }

    public static final class InitialHit {
        public final BitSet first_fp;
        public final FragType fragment_type;
        public final FragId primary_fragment_id;
        //public final List<BitSet> remaining_fps;

        public final Map<FragType,BitSet> remaining_fps;
        public final Map<FragType,StereoMolecule> remaining_frags;

        //public final Set<Integer> cut_set;
        public final String cut_set_hash;
        // index of cutset fragment idcode maps to the frag-defining BitSet
        // this is defined "reverse" because the BitSets may not be unique!
        public final Map<Integer,BitSet> target_fragments_in_cutset;

        //public InitialHit(BitSet first_fp, String rxn, int frag, List<BitSet> remaining, Set<Integer> cut_set, Map<Integer,BitSet> target_fragments_in_cutset) {
        //public InitialHit(BitSet first_fp, String rxn, int frag, FragId primary_fragment_id, List<BitSet> remaining, String cut_set_hash, Map<Integer,BitSet> target_fragments_in_cutset) {
        public InitialHit(BitSet first_fp, String rxn, int frag, FragId primary_fragment_id, Map<FragType,BitSet> remaining, Map<FragType,StereoMolecule> remaining_frags, String cut_set_hash, Map<Integer,BitSet> target_fragments_in_cutset) {
            this.first_fp      = first_fp;
            this.fragment_type = new FragType(rxn,frag);
            this.primary_fragment_id = primary_fragment_id;
            this.remaining_fps   = remaining;
            this.remaining_frags = remaining_frags;
            //this.cut_set       = cut_set;
            this.cut_set_hash  = cut_set_hash;
            this.target_fragments_in_cutset = target_fragments_in_cutset;
        }
    }



    /**
     * NOTE!!
     * Only the BitSet for the PRIMARY FRAGMENT is a PathFP of
     * a building block!
     * The other BitSets indicate only the required bits, i.e. they have to be expanded with
     * help of the BitSet tree of the reaction fragment!!!
     *
     */
    public static final class Hit {

        public final String rxn;
        public final Map<Integer,BitSet> fragments;
        public final FragId primary_fragment_id;
        public final int primary_fragment;

        //public final Set<Integer> cut_set;
        public final String cut_set_hash;

        // index of cutset fragment idcode maps to the frag-defining BitSet
        // this is defined "reverse" because the BitSets may not be unique!
        public final Map<Integer,BitSet> target_fragments_in_cutset;

        //public Hit(String rxn, Map<Integer,BitSet> fragments, int primary_fragment, Set<Integer> cutset, Map<Integer,BitSet> target_fragments_in_cutset) {
        public Hit(String rxn, Map<Integer,BitSet> fragments, FragId primary_fragment_id, int primary_fragment, String cutset_hash, Map<Integer,BitSet> target_fragments_in_cutset) {
            this.rxn = rxn;
            this.fragments = fragments;
            this.primary_fragment_id = primary_fragment_id;
            this.primary_fragment = primary_fragment;
            //this.cut_set = cutset;
            this.cut_set_hash = cutset_hash;
            this.target_fragments_in_cutset = target_fragments_in_cutset;
        }

    }

    /**
     * This represents a combinatorial hit, where all the stored FragIds themselves
     * passed the synthon substructure check.
     *
     */
    public static final class CombinatorialHit implements Serializable {
        public final String rxn;
        public final Map<FragType,List<FragId>> hit_fragments;
        public final SynthonShredder.SplitResult sri;
        public final Map<Integer,Pair<FragType,BitSet>> mapping;
        public CombinatorialHit(String rxn, Map<FragType,List<FragId>> hit_fragments, SynthonShredder.SplitResult split_result, Map<Integer,Pair<FragType,BitSet>> mapping) {
            this.rxn = rxn;
            this.hit_fragments = hit_fragments;
            this.sri = split_result;
            this.mapping = mapping;
        }
    }

    /**
     * This represents an expanded hit with the full information about all fragments.
     */
    public static final class ExpandedHit {
        public final String rxn;
        public final Map<Integer,FragId> hit_fragments;
        public ExpandedHit(String rxn, Map<Integer,FragId> fragments) {
            this.rxn = rxn;
            this.hit_fragments = fragments;
        }

        public int hashCode() {
            return  hit_fragments.hashCode() ^ ( this.rxn.hashCode() << 2 );
        }
        public boolean equals(Object o) {
            if(!(o instanceof ExpandedHit)) {return false;}
            return this.rxn.equals( ((ExpandedHit)o).rxn ) && this.hit_fragments.equals( ((ExpandedHit)o).hit_fragments );
        }

    }

    /**
     * Contains all fragments. Fragments are sorted by the connectors they contain. Connectors are integers
     * of the connector atoms. I.e. a fragment containing a U and Np connector must contain 92 and 93.
     *
     */
    //Map<BitSet,Set<FragId>> fragments = new HashMap<>();
    //protected Map< Set<Integer> , Map<BitSet,Set<FragId>> > fragments_by_connectors = new HashMap<>();
    protected Map< BitSet , Map<BitSet,Set<FragId>> > fragments_by_connectors = new HashMap<>();


    /**
     * NOTE: the connector counts list is SORTED DESCENDING!
     */
    protected Map<List<Integer>,List<String>> rxns_by_connector_counts = new HashMap<>();


    protected Map<String,RxnId> rxns = new HashMap<>();

    public RxnId getRxn(String id) {
        return this.rxns.get(id);
    }

    public List<String> getRxnIds() {
        return new ArrayList<>(this.rxns.keySet());
    }

    /**
     * Returns all valid connector configurations in a set.
     * @return
     */
    //public Set<Set<Integer>> getFragmentsByConnectorsKeys() {
    public Set<BitSet> getFragmentsByConnectorsKeys() {
        return fragments_by_connectors.keySet();
    }

    /**
     * Returns empty hashset if connector set is not found.
     *
     * @param connectors
     * @return
     */
    //public Map<BitSet,Set<FragId>> getAllFragments_byConnectors(Set<Integer> connectors) {
    public Map<BitSet,Set<FragId>> getAllFragments_byConnectors(BitSet connectors) {
        Map<BitSet,Set<FragId>> all_frags = fragments_by_connectors.get(connectors);
        if(all_frags==null){all_frags = new HashMap<>();}
        return all_frags;
    }



    /**
     * Datastructure for connector-proximal analysis
     * Map:  "frag type" -> ( BitSet Tree containing the FP-BitSets  )
     */
    public Map<FragType,BitSetTree> connector_fps_sorted_by_fragtype = new HashMap<>();


    /**
     * Datastructure for connector-proximal analysis
     * Map:  "frag type" -> ( BitSet Tree containing the FP-BitSets  )
     */
    public Map<FragType,Map<BitSet,List<FragId>>> frags_sorted_by_connector_fp = new HashMap<>();

    /**
     * Datastructure for connector-proximal analysis
     * Map:  frag type -> ( fp_conn_bitset -> FastSubstrutureSearcher )
     */
    public Map<FragType,Map<BitSet,FastSubstructureSearcher>> substructure_searchers_sorted_by_connector_fp = new HashMap<>();

    /**
     * FragIds, by their String, this is used in the FastSubstructureSearcher as ID.
     * Map:  String -> FragId ; string is from FragId.toString function
     */
    public Map<String,FragId> frags_by_id = new ConcurrentHashMap<>();



    public Map<Integer,FragType> getFragTypes(String rxn) {
        return this.fragment_type_map.get(rxn);
    }

    //Map<String,BitSetTree> ffps_sorted_by_rxn_BT = new HashMap<>();
    //Map<String,List<BitSet>> ffps_sorted_by_rxn_L = new HashMap<>();
    //Map<String,Map<Integer,List<BitSet>>> ffps_sorted_by_rxn_and_frag_L = new HashMap<>();
    Map<String,Map<Integer,BitSetTree>> ffps_sorted_by_rxn_and_frag_BT = new HashMap<>();

    /**
     * These two datastructures contain all fragments / all fragment types sorted by rxn.
     */
    Map<String,Map<Integer,List<FragId>>> fragment_map = null;
    Map<String,Map<Integer,FragType>> fragment_type_map = null;

    Map<String,Map<Integer,Map<BitSet,List<FragId>>>> fragment_map_2 = null;

    //BitSetTree bst = null;
    //Map<Set<Integer>,BitSetTree> bsts_by_connectors = new HashMap<>();
    Map<BitSet,BitSetTree> bsts_by_connectors = new HashMap<>();

//    /**
//     * a fragment that contains u, np, pu will be in all bsts for: (u) , (np) , (pu) , (u,np) , (u,pu) , (np,pu) , (u,np,pu).
//     */
//    Map<BitSet,BitSetTree> bsts_by_contained_connector_subsets = new HashMap<>();
    Map<String,BitSetTree> bsts_labeled_by_rxn = new HashMap<>();


    /**
     * The key strings are generated via: encodeConnectorConfig( computeConnectorConfig( rxn ))
     */
    Map<String,List<String>> rxns_by_connector_config = new HashMap<>();


    /**
     * From https://stackoverflow.com/questions/27331175/java-bitset-comparison
     * @param lhs
     * @param rhs
     * @return
     */
    public static int compareBitsets(BitSet lhs, BitSet rhs) {
        if (lhs.equals(rhs)) return 0;
        BitSet xor = (BitSet)lhs.clone();
        xor.xor(rhs);
        int firstDifferent = xor.length()-1;
        if(firstDifferent==-1)
            return 0;
        return rhs.get(firstDifferent) ? 1 : -1;
    }

    /**
     *
     * @param connector_sets
     * @return
     */
    public static String encodeConnectorConfig(List<BitSet> connector_sets) {
        connector_sets.sort( (ba,bb) -> compareBitsets(ba,bb) );
        List<String> ssets = connector_sets.stream().map( bi -> createBitSetString(bi,MAX_CONNECTORS) ).collect(Collectors.toList());
        //List<String> ssets = connector_sets.stream().map( bi -> bi.cardinality() ).sorted( (ix,iy) -> -Integer.compare(ix,iy) ).map(ci -> ""+ci).collect(Collectors.toList());
        return String.join(";",ssets);
    }

    public static List<BitSet> decodeConnectorConfig(String cc) {
        String splits[] = cc.split(";");
        List<BitSet> result = new ArrayList<>();
        for( String si : splits ) {
            BitSet bsi = parseBitSetString(si);
            result.add(bsi);
        }
        return result;
    }

//    public static String encodeConnectorCounts() {
//    }
//
//    public static List<Integer> decodeConnectorCounts() {
//    }

//    public static String encodeConnectorConfig(List<Set<Integer>> connector_sets) {
//        List<String> sets = new ArrayList<>();
//
//        for(Set<Integer> ssi : connector_sets) {
//            List<Integer> csi = new ArrayList<>(ssi);
//            Collections.sort(csi);
//            String si = csi.stream().map( ci -> Integer.toString(ci) ).collect(Collectors.joining(","));
//            sets.add(si);
//        }
//        Collections.sort(sets);
//        return String.join(";",sets);
//    }
//    public static List<List<Integer>> decodeConnectorConfig(String cc) {
//        String splits[] = cc.split(";");
//        List<List<Integer>> result = new ArrayList<>();
//        for( String si : splits ) {
//            List<Integer> ri = new ArrayList<>();
//            String ints[] = si.split(",");
//            for(String int_i : ints) {
//                ri.add(Integer.parseInt(int_i));
//            }
//            Collections.sort(ri);
//            result.add(ri);
//        }
//        return result;
//    }

    //Map<String,Map<Set<Integer>,List<FragType>>> fragment_types_by_rxn_and_connector_config = new HashMap<>();
    Map<String,Map<BitSet,List<FragType>>> fragment_types_by_rxn_and_connector_config = new HashMap<>();

    public List<FragId> getSynthonSet(String rxn, int frag) {
        return this.fragment_map.get(rxn).get(frag);
    }

    public SynthonSpace() {

    }

    /**
     * !! Guards write access and manipulations of:
     *
     * this.fragments_by_connectors
     *
     */
    private Object lock_WriteA = new Serializable(){};



    /**
     * !! Guards write access and manipulations of:
     *
     * this.ffps_sorted_by_rxn_and_frag_BT
     * this.rxns_by_connector_counts
     * this.rxns
     */
    private Object lock_WriteB = new Serializable(){};


    /**
     * !! Guards write access and manipulations of:
     *
     * this.frags_sorted_by_connector_fp
     * this.connector_fps_sorted_by_fragtype
     * this.substructure_earchers_sorted_by_connector_fp
     *
     */
    private Object lock_WriteC = new Serializable(){};


    /**
     * Checks:
     *
     * 1. that all synthons in each set have the same connectors
     * 2. that the connectors are correct, i.e. that a they are starting at 92 and
     *    continue continuosly, e.g. that a 2 connector rxn has 92 and 93 connectors.
     * 3. that the bonds inbetween connectors are correct (each connector pair must have same bonds)
     *
     * NOTE! The list objects in the molecules parameter must either contain StereoMolecule objects
     * or String objects containing the idcode!!
     *
     * @param molecules
     * @return
     * @throws Exception
     */
    public static boolean ensureReactionIsOk(Map<Integer, List<Object>> molecules) throws Exception {

        int connector_counts[]                      = new int[200];
        Map<Integer,Integer>   connector_bond_types = new HashMap<>();

        IDCodeParser icp = new IDCodeParser();

        for(Integer ki : molecules.keySet()) {
            Set<Integer> connector_set = null;

            StereoMolecule template_i = null;
            if(molecules.get(ki).iterator().next() instanceof StereoMolecule) {
                template_i = (StereoMolecule) molecules.get(ki).iterator().next();
            }
            else if (molecules.get(ki).iterator().next() instanceof String) {
                template_i = new StereoMolecule();
                icp.parse(template_i, (String) molecules.get(ki).iterator().next() );
            }
            //molecules.get(ki).iterator().next();

            connector_set = computeConnectorSet(template_i);
            for(Integer ci : connector_set) { connector_counts[ci.intValue()]++;}

            StereoMolecule mi = new StereoMolecule();
            for(int mzi=0; mzi<molecules.get(ki).size(); mzi++) {
                StereoMolecule mi_pre = null;
                if( molecules.get(ki).get(mzi) == null ) {
                    System.out.println("Null molecule object supplied -> skip");
                    continue;
                }
                else if( molecules.get(ki).get(mzi) instanceof StereoMolecule) {
                    mi_pre = (StereoMolecule) molecules.get(ki).get(mzi);
                }
                else if(molecules.get(ki).get(mzi) instanceof String) {
                    mi_pre = new StereoMolecule();
                    icp.parse(mi_pre, (String) molecules.get(ki).get(mzi) );
                }
                else {
                    System.out.println("Incmpatible molecule object supplied -> skip");
                    continue;
                }
                mi = mi_pre;

            //for(StereoMolecule mi : molecules.get(ki)) {
                Set<Integer> csi = computeConnectorSet(mi);
                if(!connector_set.equals(csi)) {
                    throw new Exception("Connector Set " + ki + " has inconsistent connectors");
                }
                Map<Integer,Integer> connector_pos = computeConnectorPositions(mi);
                for(Integer cti : connector_pos.keySet()) {
                    int cpi = connector_pos.get(cti);
                    // check degree 1:
                    if(mi.getConnAtoms(cpi)!=1) { throw new Exception( "Connector with degree != 1" ); }
                    // check binding type:
                    //int cbt = mi.getBondType( mi.getBond( cpi , mi.getConnAtom(cpi,0) ) );
                    // check bond order:
                    int cb_order = mi.getBondOrder( mi.getBond( cpi , mi.getConnAtom(cpi,0) ) );

                    if( !connector_bond_types.containsKey( cti ) ) {
                        connector_bond_types.put(cti , cb_order );
                    }
                    else {
                        if( cb_order != connector_bond_types.get( cti )) {
                            throw new Exception("Inconsistent bond order detected! "+ cb_order +" != "+connector_bond_types.get(cti));
                        }
                    }
                }
            }
        }

        int num_connectors = connector_bond_types.keySet().size();

        // check that we have two connectors of each, and that they are starting at 92 and are ascending:
        List<Integer> connectors = connector_bond_types.keySet().stream().sorted().collect(Collectors.toList());

        if( connectors.get(0).intValue() != 92) { throw new Exception ("Smallest connector is not 92"); }
        for(int zi=0;zi<connectors.size();zi++) {
            if(connectors.get(zi) != 92+zi) {
                throw new Exception( "Connectors are not ascending from 92: "+connectors.toString() );
            }
        }

        for( int ci : connectors) {
            if(connector_counts[ci] != 2) { throw new Exception("Found "+connector_counts[ci]+" connectors of type "+ci+ " instead of 2"); }
        }

        return true;
    }

    /**
     * Supports synthon reactions with up to MAX_CONNECTOR different connectors.
     * The connector pairs must be (1.) U, (2.) Np , (3.) Pu, .. up to atomic no 99.
     *
     * This is now thread-safe! The only parallel write operations on the SynthonSpace object are:
     * 1. putting the fragments_by_connectors
     *
     * NOTE! The list objects in the molecules parameter must either contain StereoMolecule objects
     * or String objects containing the idcode!!
     *
     * @param rxn_id
     * @param molecules
     * @param idcode_to_identifier
     */
    public void addReaction(String rxn_id, Map<Integer, List<Object>> molecules , Map<String,String> idcode_to_identifier , CachedDescriptorProvider cdp)
    {
        //DescriptorHandlerLongPFP512 dh = new DescriptorHandlerLongPFP512();
        //FingerPrintGenerator fpgen = new FingerPrintGenerator();

        try {
            ensureReactionIsOk(molecules);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Reaction "+rxn_id+" -> detected problem, SKIP!");
            return;
        }

        Map<Integer,BitSetTree> all_bitsets_sorted = new HashMap<>();
        ArrayList<BitSet>         all_bitsets        = new ArrayList<>();


        IDCodeParser icp_a = new IDCodeParser();
        //Map<Integer,List<BitSetTree>> connector_fps_sorted_by_num_connectors = new HashMap<>();
        //Map<BitSet,List<FragId>>
        List<Integer> count_connectors   = new ArrayList<>(); // store the number of connectors for each synthon.
        List<BitSet>  connector_configs  = new ArrayList<>(); // store the connector config bitstrings for each synthon.

        for( Integer ki : molecules.keySet() ) {
            Map<BitSet,String> descriptors = new HashMap<>();
            List<BitSet> descriptors_i = new ArrayList<>();

            // This is for the connector proximal region analysis:
            Map<BitSet,List<FragId>> connector_proximal_sorted_fragments = new HashMap<>();
            int connector_region_size = CONNECTOR_REGION_SIZE; // max dist that is considered

            //for( StereoMolecule mi_pre : molecules.get(ki) ) {
            for(int mzi=0; mzi<molecules.get(ki).size(); mzi++) {
                StereoMolecule mi_pre = null;
                if( molecules.get(ki).get(mzi) == null ) {
                    System.out.println("Null molecule object supplied -> skip");
                    continue;
                }
                else if( molecules.get(ki).get(mzi) instanceof StereoMolecule) {
                    mi_pre = (StereoMolecule) molecules.get(ki).get(mzi);
                }
                else if(molecules.get(ki).get(mzi) instanceof String) {
                    mi_pre = new StereoMolecule();
                    icp_a.parse(mi_pre, (String) molecules.get(ki).get(mzi) );
                }
                else {
                    System.out.println("Incmpatible molecule object supplied -> skip");
                    continue;
                }

                String original_id_code = ""+new String( mi_pre.getIDCode() );
                //mi.setFragment(true);
                //mi = (new Canonizer(mi)).getCanMolecule();
                //StereoMolecule mi_2 = (new Canonizer(mi_pre)).getCanMolecule();
                //Canonizer cmi = new Canonizer(mi_pre);
                //cmi.
                mi_pre.ensureHelperArrays(Molecule.cHelperCIP);
                StereoMolecule mi = new StereoMolecule(mi_pre);
                mi.ensureHelperArrays(StereoMolecule.cHelperCIP);

                //BitSet bs = BitSet.valueOf( dh.createDescriptor(mi) );
                //BitSet bs = fpgen.getFingerprint(mi,BITS);

                //BitSet bs = getFP( getDescriptorHandler_PerThread() , mi );
                BitSet bs = null;
                if(cdp==null) {
                    bs = getFP( getDescriptorHandler_PerThread() , mi );
                }
                else{
                    bs = cdp.getFP_cached(mi);
                }

                // create the fp for the connector with non-unique connectors:
                StereoMolecule mi_non_unique_connis = new StereoMolecule(mi);
                mi_non_unique_connis.ensureHelperArrays(Molecule.cHelperNeighbours);
                for(int zi=0;zi<mi_non_unique_connis.getAtoms();zi++) {
                    if(mi_non_unique_connis.getAtomicNo(zi) >= 92 ) {
                        mi_non_unique_connis.setAtomicNo(zi,92);
                    }
                }
                //BitSet bs_non_unique = getFP( getDescriptorHandler_PerThread() , mi_non_unique_connis );
                BitSet bs_non_unique = null;
                if(cdp==null) {
                    bs_non_unique = getFP( getDescriptorHandler_PerThread() , mi_non_unique_connis );
                }
                else{
                    bs_non_unique = cdp.getFP_cached(mi_non_unique_connis);
                }

                descriptors.put( bs , ""+new String( mi.getIDCode() ) );
                descriptors_i.add(bs);
                all_bitsets.add(bs);

                FragId frag = new FragId( rxn_id,ki, mi , idcode_to_identifier.get( original_id_code ) , bs , bs_non_unique );
                this.frags_by_id.put(frag.toString(),frag);

                //Set<Integer> connectors = frag.getConnectors();
                BitSet connectors = frag.getConnectors();

                synchronized(lock_WriteA) {
                    if (!fragments_by_connectors.containsKey(connectors)) {
                        fragments_by_connectors.put(connectors, new HashMap<>());
                    }
                    Map<BitSet, Set<FragId>> fragments = getAllFragments_byConnectors(connectors);
                    if (!fragments.containsKey(bs)) {
                        fragments.put(bs, new HashSet<>());
                    }
                    //fragments.get(bs).add(new FragId( rxn_id,ki, mi.getIDCode() , idcode_to_identifier.get(mi.getIDCode() ) ));
                    //fragments.get(bs).add(new FragId( rxn_id,ki, mi , idcode_to_identifier.get( mi.getIDCode() ) ));
                    fragments.get(bs).add(new FragId(rxn_id, ki, mi, idcode_to_identifier.get(original_id_code), bs,bs_non_unique));
                }

                // hmm..
                //mi = null;


                // analyze the fragment connector regions:
                // CONVENTION: ALL CONNECTORS IN THE FRAGMENTS ARE JUST Uranium (i.e. atom no 92)!!

                StereoMolecule mi_connector_proximal = createConnectorProximalFragment(mi_pre,connector_region_size);
                //System.out.println("cut fragment idcode: "+mi_connector_proximal.getIDCode());

                //BitSet bci = this.getFP(mi_connector_proximal);
                BitSet bci = null;
                if(cdp==null) {
                    bci = getFP( getDescriptorHandler_PerThread() , mi_connector_proximal );
                }
                else{
                    bci = cdp.getFP_cached(mi_connector_proximal);
                }

                if(!connector_proximal_sorted_fragments.containsKey(bci)){connector_proximal_sorted_fragments.put(bci,new ArrayList<>());}
                //FragId frag_b = new FragId( rxn_id,ki, mi , idcode_to_identifier.get( original_id_code ) , bs );
                connector_proximal_sorted_fragments.get(bci).add(frag);

            }
            BitSetTree bst_x = BitSetTree.createTree( new HashSet<>(descriptors_i) ,BITS,BITTREE_BIN_SIZE);
            all_bitsets_sorted.put(ki, bst_x);

            // for connector-proximity analysis:
            // create connector fp bitset tree:
            BitSetTree conn_bst_x = BitSetTree.createTree( new HashSet<>(connector_proximal_sorted_fragments.keySet()) , BITS,BITTREE_BIN_SIZE );

            // then, to exploit the connector-proximity initial hits, we also need the fast substructure searcher for the different conn fps:
            Map<BitSet,FastSubstructureSearcher> sorted_substructure_searchers = new HashMap<>();
            for( BitSet bffi : connector_proximal_sorted_fragments.keySet()) {
                FastSubstructureSearcher fss = new FastSubstructureSearcher();
                //fss.setFP(this.mFP,this.BITS,true); // !! SET ON THE FLY, else all the descriptors take a lot of memory..
                fss.setFP(this.mFP,this.BITS,false); // !! SET ON THE FLY, else all the descriptors take a lot of memory..
                Map<String,String> idcodes_to_ids = new HashMap<>();
                connector_proximal_sorted_fragments.get(bffi).stream().forEach( fi -> idcodes_to_ids.put(fi.idcode,fi.toString()) );
                Map<String,BitSet> descriptors_map = new HashMap<>();
                connector_proximal_sorted_fragments.get(bffi).stream().forEach( fi -> descriptors_map.put(fi.idcode,fi.fp) );
                List<String> idcodes_i = new ArrayList<>(idcodes_to_ids.keySet());
                try {       // NOTE: TODO: using FFP_plus_1024 is not a good idea..
                    fss.initWithIDCodes(idcodes_i,idcodes_to_ids,descriptors_map,1,null,false);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                fss.clearCaches();
                sorted_substructure_searchers.put(bffi,fss);
            }

            synchronized (lock_WriteC) {
                frags_sorted_by_connector_fp.put(new FragType(rxn_id, ki), connector_proximal_sorted_fragments);
                connector_fps_sorted_by_fragtype.put(new FragType(rxn_id, ki),conn_bst_x);
                substructure_searchers_sorted_by_connector_fp.put(new FragType(rxn_id, ki),sorted_substructure_searchers);
            }
            System.out.println("Added frag type "+ (rxn_id+":"+ki) +" -> number of connector-proximity types: "+sorted_substructure_searchers.keySet().size());

            int num_connectors = connector_proximal_sorted_fragments.values().iterator().next().iterator().next().connectors.cardinality();
            count_connectors.add(num_connectors);

            BitSet connector_config_i = connector_proximal_sorted_fragments.values().iterator().next().iterator().next().connectors;
            connector_configs.add(connector_config_i);
//            if( !connector_fps_sorted_by_num_connectors.containsKey(num_connectors) ){connector_fps_sorted_by_num_connectors.put(num_connectors,new ArrayList<>());}
//            connector_fps_sorted_by_num_connectors.get(num_connectors).add(bst_x);

        }

        //ffps_sorted_by_rxn_L.put(rxn_id,all_bitsets);
        //ffps_sorted_by_rxn_and_frag_L.put(rxn_id,all_bitsets_sorted);

        //if(USE_BITTREES_FOR_RXNSETS) {
        //ffps_sorted_by_rxn_BT.put(rxn_id, BitSetTree.createTree( new HashSet<>(all_bitsets), BITS, BITTREE_BIN_SIZE));

        count_connectors.sort((x,y) -> -Integer.compare(x,y));

        int num_connectors = count_connectors.stream().mapToInt( ci -> ci.intValue() ).sum() / 2; // can be computed by dividing the total number of connectors by two.

        synchronized(lock_WriteB) {
            ffps_sorted_by_rxn_and_frag_BT.put(rxn_id, all_bitsets_sorted);
            if(!rxns_by_connector_counts.containsKey(count_connectors)) { rxns_by_connector_counts.put(count_connectors,new ArrayList<>()); }
            rxns_by_connector_counts.get(count_connectors).add(rxn_id);
            rxns.put(rxn_id,new RxnId(rxn_id,num_connectors,count_connectors,connector_configs));
        }
    }

    /**
     * Creates the connector proximal fragments.
     *
     * NOTE!! In these fragments, ALL CONNECTORS, i.e. all
     * atoms with atomic no inbetween 92 and 99 WILL BE
     * REPLACED BY uranium, i.e. by 92 !!
     *
     * NOTE!! In these fragments, all NARROWING query features will
     *        be removed. (the rationale is, NARROWING query features
     *        can only RESTRICT the set of compounds, so to be
     *        safe we do not want to have this in a pruning
     *        step)
     *
     * @param mi_pre
     * @param connector_region_size
     * @return
     */
    public static StereoMolecule createConnectorProximalFragment(StereoMolecule mi_pre, int connector_region_size) {
        StereoMolecule mi_conn = new StereoMolecule(mi_pre);

        //mi_conn.removeQueryFeatures(); // !! Remove query features
        // REMOVE ONLY NARROWING QUERY FEATURES:
        QueryFeatureHelper.removeNarrowingQueryFeatures(mi_conn);


        mi_conn.ensureHelperArrays(Molecule.cHelperCIP);
        //mi_conn.getAtoms();
        // create the connector-near fragment:
        // 1. set all connector atoms to uranium, and store the positions:
        List<Integer> connector_positions = new ArrayList<>();
        for(int zi=0;zi<mi_conn.getAtoms();zi++) {
            int an = mi_conn.getAtomicNo(zi);
            if(an>=92&&an<92+MAX_CONNECTORS) {
                connector_positions.add(zi);
                mi_conn.setAtomicNo(zi,92);
            }
        }
        // 2. cut out the connector region:
        boolean keep_atoms[] = new boolean[mi_conn.getAtoms()];
        for(int zi=0;zi<mi_conn.getAtoms();zi++) {
            for(int ci : connector_positions) {
                // NOTE! this returns -1 if no path is found within connector_region_size (I think..)
                int path_length = mi_conn.getPathLength(ci,zi,connector_region_size,null);
                if ( path_length>=0 ){
                    if(path_length<=connector_region_size) {
                        keep_atoms[zi] = true;
                    }
                }
            }
        }
        StereoMolecule mi_cut = new StereoMolecule();
        mi_conn.copyMoleculeByAtoms(mi_cut,keep_atoms,true,null);
        mi_cut.ensureHelperArrays(Molecule.cHelperCIP);

        if(false) {
            System.out.println("CPF: " + mi_cut.getIDCode());
        }

        return mi_cut;
    }

    public void reinitBitTree() {
        //for(Set<Integer> conni : fragments_by_connectors.keySet()) {
        synchronized(lock_WriteA){
            synchronized (lock_WriteB) {
                synchronized (lock_WriteC) { // all this synch. should not be needed?..
                    for(BitSet conni : fragments_by_connectors.keySet()) {
                        BitSetTree bst = BitSetTree.createTree(fragments_by_connectors.get(conni).keySet(), BITS, BITTREE_BIN_SIZE);
                        bsts_by_connectors.put(conni, bst);
                    }
                    for(String rxi : rxns.keySet()) {
                        Set<BitSet> all_ffps_labeled_for_rxn = new HashSet<>();
                        for( int fti : fragment_map.get(rxi).keySet()) {
                            all_ffps_labeled_for_rxn.addAll( fragment_map.get(rxi).get(fti).stream().map( fi -> fi.fp ).collect(Collectors.toList()));
                        }
                        BitSetTree bst = BitSetTree.createTree(all_ffps_labeled_for_rxn, BITS, BITTREE_BIN_SIZE);
                        bsts_labeled_by_rxn.put(rxi,bst);
                    }
//                    // and create bsts for all subsets..
//                    // 1. determine all bits that we have as connectors
//                    BitSet all_bits = fragments_by_connectors.keySet().stream().reduce( (x,y) -> {BitSet z = ((BitSet)x.clone()); z.or(y); return z;} ).get();
//                    // 2. determine all subsets of these bits..
//                    // 2.1 check that we have "consecutive" connectors
//                    int all_bits_cardinality = all_bits.cardinality();
//                    for(int zi=0;zi<all_bits_cardinality;zi++) { if(!all_bits.get(zi)){throw new RuntimeException("Non consecutive connectors? This is not allowed..");} }
                }
            }
        }

    }

    /**
     *
     * Assumes that the following fields exist and inits all maps from these fields:
     * fragment_type_map
     * fragments_by_connectors
     *
     */
    public void reinitHelperMaps() {
        // also init rxns_by_connector_config
        rxns_by_connector_config = new HashMap<>();
        for(String rxn : this.fragment_type_map.keySet()) {
            //List<Set<Integer>> rxn_connector_config = new ArrayList<>();
//            List<BitSet> rxn_connector_config = new ArrayList<>();
//            for( FragType fti : this.fragment_type_map.get(rxn).values() ) {
//                rxn_connector_config.add( getConnectorSetForFragType(fti) );
//            }

            List<BitSet> rxn_connector_config = computeRxnConnectorConfig(rxn);

            if(!rxns_by_connector_config.containsKey(encodeConnectorConfig(rxn_connector_config))) {
                rxns_by_connector_config.put( encodeConnectorConfig(rxn_connector_config) , new ArrayList<>());
            }

            // todo: here we could check, if the synthon sets make sense for this rxn or if maybe one is missing..
            rxns_by_connector_config.get( encodeConnectorConfig(rxn_connector_config) ).add(rxn);
        }

        // also init fragment_types_by_rxn_and_connector_config
        fragment_types_by_rxn_and_connector_config = new HashMap<>();
        for(String rxn : this.fragment_type_map.keySet()) {
            //Map<Set<Integer>,List<FragType>> rxn_i_map = new HashMap<>();
            Map<BitSet,List<FragType>> rxn_i_map = new HashMap<>();
            for( FragType fti : this.fragment_type_map.get(rxn).values() ) {
                //Set<Integer> cci = getConnectorSetForFragType(fti);
                BitSet cci = getConnectorSetForFragType(fti);
                if(!rxn_i_map.containsKey(cci)){rxn_i_map.put(cci,new ArrayList<>());}
                rxn_i_map.get(getConnectorSetForFragType(fti)).add(fti);
            }
            fragment_types_by_rxn_and_connector_config.put(rxn,rxn_i_map);
        }

        // also init fragment_map and fragment_map_2
        fragment_map   = new HashMap<>();
        fragment_map_2 = new HashMap<>();

        //for(Set<Integer> si : fragments_by_connectors.keySet()) {
        for(BitSet si : fragments_by_connectors.keySet()) {
            for( Set<FragId> fid_set : fragments_by_connectors.get(si).values() ) {
                for( FragId fid : fid_set ) {

                    if( !fragment_map.containsKey(fid.rxn_id) ) { fragment_map.put(fid.rxn_id,new HashMap<>()); }
                    if( !fragment_map.get(fid.rxn_id).containsKey(fid.frag) ) { fragment_map.get(fid.rxn_id).put(fid.frag,new ArrayList()); }

                    if( !fragment_map_2.containsKey(fid.rxn_id) ) { fragment_map_2.put(fid.rxn_id,new HashMap<>()); }
                    if( !fragment_map_2.get(fid.rxn_id).containsKey(fid.frag) ) { fragment_map_2.get(fid.rxn_id).put(fid.frag,new HashMap<>()); }

                    if( !fragment_map_2.get(fid.rxn_id).get(fid.frag).containsKey(fid.fp) ) { fragment_map_2.get(fid.rxn_id).get(fid.frag).put(fid.fp,new ArrayList<>()); }

                    fragment_map.get(fid.rxn_id).get(fid.frag).add(fid);
                    fragment_map_2.get(fid.rxn_id).get(fid.frag).get(fid.fp).add(fid);
                }
            }
        }
    }

    public List<BitSet> computeRxnConnectorConfig(String rxn) {
        List<BitSet> rxn_connector_config = new ArrayList<>();
        for( FragType fti : this.fragment_type_map.get(rxn).values() ) {
            rxn_connector_config.add( getConnectorSetForFragType(fti) );
        }
        return rxn_connector_config;
    }

    public List<InitialHit> screen(CachedDescriptorProvider cdh, StereoMolecule mol, int num_splits , Set<String> rxns_to_omit , int max_fragments, int max_initial_hits) {
        return screen(cdh,mol,num_splits,rxns_to_omit,max_fragments,max_initial_hits,4);
    }

    /**
     *
     *
     * @param mol
     * @param num_splits
     * @param rxns_to_omit
     * @param max_initial_hits if we reach this number of initial hits then we stop looking for more
     * @param max_fragments maximum number of synthon sets in a single synthon reaction
     * @return
     */
    public List<InitialHit> screen(CachedDescriptorProvider cdh, StereoMolecule mol, int num_splits , Set<String> rxns_to_omit , int max_fragments, int max_initial_hits, int threads) {

        // change this if you have a more complicated virtual space!!
        final int MAX_FRAGMENTS = max_fragments;

        //CachedDescriptorProvider cdh = new CachedDescriptorProvider(this);

        System.out.println("[SS_SCREEN] Start screening, splits = "+num_splits);
        Map<FragType,List<InitialHit>> initial_hits = findInitialHits(this,cdh,mol,num_splits,MAX_FRAGMENTS, max_initial_hits, threads);

        //List<Hit> discovered_hits = new ArrayList<>();
        List<InitialHit> discovered_hits = new ArrayList<>();

        // sort initial hits by rxn:
        Map<String,Map<FragType,List<InitialHit>>> sorted_initial_hits = new HashMap<>();
        for(FragType fti : initial_hits.keySet()) {
            String rxn = fti.rxn_id;
            if(!sorted_initial_hits.containsKey(rxn)){sorted_initial_hits.put(rxn,new HashMap<FragType,List<InitialHit>>());}
            sorted_initial_hits.get(rxn).put(fti,initial_hits.get(fti));
        }

        for(String rxn : sorted_initial_hits.keySet()) {

            //for(FragType ft : initial_hits.keySet()) {
            for (FragType ft : sorted_initial_hits.get(rxn).keySet()) {
                if (discovered_hits.size() >= max_initial_hits) {
                    System.out.println("BREAK EARLY!! SynthonSpace::screen(..) -> reached max number of initial hits: " + max_initial_hits);
                    break;
                }
                //System.out.println("Analyze initial hits of frag type " + ft.toString() );

                List<InitialHit> candidates = initial_hits.get(ft);
                for (int zi = 0; zi < candidates.size(); zi++) {
                    // check if we have enough..
                    if (discovered_hits.size() >= max_initial_hits) {
                        break;
                    }
                    InitialHit hi = candidates.get(zi);

                    // check if rxn is compatible with connector configuration is NOT necessary, we already do this when computing the initial hits!

                    try {
                        if (rxns_to_omit.contains(hi.fragment_type.rxn_id)) {
                            //System.out.println("Skip, old one "+hi.fragment_type.rxn_id);
                            continue;
                        }
                    }
                    catch(Exception ex) {
                        System.out.println("[ERR] wut?? -> we cannot really end up here?!?");
                        continue;
                    }

                    // sort remaining fps by cardinality :) (descending)
                    List<Map.Entry<FragType, BitSet>> remaining = new ArrayList<>(hi.remaining_fps.entrySet());
                    remaining.sort((hta, htb) -> -Integer.compare(hta.getValue().cardinality(), htb.getValue().cardinality()));

                    // First determine all possible fragment to bitset mappings,
                    // then verify subsets that are ok (and verify that all are ok).
                    boolean still_ok = true;
                    for (Map.Entry<FragType, BitSet> fi : remaining) {
                        BitSetTree.Node result_node[] = new BitSetTree.Node[1];
                        //Set<Integer> connectors_fi = fi.getKey().;
                        BitSet       fp_fi         = fi.getValue();
                        // find fragment type:
                        FragType fti = fi.getKey(); //fragment_types_by_rxn_and_connector_config.get(rxn).get(connectors_fi);
                        //boolean found = this.ffps_sorted_by_rxn_and_frag_BT.get(rxn).get(fi).testSubset(q,result_node);
                        boolean found = this.ffps_sorted_by_rxn_and_frag_BT.get(rxn).get( fti.frag ).testSubset(fp_fi, result_node);
                        still_ok &= found;
                        if(!still_ok){break;}
                    }

                    if(still_ok) {
                        //PrimaryHit pi = new SynthonSimilaritySpace.PrimaryHit()
                        discovered_hits.add( hi );
                    }
                }
            }
        }
        return discovered_hits;
    }


    /**
     * This expands the hits, i.e. it collects all possible fragments given by the BitSets for the
     * non-primary fragments.
     *
     * The compatible fragments are computed in two steps: first, the PathFP bits are considered, then
     * a matching is performed against the actual fragment.
     *
     */
    public Map<InitialHit,List<ExpandedHit>> expandHits_New(List<InitialHit> hits, int max_per_frag, int threads) {

        // We create one pool of fixed size, and then we submit tasks and later also nested tasks to this pool.
        ExecutorService main_pool = Executors.newFixedThreadPool(threads);

        // it may happen that we perform the same substructure searches multiple times.
        Map<String,Boolean> substructureTestCache = new ConcurrentHashMap<>();

        Map<InitialHit,Map<FragType, List<FragId>>> all_found_fragments = new ConcurrentHashMap<>();
        Map<InitialHit,List<ExpandedHit>> result = new ConcurrentHashMap<>();

        List<Future> tasks_supersets = new ArrayList<>();
        for(InitialHit hi : hits) {
            Runnable hi_r = new Runnable(){
                @Override
                public void run() {
                    expandHits_New_findSupersets(hi, all_found_fragments);
                }
            };
            tasks_supersets.add( main_pool.submit(hi_r) );
        }

        // wait for all tasks to finish:
        if(logLevel_hitExpansion>0){
            System.out.println("[EXH] : Find Supersets:\n");
        }
        int cnt_tsup = 0;
        for(Future fi : tasks_supersets) {
            try {
                fi.get();
                cnt_tsup++;
                if(logLevel_hitExpansion>0) {
                    System.out.print(".");
                    if (cnt_tsup % 60 == 0) {
                        System.out.print("\n");
                    }
                }
            }
            catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        if(logLevel_hitExpansion>0){
            System.out.println("[EXH] : Find Supersets: done!\n");
        }


        for(InitialHit hi : hits) {

            Map<FragType, List<FragId>> found_fragments = all_found_fragments.get(hi);

            if(logLevel_hitExpansion>0){
                System.out.println("[EXH] : Hit "+hi.toString()+" second stage -> Candidates: "+found_fragments.values().stream().mapToInt( vi -> vi.size() ).sum());
            }

            // now filter by matching..
            Map<Integer, Set<FragId>> matching_fragments = new HashMap<>();

            // add initial fragment
            Set<FragId> initial_frag_set = new HashSet<>();
            initial_frag_set.add(hi.primary_fragment_id);
            matching_fragments.put(hi.primary_fragment_id.frag,initial_frag_set);

            IDCodeParser idp = new IDCodeParser();
            if (logLevel_hitExpansion > 1) {
                System.out.println("Match candidate fragments: ");
            }

            // we queue tasks in the inner loop..
            List<Future> tasks_ss = new ArrayList<>();

            for (FragType fti : found_fragments.keySet()) {

//                if (ki.equals(hi.primary_fragment)) {
//                    //matching_fragments.put(hi.primary_fragment,new HashSet<>(found_fragments.get(ki)));
//                    HashSet<FragId> primary_fragment_set = new HashSet<>();
//                    primary_fragment_set.add(hi.primary_fragment_id);
//                    matching_fragments.put(hi.primary_fragment, primary_fragment_set);
//                    continue;
//                }

//                if(ki.equals( hi.primary_fragment)) {
//                    System.out.println("Primary fragment!");
//                }

                if (logLevel_hitExpansion > 1) {
                    System.out.println("Synthon " + fti.toString());
                }
                if (logLevel_hitExpansion > 1) {
                    System.out.println("Perform " + found_fragments.get(fti).size() + " substructure searches for FragType " + hi.fragment_type.rxn_id + " F:" + fti.toString());
                }


                List<FragId> candidate_frags = found_fragments.get(fti);

                // match
                StereoMolecule target_frag = hi.remaining_frags.get(fti);
                target_frag.setFragment(true);

                String target_idcode = target_frag.getIDCode();

                //SSSearcher sss = new SSSearcher();
                //sss.setFragment(target_frag);
                SSSearcherProvider ssvp = new SSSearcherProvider(target_frag);

                //List<FragId> matching_frags_i = new ArrayList<>();
                Set<FragId> matching_frags_i = Collections.newSetFromMap(new ConcurrentHashMap<>());

                //StereoMolecule lib_frag = new StereoMolecule();

                for (FragId fid : found_fragments.get(fti)) {
                    //expandHits_New_matchFragment(substructureTestCache,sss,fid,matching_frags_i,target_idcode,max_per_frag);
                    //expandHits_New_matchFragment(substructureTestCache,ssvp,fid,matching_frags_i,target_idcode,max_per_frag);
                    Runnable r_matching = new Runnable() {
                        @Override
                        public void run() {
                            expandHits_New_matchFragment(substructureTestCache,ssvp,fid,matching_frags_i,target_idcode,max_per_frag);
                        }
                    };
                    tasks_ss.add( main_pool.submit(r_matching) );
                }
                matching_fragments.put(fti.frag, matching_frags_i);

//                if (matching_frags_i.isEmpty()) {
//                    // in this case we can continue, this Hit does not work..
//                    break;
//                }
            }

            // wait for pool:
            // wait for all tasks to finish:
            if(logLevel_hitExpansion>0){
                System.out.println("[EXH] : Match fragment:\n");
            }
            int cnt_mf = 0;
            for(Future fi : tasks_ss) {
                try {
                    fi.get();
                    cnt_mf++;
                    if(logLevel_hitExpansion>0) {
                        System.out.print(".");
                        if (cnt_mf % 60 == 0) {
                            System.out.print("\n");
                        }
                    }
                }
                catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                }
            }
            if(logLevel_hitExpansion>0){
                System.out.println("[EXH] : fragment: done!\n");
            }

//            // apply max per frag..
//            for(Integer ki : found_fragments.keySet()) {
//                Collections.shuffle(found_fragments.get(ki),r);
//                found_fragments.put(ki,found_fragments.get(ki).subList(0,Math.min(found_fragments.get(ki).size(),max_per_frag)));
//            }

            // add combinations..
            // quick hack (without recursion)
            //List<Map<Integer,FragId>> frags_combined = combine3(found_fragments);
            boolean is_empty = false;
            Map<Integer, List<FragId>> matching_fragments_l = new HashMap<>();
            for (Integer ki : matching_fragments.keySet()) {
                matching_fragments_l.put(ki, new ArrayList<>(matching_fragments.get(ki)));
                if (matching_fragments.get(ki).isEmpty()) {
                    is_empty = true;
                }
            }

            if (!is_empty) {
                List<Map<Integer, FragId>> frags_combined = combine3(matching_fragments_l);
                List<ExpandedHit> exp_hits = new ArrayList<>(frags_combined.stream().map(mi -> new ExpandedHit(hi.fragment_type.rxn_id, mi)).collect(Collectors.toList()));
                result.put(hi, exp_hits);
            }
        }

        // required?
        main_pool.shutdown();

        return result;
    }

    private void expandHits_New_findSupersets(InitialHit hi, Map<InitialHit,Map<FragType, List<FragId>>> all_found_fragments) {
        Map<FragType, List<FragId>> found_fragments = new HashMap<>();
        List<FragType> missing_frags = new ArrayList<>(hi.remaining_frags.keySet());

        long combis = 1;
        for (FragType fti : missing_frags) {

            List<BitSet> fp_supersets = new ArrayList<>();
            this.ffps_sorted_by_rxn_and_frag_BT.get(fti.rxn_id).get(fti.frag).root.collectSuperSets(hi.remaining_fps.get(fti), fp_supersets);

            FragType ft = fti;//this.fragment_type_map.get(hi.rxn).get(mfi);
            Map<BitSet, Set<FragId>> fragments_i = this.fragments_by_connectors.get(getConnectorSetForFragType(ft));

            // now resolve the fp bitset to the actual fragment:
            List<FragId> frags = new ArrayList<>();
            for (BitSet bsi : fp_supersets) {
                //List<FragId> frags_i = this.fragments.get(bsi).stream().filter( fi -> fi.frag==mfi && fi.rxn_id.equals(hi.rxn) ).collect(Collectors.toList());
                List<FragId> frags_i = new ArrayList<>();
                try {
                    frags_i = fragments_i.get(bsi).stream().filter(fi -> fi.frag == fti.frag && fi.rxn_id.equals(fti.rxn_id)).collect(Collectors.toList());
                } catch (Exception ex) {
                    System.out.println("[Exception] Could not find frag for frag fp.. this is strange..");
                }

                frags.addAll(frags_i);
            }

            combis *= frags.size();
            found_fragments.put(fti, frags);
        }
        all_found_fragments.put(hi,found_fragments);
    }


    /**
     * Thread-safe sss provider
     */
    static class SSSearcherProvider {
        StereoMolecule mM                             = null;
        Map<Thread,SSSearcher> mThreadSafeSSSearchers = new HashMap<>();

        public SSSearcherProvider(StereoMolecule target_frag) {
            this.mM = target_frag;
        }

        public synchronized SSSearcher getSSS() {
            Thread tci = Thread.currentThread();
            if(mThreadSafeSSSearchers.containsKey(tci)) { return mThreadSafeSSSearchers.get(tci); }
            SSSearcher ssi = new SSSearcher();
            ssi.setFragment(mM);
            mThreadSafeSSSearchers.put(tci,ssi);
            return ssi;
        }
    }

    private static void expandHits_New_matchFragment(Map<String,Boolean> substructureTestCache, SSSearcherProvider sssvp, FragId fid, Set<FragId> matching_frags_i, String target_idcode, int max_per_frag) {

        SSSearcher sss = sssvp.getSSS();

        IDCodeParser idp = new IDCodeParser();
        StereoMolecule lib_frag = new StereoMolecule();

        if (matching_frags_i.size() >= max_per_frag) {
            System.out.println("BREAK EARLY!! Reached " + max_per_frag + " extension hits for fragment " + fid.toString());
            //continue;
            return;
        }

        String match_key = fid.idcode + "_sstest_in_" + target_idcode;

        boolean is_subset = false;
        if (substructureTestCache.containsKey(match_key)) {
            is_subset = substructureTestCache.get(match_key);
            idp.parse(lib_frag, fid.idcode);
        } else {
            idp.parse(lib_frag, fid.idcode);
            lib_frag.setFragment(false);
            lib_frag.ensureHelperArrays(StereoMolecule.cHelperCIP);
            sss.setMolecule(lib_frag);
            is_subset = sss.isFragmentInMolecule();
            substructureTestCache.put(match_key, is_subset);
            if (is_subset) {
                System.out.println("ok..");
            }
        }
        if (logLevel_hitExpansion > 1) {
            System.out.println("Match fragment / synthon: " +
                    //Hyperspace.idcodeToSmiles(target_frag.getIDCode())+"."+Hyperspace.idcodeToSmiles((lib_frag.getIDCode())));
                    //Hyperspace_old.idcodeToSmiles(target_frag.getIDCode()) + "." + Hyperspace_old.idcodeToSmiles((fid.idcode)));
                    HyperspaceUtils.idcodeToSmiles(target_idcode) + "." + HyperspaceUtils.idcodeToSmiles((fid.idcode)));

            System.out.println(SynthonAssembler.idcodeToSmiles(fid.idcode) + "." + SynthonAssembler.idcodeToSmiles(fid.idcode));
            //DebugOutput.plotMolecules("Match"+fid.idcode,new String[]{target_frag.getIDCode(),lib_frag.getIDCode()},1,2);
//                        try {
//                            Thread.sleep(6000);
//                        } catch (InterruptedException e) {
//                            e.printStackTrace();
//                        }
        }

        if (is_subset) {
            //System.out.println("Frag OK!");
            matching_frags_i.add(fid);
        } else {
            //System.out.println("Frag REJECTED!");
        }


        if (logLevel_hitExpansion > 1) {
            System.out.println((is_subset) ? "MATCH" : "no match");
            if (is_subset) {
                System.out.println(target_idcode + "." + HyperspaceUtils.idcodeToSmiles(fid.idcode));
//                            DebugOutput.plotMolecules("SS:" + ((is_subset) ? "MATCH" : "no match"), new String[]{target_frag.getIDCode(), fid.idcode}, 1, 2);
//                            try {
//                                Thread.sleep(6000);
//                            } catch (InterruptedException e) {
//                                e.printStackTrace();
//                            }
            }
        }
        if (logLevel_hitExpansion > 1) {
            System.out.println(target_idcode + " , " + HyperspaceUtils.idcodeToSmiles(lib_frag.getIDCode()) + "," + ((is_subset) ? 1 : 0));
        }
    }

    private Integer findInverse(Map<Integer,BitSet> d , BitSet t ) {
        Map.Entry<Integer,BitSet> entry = d.entrySet().stream().filter(ei -> ei.getValue().equals(t) ).findAny().get();
        return entry.getKey();
    }

    // hack (although, maybe faster than recursion?)
    private static List<Map<Integer,FragId>> combine3(Map<Integer,List<FragId>> f) {
        ArrayList<Map<Integer,FragId>> combined = new ArrayList<>();
        if(f.keySet().size()==1) {
            Integer key_a = f.keySet().iterator().next();
            for(int zi=0;zi<f.get(key_a).size();zi++) {
                Map<Integer, FragId> ma = new HashMap<>();
                ma.put(key_a, f.get(key_a).get(zi));
                combined.add(ma);
            }
        }
        if(f.keySet().size()==2) {
            Iterator<Integer> it = f.keySet().iterator();
            Integer key_a = it.next();
            Integer key_b = it.next();
            for(int zi=0;zi<f.get(key_a).size();zi++) {
                for(int zj=0;zj<f.get(key_b).size();zj++) {
                    Map<Integer, FragId> ma = new HashMap<>();
                    ma.put(key_a, f.get(key_a).get(zi));
                    ma.put(key_b, f.get(key_b).get(zj));
                    combined.add(ma);
                }
            }
        }
        if(f.keySet().size()==3) {
            Iterator<Integer> it = f.keySet().iterator();
            Integer key_a = it.next();
            Integer key_b = it.next();
            Integer key_c = it.next();
            for(int zi=0;zi<f.get(key_a).size();zi++) {
                for(int zj=0;zj<f.get(key_b).size();zj++) {
                    for(int zk=0;zk<f.get(key_c).size();zk++) {
                        Map<Integer, FragId> ma = new HashMap<>();
                        ma.put(key_a, f.get(key_a).get(zi));
                        ma.put(key_b, f.get(key_b).get(zj));
                        ma.put(key_c, f.get(key_c).get(zk));
                        combined.add(ma);
                    }
                }
            }
        }
        return combined;
    }


    /**
     * Here we store the idcodes of the fragments resulting from a cut pattern.
     */
    //private Map<Set<Integer>,List<String>> cut_fragments = new HashMap<>();
    private Map<String,List<String>> cut_fragments = new HashMap<>();

//    /**
//     * Here we store the split results, based on the
//     */
//    private Map<Integer,List<String>> split_results = new HashMap<>();

    /**
     *
     * @param mol
     * @param num_splits
     * @param max_fragments if we cut the molecule into more than max_fragments we do not consider the cutset further.
     * @return
     */
    public static Map<FragType,List<InitialHit>> findInitialHits(SynthonSpace space, CachedDescriptorProvider cdh , StereoMolecule mol, int num_splits, int max_fragments, int max_hits, int threads) {

        List<InitialHit> initial_hits                      = Collections.synchronizedList(new ArrayList<>() );
        Map<FragType,List<InitialHit>> initial_hits_sorted = new ConcurrentHashMap<>();

        int nb = mol.getBonds();

        List<int[]> splits  = null;
        if(num_splits==0) {
            splits = new ArrayList<>();
            splits.add(new int[0]);
        }
        else {
            List<int[]> combi_list = CombinationGenerator.getAllOutOf(nb, num_splits);
            // !! returns null if b > a..
            if(combi_list==null) { return new HashMap<>(); }
            splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());
        }

        ExecutorService main_pool = Executors.newFixedThreadPool(threads);

        List<Future> tasks_initialhits = new ArrayList<>();

        for(int[] si : splits) {
            final int num_connectors = num_splits+1; // this we should also fix at the beginning..
            //findInitialHits_forSplitPattern(space,cdh,mol,si,num_connectors,max_fragments,max_hits,initial_hits,initial_hits_sorted);

            Runnable r_initialhits = new Runnable() {
                @Override
                public void run() {
                    findInitialHits_forSplitPattern(space,cdh,mol,si,num_connectors,max_fragments,max_hits,initial_hits,initial_hits_sorted);
                }
            };
            tasks_initialhits.add(main_pool.submit(r_initialhits));
        }

        // wait for pool:
        if(logLevel_findCandidates>0) {
            System.out.println("[FCH] Start screening for candidates:+\n");
        }
        int cnt_fc = 0;
        for(Future f_i : tasks_initialhits) {
            try {
                f_i.get();
                if(logLevel_findCandidates>0) {
                    System.out.print(".");
                    if ( (++cnt_fc) % 60 == 0) {
                        System.out.print("\n");
                    }
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        if(logLevel_findCandidates>0) {
            System.out.println("\n[FCH] Start screening for candidates: Done!");
        }

        // required?
        main_pool.shutdown();

        return initial_hits_sorted;
    }


    public static List<CombinatorialHit> findExpandedHits_withConnProximityMatching(SynthonSpace space, CachedDescriptorProvider cdh ,
                                                                                    StereoMolecule mol, int num_splits, int max_fragments,
                                                                                    Set<String> rxns_to_omit ,
                                                                                    int max_hits, int threads
                                                                                    )
    {
        return findExpandedHits_withConnProximityMatching(space,cdh,mol,num_splits,max_fragments,
                rxns_to_omit,max_hits,threads,null);
    }

    /**
     *
     *
     * @param space
     * @param cdh
     * @param mol
     * @param num_splits
     * @param max_fragments
     * @param rxns_to_omit
     * @param max_hits maximum number of combinatorial hits that should be examined
     * @param threads
     * @param output_statistics, can be null, then no statistics are returned
     * @return
     */
    public static List<CombinatorialHit> findExpandedHits_withConnProximityMatching(SynthonSpace space, CachedDescriptorProvider cdh ,
                                                                  StereoMolecule mol, int num_splits, int max_fragments,
                                                                  Set<String> rxns_to_omit ,
                                                                  int max_hits, int threads,
                                                                  List<SplitPatternWithConnectorProximityPruningStatistics> output_statistics)  {

       List<CombinatorialHit> expanded_hits                      = Collections.synchronizedList(new ArrayList<>() );
       //Map<FragType,List<InitialHit>> initial_hits_sorted   = new ConcurrentHashMap<>();

       int nb = mol.getBonds();

       List<int[]> splits  = null;
       if(num_splits==0) {
           splits = new ArrayList<>();
           splits.add(new int[0]);
       }
       else {
           List<int[]> combi_list = CombinationGenerator.getAllOutOf(nb, num_splits);
           // !! returns null if b > a..
           if(combi_list==null) { return new ArrayList<>(); }
           splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());
       }

       ExecutorService main_pool = Executors.newFixedThreadPool(threads);

       List<Future> tasks_initialhits = new ArrayList<>();

       for(int[] si : splits) {
           final int num_connectors = num_splits+1; // this we should also fix at the beginning..
           //findInitialHits_forSplitPattern(space,cdh,mol,si,num_connectors,max_fragments,max_hits,initial_hits,initial_hits_sorted);

           Runnable r_initialhits = new Runnable() {
               @Override
               public void run() {
                   findExpandedHits_forSplitPattern_withConnProximityMatching(space,cdh,mol,si,num_connectors,max_fragments,max_hits,rxns_to_omit,expanded_hits,output_statistics);
               }
           };
           tasks_initialhits.add(main_pool.submit(r_initialhits));
       }

       // wait for pool:
       if(logLevel_findCandidates>0) {
           System.out.println("[FCH] Start screening for candidates:+\n");
       }
       int cnt_fc = 0;
       for(Future f_i : tasks_initialhits) {
           try {
               f_i.get();
               if(logLevel_findCandidates>0) {
                   //System.out.print(".");
                   if ( (++cnt_fc) % 60 == 0) {
                   //    System.out.print("\n");
                   }
               }
           } catch (InterruptedException e) {
               e.printStackTrace();
           } catch (ExecutionException e) {
               e.printStackTrace();
           }
       }
       if(logLevel_findCandidates>0) {
           System.out.println("\n[FCH] Start screening for candidates: Done!");
       }

       // required?
       main_pool.shutdown();


       // if we used output_statistics, print some statistics:
       if(true) {
           if(output_statistics!=null) {
               System.out.println("Search complete!"+" [Splits="+num_splits+"]"+"Query= "+mol.getIDCode());
               System.out.println("Statistics:");
               System.out.println("Splits: "+output_statistics.size()+" Valid: "+output_statistics.stream().filter( si -> si.valid_split ).count());
               int total_possible_rxn_mappings    = output_statistics.stream().flatMapToInt(si -> si.possible_rxn_mappings.values().stream().mapToInt( mi -> mi.size() )).sum();
               int total_processed_labeled_splits = output_statistics.stream().mapToInt(si -> si.total_labeled_connector_splits_processed).sum();
               int total_enumerated_sss           = output_statistics.stream().mapToInt(si -> si.total_sss_performed).sum();
               System.out.println("Possible Rxn Maps: "+total_possible_rxn_mappings);
               System.out.println("Processed labeled splits: "+total_processed_labeled_splits);
               System.out.println("Synthon substructure searches: "+total_enumerated_sss);
           }
       }

       return expanded_hits;
    }


    public static interface CombinatorialHitReceiver {
        public void addCombinatorialHits(List<CombinatorialHit> hi, List<SplitPatternWithConnectorProximityPruningStatistics> stats);
    }

    /**
     * Same as function findExpandedHits_withConnProximityMatching, but search is not performed in stages but
     * each hit candidate hit is directly processed, and hits are streamed to the output during the computation.
     *
     * @param space
     * @param cdh
     * @param mol
     * @param num_splits
     * @param max_fragments
     * @param rxns_to_omit
     * @param max_hits maximum number of combinatorial hits that should be examined
     * @param threads
     * @param receiver
     * @return
     */
    public static void findExpandedHits_withConnProximityMatching_streaming(SynthonSpace space, CachedDescriptorProvider cdh ,
                                                                                    StereoMolecule mol, int num_splits, int max_fragments,
                                                                                    Set<String> rxns_to_omit ,
                                                                                    int max_hits, int threads,
                                                                                    CombinatorialHitReceiver receiver) {


        int nb = mol.getBonds();

        List<int[]> splits  = null;
        if(num_splits==0) {
            splits = new ArrayList<>();
            splits.add(new int[0]);
        }
        else {
            List<int[]> combi_list = CombinationGenerator.getAllOutOf(nb, num_splits);
            // !! returns null if b > a..
            if(combi_list==null) { return; }
            splits = combi_list.stream().filter(ci -> ci.length == num_splits).collect(Collectors.toList());
        }

        ExecutorService main_pool = Executors.newFixedThreadPool(threads);

        List<Future> tasks_sss = new ArrayList<>();

        for(int[] si : splits) {
            final int num_connectors = num_splits+1; // this we should also fix at the beginning..
            //findInitialHits_forSplitPattern(space,cdh,mol,si,num_connectors,max_fragments,max_hits,initial_hits,initial_hits_sorted);

            Runnable r_initialhits_and_expansion = new Runnable() {
                @Override
                public void run() {
                    List<CombinatorialHit> expanded_hits_i                    = Collections.synchronizedList(new ArrayList<>() ); // synchonrized not needed..
                    List<SplitPatternWithConnectorProximityPruningStatistics> output_statistics_i = Collections.synchronizedList(new ArrayList<>() );  // synchonrized not needed..
                    findExpandedHits_forSplitPattern_withConnProximityMatching(space,cdh,mol,si,num_connectors,max_fragments,max_hits,rxns_to_omit,expanded_hits_i,output_statistics_i);
                    if(!expanded_hits_i.isEmpty()) {
                        receiver.addCombinatorialHits(expanded_hits_i, output_statistics_i);
                    }
                }
            };
            tasks_sss.add(main_pool.submit(r_initialhits_and_expansion));
        }

        // wait for pool:
        if(logLevel_findCandidates>0) {
            System.out.println("[FCH] Substructure search: Start\n");
        }
        int cnt_fc = 0;
        for(Future f_i : tasks_sss) {
            try {
                if(Thread.interrupted()) {
                    f_i.cancel(true);
                }
                f_i.get();
                if(logLevel_findCandidates>0) {
                    //System.out.print(".");
                    if ( (++cnt_fc) % 60 == 0) {
                        //    System.out.print("\n");
                    }
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        if(logLevel_findCandidates>0) {
            System.out.println("\n[FCH] Substructure search: Done!");
        }

    }


    public static class SplitPatternWithConnectorProximityPruningStatistics {
        public final boolean valid_split;
        public final Map<String,List<Map<Integer,Pair<FragType,BitSet>>>> possible_rxn_mappings;
        public final int total_labeled_connector_splits_processed;
        public final int total_sss_performed;

        public SplitPatternWithConnectorProximityPruningStatistics(boolean valid_split, Map<String, List<Map<Integer, Pair<FragType, BitSet>>>> possible_rxn_mappings, int total_labeled_connector_splits_processed, int total_sss_performed) {
            this.valid_split = valid_split;
            this.possible_rxn_mappings = possible_rxn_mappings;
            this.total_labeled_connector_splits_processed = total_labeled_connector_splits_processed;
            this.total_sss_performed = total_sss_performed;
        }
    }

    /**
     * This function works differently than the original findInitialHits:
     * We first try to map the fragments to connector-proximal fps of specific reactions.
     *
     * If we find any reactions that can match, then we enumerate all unique-connector variants
     * of the split synthons, and check if we find a fingerprint hit for the biggest fragment
     * of the split.
     *
     * It this test is passed, then we start matching the unique connector synthons from the split
     * against the synthon reaction synthons, via substructure search.
     *
     * @param space
     * @param cdh
     * @param mol
     * @param split_pattern
     * @param num_connectors
     * @param max_fragments
     * @param max_hits
     * @param expanded_hits
     */
    public static void findExpandedHits_forSplitPattern_withConnProximityMatching(SynthonSpace space, CachedDescriptorProvider cdh , StereoMolecule mol , int[] split_pattern , int num_connectors, int max_fragments, int max_hits, Set<String> rxns_to_omit, List<CombinatorialHit> expanded_hits , List<SplitPatternWithConnectorProximityPruningStatistics> out_statistics ) {

        // check if we are done..
        if(expanded_hits.size() >= max_hits ){
            System.out.println("BREAK EARLY!! findInitialHits : reached max initial hits: "+expanded_hits.size());
            //break;
            return;
        }

        StereoMolecule mi = new StereoMolecule(mol);
        mi.ensureHelperArrays(StereoMolecule.cHelperCIP);

        SynthonShredder.SplitResult split_result = SynthonShredder.trySplit(mi,split_pattern,max_fragments);
        if(split_result==null) {
            if(out_statistics!=null) {
                out_statistics.add(new SplitPatternWithConnectorProximityPruningStatistics(false,new HashMap<>(),0,0));
            }
            return;
        }

        if(false) {
            System.out.println("Test Split: "+split_result.toString());
        }

        //List<SynthonShredder.SplitResult> split_variants = new ArrayList<>();
        //List<SynthonShredder.SplitResult> split_variants = SynthonShredder.getAllSplitVariations(split_result, num_connectors);

        // compute initial matchings
        //Map<String,List<Map<FragType,BitSet>>> possible_rxn_mappings = computeConnectorProximityPruningMatchedVariants(space,cdh,split_result);
        Map<String,List<Map<Integer,Pair<FragType,BitSet>>>> possible_rxn_mappings = computeConnectorProximityPruningMatchedVariants(space,cdh,split_result,rxns_to_omit);

        if(true) {
            if(possible_rxn_mappings.values().stream().mapToInt(vi -> vi.size()).sum() > 0 ) {
                System.out.println("Rxn Matchings found!");
                System.out.println("Split: "+split_result.toString());
                System.out.println("Rxns: "+possible_rxn_mappings.keySet().stream().filter( ki -> !possible_rxn_mappings.get(ki).isEmpty() ).collect(Collectors.joining(",")));
                System.out.println("");
            }
            else if(false) {
                System.out.println("No Rxn Matchings found..");
                System.out.println("Split: "+split_result.toString());
            }
        }

        // check if we find initial hits, if this is the case, we compute all splits with unique connectors.
        // !!! NOTE: we must search through different unique connector configurations, depending on the
        // !!!       number of connectors in the reaction that we consider
        //Map<Integer,List<StereoMolecule[]>> unique_splits_for_reactions_with_n_connectors = new HashMap<>(); // Maps from number of connectors to the list of molecules to check
        Map<String,List<StereoMolecule[]>> unique_splits_for_reaction = new HashMap<>();

        // remove all possible_rxn_mappings without values:
        List<String> rxns_to_remove = new ArrayList<>();
        for(String ri : possible_rxn_mappings.keySet()) {
            if(possible_rxn_mappings.get(ri).isEmpty()) {rxns_to_remove.add(ri);}
        }
        for(String ri : rxns_to_remove){ possible_rxn_mappings.remove(ri); }

        if( possible_rxn_mappings.values().size() > 0 ) {
            // determine for each reaction the number of connectors.
            List<String> possible_rxns = new ArrayList<>(possible_rxn_mappings.keySet());
            Map<List<Integer>,List<StereoMolecule[]>> cachedLabeledFragsForConnectorSet = new HashMap<>();
            //int[] numbers_of_connectors = possible_rxns.stream().mapToInt(
            //                ri -> space.getRxn(ri).num_connectors ).distinct().toArray();
            for(int zi=0;zi<possible_rxns.size();zi++) {
                //int nc = numbers_of_connectors[zi];
                String ri = possible_rxns.get(zi);
                int nc = space.getRxn(ri).num_connectors;
                List<Integer> connis_to_use = new ArrayList<>(); for(int zn=0;zn<nc;zn++){ connis_to_use.add(92+zn); }
                if(!cachedLabeledFragsForConnectorSet.containsKey(connis_to_use)) {
                    List<StereoMolecule[]> unique_connector_splits = split_result.getAllSplitsWithUniqueConnectors(connis_to_use);
                    cachedLabeledFragsForConnectorSet.put(connis_to_use,unique_connector_splits);
                }
                List<StereoMolecule[]> unique_connector_splits = new ArrayList<>( cachedLabeledFragsForConnectorSet.get(connis_to_use) );

                // here we perform the "largest fragment existance" pruning!
                // NOTE!! this version did not always work, because the labeled fragment can contain supersets of the
                // matched labeled fragment. I.e. we would have to check not only frag_0_connectorset, but in fact
                // the connectorset of the synthon of the synthon reaction, which in some cases contains additional
                // connectors. We simplified this and just run the largest labeled fragment check against all
                // synthons in the reaction, then this should automatically work.
                if(false) {
                    List<StereoMolecule[]> to_prune = new ArrayList<>();
                    for(StereoMolecule[] split_i : unique_connector_splits) {
                        //unique_connector_splits;
                        BitSet frag_0_connectorset = determineConnectorsInMolecule(split_i[0]);
                        BitSet frag_0_fp = space.getFP(split_i[0]);
                        boolean test = false; // do we find initial hit?
                        BitSetTree bst_frag0 = null;
                        if (space.bsts_by_connectors.containsKey(frag_0_connectorset)) {
                            bst_frag0 = space.bsts_by_connectors.get(frag_0_connectorset);
                            BitSetTree.Node[] result = new BitSetTree.Node[1];
                            test = bst_frag0.testSubset(frag_0_fp, result);
                        }
                        if(!test) {
                            // remove this one
                            to_prune.add(split_i);
                        }
                    }
                    unique_connector_splits.removeAll(to_prune);
                    System.out.println("LargestFrag Heuristic removed "+to_prune.size()+" / "+(to_prune.size()+unique_connector_splits.size())+" UniqueConnector Splits");
                }
                // new (simpler) largest fragment heuristic, just check all labeled synthons of the matched synthon rxn.
                if(true) {
                    List<StereoMolecule[]> to_prune = new ArrayList<>();
                    for(StereoMolecule[] split_i : unique_connector_splits) {
                        //unique_connector_splits;
                        BitSet frag_0_fp = space.getFP(split_i[0]);
                        boolean test = false; // do we find initial hit?
                        BitSetTree bst_frag0 = space.bsts_labeled_by_rxn.get(possible_rxns.get(zi));
                        BitSetTree.Node[] result = new BitSetTree.Node[1];
                        test = bst_frag0.testSubset(frag_0_fp, result);
                        if(!test) {
                            // remove this one
                            to_prune.add(split_i);
                        }
                    }
                    unique_connector_splits.removeAll(to_prune);
                    System.out.println("LargestFrag Heuristic removed "+to_prune.size()+" / "+(to_prune.size()+unique_connector_splits.size())+" UniqueConnector Splits");
                }

                //unique_splits_for_reactions_with_n_connectors.put( nc , unique_connector_splits );
                unique_splits_for_reaction.put(ri,unique_connector_splits);
            }
              //unique_splits = split_result.getAllSplitsWithUniqueConnectors();
        }

        int total_labeled_connector_splits_processed = 0;
        int total_enumerated_sss_performed = 0;

        for(String rxn : possible_rxn_mappings.keySet()) {

            if(true) { // debug output:
                System.out.println("All Matchings:");
                for(Map<Integer,Pair<FragType,BitSet>> mpi : possible_rxn_mappings.get(rxn)) {
                    for(Integer xfi : mpi.keySet().stream().sorted( Integer::compareTo ).collect(Collectors.toList())) {
                        System.out.println(xfi+" -> "+mpi.get(xfi).getLeft()+" : "+mpi.get(xfi).getRight());
                    }
                }
            }

            int ri_num_connis = space.getRxn(rxn).num_connectors;
            for(Map<Integer,Pair<FragType,BitSet>> mpi : possible_rxn_mappings.get(rxn)) {
                // Now iterate over the unique connector splits. We require these to match the actual
                // substructures:
                //for (StereoMolecule[] frags : unique_splits_for_reactions_with_n_connectors.get(ri_num_connis)) {
                for (StereoMolecule[] frags : unique_splits_for_reaction.get(rxn)) {
                    total_labeled_connector_splits_processed++; // only for statistics
                    if(true) {
                        System.out.println("Split synthons to search: "+ Arrays.stream(frags).map( mii -> HyperspaceUtils.idcodeToSmiles(mii.getIDCode()) ).collect(Collectors.joining(".")) );
                        System.out.println("Split synthons to search (IDCODE)): "+Arrays.stream(frags).map( mii -> mii.getIDCode() ).collect(Collectors.joining(" : ")) );
                    }

                    Map<FragType, List<FragId>> matching_frag_ids = new HashMap<>();

                    // now we iterate over the fragments and see if we find matching ones   // TODO: think about sorting, biggest or smallest first??
                    for( Integer split_idx : mpi.keySet().stream().sorted().collect(Collectors.toList())) {
                        FragType fti = mpi.get(split_idx).getLeft();
                        BitSet bsi = mpi.get(split_idx).getRight();

                        // this is the fragment to find:
                        StereoMolecule ffi = frags[split_idx];

                        // !!! NOW: we have to enumerate all supersets of the bsi connector fingerprint in the
                        // corresponding bitset tree !!!
                        List<BitSet> supersets_bsi = new ArrayList<>();
                        space.connector_fps_sorted_by_fragtype.get(fti).root.collectSuperSets(bsi,supersets_bsi);

                        //List<String> hits_i = new ArrayList<>();
                        Map<String,List<String>> resolved_ids = new HashMap<>();

                        if (space.substructure_searchers_sorted_by_connector_fp.containsKey(fti)) {
                            int info_num_structures_searches = 0;
                            boolean found_frag_hit = false;
                            for(BitSet ss_bsi : supersets_bsi) {
                                if (space.substructure_searchers_sorted_by_connector_fp.get(fti).containsKey(ss_bsi)) {
                                    ffi.setFragment(true);
                                    info_num_structures_searches+=space.substructure_searchers_sorted_by_connector_fp.get(fti).get(ss_bsi).mIDs.keySet().size(); // only for info/debug
                                    List<String> hits_i = space.substructure_searchers_sorted_by_connector_fp.get(fti).get(ss_bsi).findSubstructure(cdh, ffi, 200); // TODO: parametrize max_hits parameter
                                    //resolved_ids = space.substructure_searchers_sorted_by_connector_fp.get(fti).get(ss_bsi).resolveIdcodeToIDs(hits);
                                    Map<String,List<String>> resolved_ids_i = space.substructure_searchers_sorted_by_connector_fp.get(fti).get(ss_bsi).resolveIdcodeToIDs(hits_i);
                                    resolved_ids.putAll(resolved_ids_i);
                                    if(!resolved_ids_i.isEmpty()) { found_frag_hit = true; }
                                }
                            }
                            total_enumerated_sss_performed+=info_num_structures_searches; // only for statistics
                            if(true) { // full debug output..
                                String info_a = rxn+"_idc="+ffi.getIDCode()+"_fti="+fti.toString()+"_bsi="+bsi+" tot_ss_performed="+info_num_structures_searches;
                                String str_all_info = info_a;
                                if (found_frag_hit) {
                                    System.out.println("FRAG_TEST: HIT!! -> "+str_all_info);
                                }
                                else {
                                    System.out.println("FRAG_TEST: NOHIT -> "+str_all_info);
                                }
                            }
                        }

                        if(resolved_ids.isEmpty()) {
                            // we can just break, as we do not manage to find matching frags for
                            // one of the split synthons in this split.
                            break;
                        }

                        // resolve all hits to FragId objects:
                        List<String> all_resolved_ids = resolved_ids.values().stream().flatMap( vi -> vi.stream() ).collect(Collectors.toList());
                        List<FragId> matching_frag_ids_i = new ArrayList<>();
                        all_resolved_ids.stream().forEach( vi -> matching_frag_ids_i.add( space.frags_by_id.get(vi) ) );

                        // add to results:
                        if(!matching_frag_ids.containsKey(fti)) {matching_frag_ids.put(fti,new ArrayList<>());}
                        matching_frag_ids.get(fti).addAll(matching_frag_ids_i);
                    }

                    // now add a combinatorial hit in case that all frag types were found:
                    if (matching_frag_ids.keySet().stream().filter(ki -> matching_frag_ids.get(ki).size() > 0).count() == mpi.keySet().size()) {
                        System.out.println("FULL_HIT! -> "+rxn+" -> "+matching_frag_ids);
                        expanded_hits.add(new CombinatorialHit(rxn, matching_frag_ids, split_result, mpi));
                    }

                }
            }

        }

        if(out_statistics!=null) {
            out_statistics.add(new SplitPatternWithConnectorProximityPruningStatistics(true,possible_rxn_mappings,total_labeled_connector_splits_processed,total_enumerated_sss_performed));
        }

        if(false) {
            /**
             * The following is working, but it is not efficient, and it will generate a lot
             * of false positives. This was replaced by the code above..
             */
            for (String rxn : possible_rxn_mappings.keySet()) {
                for (Map<Integer, Pair<FragType, BitSet>> mpi : possible_rxn_mappings.get(rxn)) {
                    Map<FragType, List<FragId>> matching_frag_ids = new HashMap<>();
                    // now we iterate over the fragments and see if we find matching ones   // TODO: think about sorting, biggest or smallest first??
                    for (Integer split_idx : mpi.keySet().stream().sorted().collect(Collectors.toList())) {
                        FragType fti = mpi.get(split_idx).getLeft();
                        BitSet bsi = mpi.get(split_idx).getRight();
                        // NOTE! Now we have to use all connector configurations containing the
                        // connector types from the FragType!!
                        List<Integer> connis = space.getConnectorSetForFragType(fti).stream().map(pi -> pi + 92).boxed().collect(Collectors.toList());
                        List<List<Integer>> all_conni_perms = SynthonSpace.all_permutations(connis);

                        StereoMolecule template_to_find = split_result.fragments[split_idx];
                        template_to_find.ensureHelperArrays(Molecule.cHelperNeighbours);
                        int all_connector_positions[] = new int[connis.size()];
                        Arrays.fill(all_connector_positions, -1);
                        int idx_acp = 0;
                        for (int zp = 0; zp < template_to_find.getAtoms(); zp++) {
                            if (template_to_find.getAtomicNo(zp) == 92) {
                                all_connector_positions[idx_acp++] = zp;
                            }
                        }

                        for (List<Integer> conni_perm : all_conni_perms) {
                            StereoMolecule ffi = new StereoMolecule(template_to_find);
                            ffi.ensureHelperArrays(Molecule.cHelperNeighbours);
                            // put connector atom permutation:
                            for (int zp = 0; zp < all_connector_positions.length; zp++) {
                                ffi.setAtomicNo(all_connector_positions[zp], conni_perm.get(zp));
                            }
                            ffi.ensureHelperArrays(Molecule.cHelperCIP);
                            List<String> hits = new ArrayList<>();
                            Map<String, List<String>> resolved_ids = new HashMap<>();
                            if (space.substructure_searchers_sorted_by_connector_fp.containsKey(fti)) {
                                if (space.substructure_searchers_sorted_by_connector_fp.get(fti).containsKey(bsi)) {
                                    hits = space.substructure_searchers_sorted_by_connector_fp.get(fti).get(bsi).findSubstructure(cdh, ffi, 200); // TODO: parametrize max_hits parameter
                                    resolved_ids = space.substructure_searchers_sorted_by_connector_fp.get(fti).get(bsi).resolveIdcodeToIDs(hits);
                                    ;
                                }
                            }

                            // resolve all hits to FragId objects:
                            List<String> all_resolved_ids = resolved_ids.values().stream().flatMap(vi -> vi.stream()).collect(Collectors.toList());
                            List<FragId> matching_frag_ids_i = new ArrayList<>();
                            all_resolved_ids.stream().forEach(vi -> matching_frag_ids_i.add(space.frags_by_id.get(vi)));

                            // add to results:
                            if (!matching_frag_ids.containsKey(fti)) {
                                matching_frag_ids.put(fti, new ArrayList<>());
                            }
                            matching_frag_ids.get(fti).addAll(matching_frag_ids_i);
                        }
                    }

                    // now add a combinatorial hit in case that all frag types were found:
                    if (matching_frag_ids.keySet().stream().filter(ki -> matching_frag_ids.get(ki).size() > 0).count() == mpi.keySet().size()) {
                        expanded_hits.add(new CombinatorialHit(rxn, matching_frag_ids, split_result, mpi));
                    }

                }
            }
        }

    }


    public static void findInitialHits_forSplitPattern(SynthonSpace space, CachedDescriptorProvider cdh , StereoMolecule mol , int[] split_pattern , int num_connectors, int max_fragments, int max_hits, List<InitialHit> initial_hits , Map<FragType,List<InitialHit>> initial_hits_sorted) {

        // check if we are done..
        if(initial_hits_sorted.values().stream().mapToInt( vi -> vi.size() ).sum() >= max_hits ){
            System.out.println("BREAK EARLY!! findInitialHits : reached max initial hits: "+initial_hits_sorted.values().stream().mapToInt( vi -> vi.size() ).sum());
            //break;
            return;
        }

        StereoMolecule mi = new StereoMolecule(mol);
        mi.ensureHelperArrays(StereoMolecule.cHelperCIP);

        SynthonShredder.SplitResult split_result = SynthonShredder.trySplit(mi,split_pattern,max_fragments);
        if(split_result==null) {return;}

        //List<SynthonShredder.SplitResult> split_variants = new ArrayList<>();
        List<SynthonShredder.SplitResult> split_variants = SynthonShredder.getAllSplitVariations(split_result, num_connectors);

        // here we do the connector-proximity pruning: NO: actually, this does not really work nicely, as
        // the whole rest of the function also will work slightly differently.
//        boolean proximity_connector_pruning = true;
//        if(proximity_connector_pruning) {
//            Map<String,List<Map<FragType,BitSet>>> possible_rxn_mappings = computeConnectorProximityPruningMatchedVariants(space,cdh,split_result);
//        }
//        else {
//            split_variants = SynthonShredder.getAllSplitVariations(split_result, num_connectors);
//        }


        for(SynthonShredder.SplitResult split : split_variants) {

            StereoMolecule frags[] = split.fragments;
            BitSet bs = cdh.getFP_cached(frags[0]);

            //Set<Integer> frag_0_connectorset = computeConnectorSet(frags[0]);
            BitSet frag_0_connectorset = determineConnectorsInMolecule(frags[0]);

            BitSetTree.Node result[] = new BitSetTree.Node[1];
            //boolean test = bst.testSubset(bs, result);

            boolean test = false; // do we find initial hit?
            BitSetTree bst_frag0 = null;
            if(space.bsts_by_connectors.containsKey(frag_0_connectorset)) {
                bst_frag0 = space.bsts_by_connectors.get(frag_0_connectorset);
                test = bst_frag0.testSubset(bs, result);
            }
            int hits_cnt = 0;

            if (test) {

                if (space.logLevel_findCandidates > 0) {
                    System.out.println("cut set with candidate hit [fp]");
                    for (int zi = 0; zi < frags.length; zi++) {
                        System.out.println(HyperspaceUtils.idcodeToSmiles(frags[zi].getIDCode()));
                    }
                }

                // count the connectors (we need this then when we check for actual library fragments)
                //int frag_cnt_u = 0; int frag_cnt_np = 0; int frag_cnt_pu = 0; int frag_cnt_am = 0;
//                for(int zi=0;zi<frags[0].getAtoms();zi++) {
//                    frag_cnt_u  += (frags[0].getAtomicNo(zi)==92)?1:0;
//                    frag_cnt_np += (frags[0].getAtomicNo(zi)==93)?1:0;
//                    frag_cnt_pu += (frags[0].getAtomicNo(zi)==94)?1:0;
//                    frag_cnt_am += (frags[0].getAtomicNo(zi)==95)?1:0;
//                }
                int frag_connector_counts[] = new int[MAX_CONNECTORS];
                for(int zi=0;zi<frags[0].getAtoms();zi++) {
                    int an = frags[0].getAtomicNo(zi);
                    if( an>=92 && an < 92+MAX_CONNECTORS ) { frag_connector_counts[an-92]++;}
                }


                // store the cut pattern idcodes:
                //List<String> cut_pattern_idcodes = new ArrayList<>();
                List<String> cut_pattern_idcodes = Arrays.stream(frags).map(fi -> fi.getIDCode()).collect(Collectors.toList());
//                    List<Integer> cut_set = Arrays.stream(si).mapToObj(xi -> new Integer(xi)).collect(Collectors.toList());
                String cut_pattern_hash = split.cutset_hash;
                space.cut_fragments.put( cut_pattern_hash , cut_pattern_idcodes);

                //System.out.println("Initial Hit!");
                List<BitSet> hits_i = new ArrayList<>();
                bst_frag0.root.collectSuperSets(bs, hits_i);

                // sort frag hits:
                Map<FragType, List<FragId>> sorted_frags = new HashMap<>();
                Set<String> all_hit_rxns = new HashSet<>();
                for (BitSet bhit : hits_i) {
                    //Set<FragId> frag_hits = fragments.get(bhit);
                    Set<FragId> frag_hits = space.getAllFragments_byConnectors(frag_0_connectorset).get(bhit);

                    for (FragId fid : frag_hits) {
                        // check that reaction is compatible with number of fragments, i.e. the reaction must have at least
                        // as many fragments as the current split
                        if (space.ffps_sorted_by_rxn_and_frag_BT.get(fid.rxn_id).keySet().size() < frags.length) {
                            // no hit..
                            continue;
                        }

                        // check the connector set of the rxn and compare against this split:
                        if(space.rxns_by_connector_config.get(split.connector_config)==null) {
                            //System.out.println("no rxns for this cutset config.. skip");
                            continue;
                        }
                        if( ! space.rxns_by_connector_config.get(split.connector_config).contains(fid.rxn_id) ) {
                            //System.out.println("rxn not compatible with cutset config.. skip");
                            continue;
                        }


                        IDCodeParser icp = new IDCodeParser();
                        StereoMolecule fmi = new StereoMolecule();
                        icp.parse(fmi, fid.idcode);

                        // This check now is redundant. We leave it in for the moment, as a sanity check
                        if(true) {
                            // !! The primary fragment MUST have the same number of connectors !!
                            int lf_connector_counts[] = new int[MAX_CONNECTORS];
                            for(int zi=0;zi<fmi.getAtoms();zi++) {
                                int an = fmi.getAtomicNo(zi);
                                if( an>=92 && an < 92+MAX_CONNECTORS ) { lf_connector_counts[an-92]++;}
                            }

                            boolean difference = false;
                            for(int zi=0;zi<MAX_CONNECTORS;zi++){ difference |= (frag_connector_counts[zi]!=lf_connector_counts[zi]); }
                            if (difference) {
                                //throw new Error("We cannot end up here. Catastrophic error in findInitialHits(..)! Wrong connectors..");
                                System.out.println("[ERROR] We cannot end up here. Catastrophic error in findInitialHits(..)! Wrong connectors..");
                                System.out.println("[ERROR] This might be due to wrong connectors in synthon set.. synthon set: "+fid.rxn_id+" : "+fid.frag+" : "+ HyperspaceUtils.idcodeToSmiles(fid.idcode));
                                //continue;
                            }
                        }

                        // Perform SSS to check if we really have a hit.
                        SSSearcher sss = new SSSearcher();
                        frags[0].setFragment(true);
                        sss.setMol(frags[0],fmi);
                        boolean primary_fragment_match = sss.isFragmentInMolecule();
                        if(!primary_fragment_match) {
                            if (space.logLevel_findCandidates > 1) {
                                System.out.println("Primary fragment after SSS: no match");
                            }
                            continue;
                        }

                        FragType ft = new FragType(fid.rxn_id, fid.frag);
                        if (!sorted_frags.containsKey(ft)) {
                            sorted_frags.put(ft, new ArrayList<>());
                        }
                        sorted_frags.get(ft).add(fid);
                        all_hit_rxns.add(ft.rxn_id);

                        if (!initial_hits_sorted.containsKey(ft)) {
                            initial_hits_sorted.put(ft, new ArrayList<>());
                        }

                        // identify the cutset-fragments via the fragment bitsets:
                        Map<Integer, BitSet> map_cutset_frag_to_bitset = new HashMap<>();
                        map_cutset_frag_to_bitset.put(0, bs); // add the root fingerprint

                        // remaining fps:
                        //List<BitSet> remaining_frags = new ArrayList<>();
                        Map<Set<Integer>,List<BitSet>> remaining_frags_fp = new HashMap<>();
                        Map<Set<Integer>,List<StereoMolecule>> remaining_frags    = new HashMap<>();

                        //Map<Set<Integer>,StereoMolecule> remaining_frags = new HashMap<>();
                        List<Set<Integer>> unmatched_conn_sets = new ArrayList<>();
                        for (int zi = 1; zi < frags.length; zi++) {
                            frags[zi].setFragment(true);
                            BitSet bs_fp_i = cdh.getFP_cached(frags[zi]);

                            Set<Integer> conn_set = computeConnectorSet(frags[zi]);
                            //remaining_frags.add(bs_fp_i);
                            unmatched_conn_sets.add(conn_set);

                            if(!remaining_frags_fp.containsKey(conn_set)) {  remaining_frags_fp.put( conn_set , new ArrayList<>());}
                            if(!remaining_frags.containsKey(conn_set)) {     remaining_frags.put(    conn_set , new ArrayList<>());}

                            remaining_frags_fp.get(conn_set).add(bs_fp_i);
                            remaining_frags.get(conn_set).add(frags[zi]);

                            map_cutset_frag_to_bitset.put(zi, bs_fp_i);
                        }

                        Map<Set<Integer>,List<Integer>> possible_frags = new HashMap<>();

                        //List<List<Integer>> rxn_connector_conf = decodeConnectorConfig(split.connector_config); // we KNOW that the split has the correct conn. config for this rxn! (because we checked it before)
                        //List<BitSet> rxn_connector_conf = decodeConnectorConfig(split.connector_config); // we KNOW that the split has the correct conn. config for this rxn! (because we checked it before)
                        //List<BitSet> rxn_connector_conf = split.connector_config


                        // now compute all fragment assignments:
                        String rxn = fid.rxn_id;
                        SynthonShredder.SplitResult svi = split;

                        // check if the assignment is ambiguous (it should not be..)
                        // 1. Init the map with possibilities, 2. remove the primary hit fragment, 3. check..
                        //Map<Set<Integer>,List<FragType>> possible_mappings = new HashMap<>( space.fragment_types_by_rxn_and_connector_config.get(rxn) );
                        //possible_mappings.get(  )
                        boolean ambiguous = unmatched_conn_sets.stream().distinct().count() < unmatched_conn_sets.size();
                        if(ambiguous) {
                            System.out.println("ERROR: abmiguous mapping -> not supported, please fix your combinatorial library! Problem in rxn: "+rxn);
                        }
                        else {
                            Map<FragType,BitSet> map_extensions = new HashMap<>();
                            Map<FragType,StereoMolecule> map_extension_templates = new HashMap<>();

                            for (int zi = 1; zi < svi.fragments.length; zi++) {
                                Set<Integer> connectors_fi = computeConnectorSet(svi.fragments[zi]);
                                BitSet fp_fi = cdh.getFP_cached(svi.fragments[zi]);

                                // find fragment type:
                                FragType fti = space.fragment_types_by_rxn_and_connector_config.get(rxn).get(connectors_fi).get(0);
                                map_extensions.put(fti,fp_fi);
                                map_extension_templates.put(fti,svi.fragments[zi]);
                            }

                            InitialHit ihit = new InitialHit(bhit, ft.rxn_id, ft.frag, fid, map_extensions, map_extension_templates, cut_pattern_hash, map_cutset_frag_to_bitset);
                            initial_hits.add(ihit);
                            initial_hits_sorted.get(ft).add(ihit);
                            hits_cnt++;
                        }


//                        if( space.fragment_types_by_rxn_and_connector_config.get(rxn).values().stream().mapToInt( li -> li.size() ).allMatch( siz -> siz==1 ) ) {
//                            // âll unique connector configs, that's easy..
//                            Map<FragType,BitSet> map_extensions = new HashMap<>();
//                            Map<FragType,StereoMolecule> map_extension_templates = new HashMap<>();
//
//                            for (int zi = 1; zi < svi.fragments.length; zi++) {
//                                Set<Integer> connectors_fi = computeConnectorSet(svi.fragments[zi]);
//                                BitSet fp_fi = cdh.getFP_cached(svi.fragments[zi]);
//
//                                // find fragment type:
//                                FragType fti = space.fragment_types_by_rxn_and_connector_config.get(rxn).get(connectors_fi).get(0);
//                                map_extensions.put(fti,fp_fi);
//                                map_extension_templates.put(fti,svi.fragments[zi]);
//                            }
//
//                            InitialHit ihit = new InitialHit(bhit, ft.rxn_id, ft.frag, fid, map_extensions, map_extension_templates, cut_pattern_hash, map_cutset_frag_to_bitset);
//                            initial_hits.add(ihit);
//                            initial_hits_sorted.get(ft).add(ihit);
//                            hits_cnt++;
//
//                        }
//                        else {
//                            // shared connector configs. Here it is a bit more complicated..
//                            List<Set<Integer>> template_list = new ArrayList<>();
//
//                            for(int zx = 1; zx<svi.fragments.length ; zx++ ) {
//                                Set<Integer> cs_ext_a = computeConnectorSet(svi.fragments[zx]);
//                                template_list.add(cs_ext_a);
//                            }
//
//                            // frag types (important, remove the primary hit frag from the possibles!)
//                            Map<Set<Integer>,List<FragType>> frag_options = space.fragment_types_by_rxn_and_connector_config.get(rxn);
//                            frag_options.values().stream().forEach( li -> li.remove(ft) );
//
//                            List<List<FragType>> all_assignments = SynthonSimilaritySpace.computeAllFragTypeOrderings( frag_options , template_list);
//
//                            for(List<FragType> assi : all_assignments) {
//                                Map<FragType,BitSet> map_extensions = new HashMap<>();
//                                Map<FragType,StereoMolecule> map_extension_templates = new HashMap<>();
//
//                                for(int zx = 0; zx<assi.size() ; zx++ ) {
//                                    FragType fti = assi.get(zx);
//                                    StereoMolecule frag_i = svi.fragments[zx+1];
//                                    map_extensions.put(fti, cdh.getFP_cached(frag_i));
//                                    map_extension_templates.put(fti,frag_i);
//                                }
//                                InitialHit ihit = new InitialHit(bhit, ft.rxn_id, ft.frag, fid, map_extensions, map_extension_templates, cut_pattern_hash, map_cutset_frag_to_bitset);
//                                initial_hits.add(ihit);
//                                initial_hits_sorted.get(ft).add(ihit);
//                                hits_cnt++;
//                            }
//                        }
                    }
                }
                if (space.logLevel_findCandidates > 0) {
                    System.out.println("Found " + hits_cnt + " initial hits in " + sorted_frags.size() +
                            " fragment types and in " + all_hit_rxns.size() + " synthon sets.");
                }
            }
        }
    }

    /**
     * Computes for all synthon rxns the possible matchings. For all computed
     * possible assignments, the integer describes the idx of the split fragment,
     * and the BitSet describes the required connector fp bitset (note: the required bits!).
     *
     * NOTE!! for substructure search, the requirement for matching is superset, therefore,
     * when we search for the actual fragments, we MUST FIRST collect
     * all the supersets of the stored connector fp bitset!!!
     *
     *
     * @param cdp
     * @param split_result
     * @return
     */
    public static Map<String,List<Map<Integer, Pair<FragType,BitSet>>>> computeConnectorProximityPruningMatchedVariants(SynthonSpace space, CachedDescriptorProvider cdp, SynthonShredder.SplitResult split_result, Set<String> rxns_to_omit) {


        // NOT_TODO: filter rxns by connector counts.. --> ACTUALLY, we cannot really filter based on
        // this, as we then miss partial hits. We have to select all reactions, that contain
        // a subset of the required connector counts.
        //List<String> rxns = space.rxns_by_connector_config.get( split_result.connector_config );
        //if( rxns==null || rxns.isEmpty()) { return new HashMap<>(); }
//        List<String> rxns = space.rxns_by_connector_counts.get( split_result.connector_counts );
//        if(rxns==null || rxns.isEmpty()) {
//            return new HashMap<>();
//        }
        // Correct filtering:
        List<String> rxns = new ArrayList<>();
        for(List<Integer> rxn_cc : space.rxns_by_connector_counts.keySet() ) {
            ArrayList<Integer> rxci = new ArrayList<>(rxn_cc);
            // try to assign smallest first
            boolean is_ok = true;
            for(int zi=0;zi<split_result.connector_counts.size();zi++) {
                int sci = split_result.connector_counts.get( split_result.connector_counts.size()-1-zi );
                Optional<Integer> oi = rxci.stream().filter( ci -> ci >= sci ).sorted().findFirst();
                is_ok &= oi.isPresent();
                if(!is_ok){break;}
                rxci.remove(new Integer(sci));
            }
            if(is_ok) {
                rxns.addAll(space.rxns_by_connector_counts.get(rxn_cc));
            }
        }

        // remove reactions that we should omit:
        rxns.removeAll(rxns_to_omit);

        // try the matching: iterate over split fragments and try to match them:
        int frag_idx = 0;

        Map<Integer,BitSet> stored_conn_fps = new HashMap<>(); // from frag-idx to stored conn-fp

        Map<Integer,BitSet> stored_conn_fps_MATCHED = new HashMap<>(); // from frag-idx to the conn-fp

        boolean matching_possible = true;

        Map<String,List<Map<FragType,Integer>>> initial_partial_matchings = new HashMap<>(); // per rxn id.., map from fragtype to idx in split result
        rxns.stream().forEach( si -> initial_partial_matchings.put(si,new ArrayList<>()) );
        rxns.stream().forEach( si -> initial_partial_matchings.get(si).add(new HashMap<>())); // empty mapping to initialize

        Map<String,List<Map<FragType,Integer>>> current_partial_matchings = initial_partial_matchings;
        while( matching_possible && (frag_idx<split_result.fragments.length) ) {
            BitSet conn_fpi = cdp.getFP_cached( createConnectorProximalFragment( split_result.fragments[frag_idx] , CONNECTOR_REGION_SIZE ) );
            stored_conn_fps.put(frag_idx,conn_fpi);

            matching_possible = false; // will be set to true, if we add extended partial mapping

            Map<String,List<Map<FragType,Integer>>> new_partial_matchings = new HashMap<>();
            for(String rxi : current_partial_matchings.keySet() ) {
                // try to extend the existing mappings:
                List<Map<FragType,Integer>> extended_maps = new ArrayList<>();
                for(Map<FragType,Integer> pmap_i : current_partial_matchings.get(rxi)){
                    // compute available frags
                    Set<FragType> available = new HashSet<>(space.fragment_type_map.get(rxi).values());
                    available.removeAll(pmap_i.keySet());

                    // loop over available, and check for each if the connector region is ok:
                    for(FragType fti : available) {
                        BitSetTree.Node first_superset[] = new BitSetTree.Node[1];
                        boolean test = space.connector_fps_sorted_by_fragtype.get(fti).testSubset(conn_fpi,first_superset);
                        if(test) {
                            Map<FragType,Integer> ext_map = new HashMap<>(pmap_i);
                            ext_map.put(fti,frag_idx);
                            extended_maps.add(ext_map);
                            matching_possible |= true;
                        }
                    }
                }
                new_partial_matchings.put(rxi,extended_maps);
            }
            current_partial_matchings = new_partial_matchings;
            frag_idx++;
        }

        if(matching_possible) {
            // create the result: i.e. replace frag indexes with conn-fps.
            //Map<String,List<Map<FragType,BitSet>>> result = new HashMap<>();
            Map<String,List<Map<Integer,Pair<FragType,BitSet>>>> result = new HashMap<>();

            for(String rxi : current_partial_matchings.keySet()) {
                result.put(rxi,new ArrayList<>());
                for( Map<FragType,Integer> map_i : current_partial_matchings.get(rxi)) {
                    //Map<FragType,BitSet> map_2_i = new HashMap<>();
                    Map<Integer,Pair<FragType,BitSet>> map_2_i = new HashMap<>();
                    //map_i.keySet().stream().forEach( ki -> map_2_i.put(ki, stored_conn_fps.get(map_i.get(ki)) ) );
                    map_i.keySet().stream().forEach( ki -> map_2_i.put( map_i.get(ki) , Pair.of( ki, stored_conn_fps.get(map_i.get(ki)) ) ) );
                    result.get(rxi).add(map_2_i);
                }
            }
            return result;
        }
        else {
            return new HashMap<>();
        }
    }

    public static Set<Integer> computeConnectorSet(StereoMolecule m) {
        m.ensureHelperArrays(Molecule.cHelperNeighbours);
        Set<Integer> ch = new HashSet<>();
        for(int zi=0;zi<m.getAtoms();zi++) {
            int an = m.getAtomicNo(zi);
            if(an>=92 && an < 92+MAX_CONNECTORS) {
                ch.add(an);
            }
        }
        return ch;
    }
    public static BitSet computeConnectorBitSet(StereoMolecule m) {
        BitSet bsi = new BitSet(16);
        m.ensureHelperArrays(Molecule.cHelperNeighbours);
        for(int zi=0;zi<m.getAtoms();zi++) {
            int an = m.getAtomicNo(zi);
            if(an>=92 && an < 92+MAX_CONNECTORS) {
                bsi.set(an-92);
            }
        }
        return bsi;
    }

    /**
     * Stores for each connector type the atom position.
     *
     * @param m
     * @return
     */
    public static Map<Integer,Integer> computeConnectorPositions(StereoMolecule m) {
        m.ensureHelperArrays(Molecule.cHelperNeighbours);
        Map<Integer,Integer> cm = new HashMap<>();
        for(int zi=0;zi<m.getAtoms();zi++) {
            int an = m.getAtomicNo(zi);
            if(an>=92 && an < 92+MAX_CONNECTORS) {
                cm.put(an,zi);
            }
        }
        return cm;
    }

    public static List<SynthonSpace.ExpandedHit> screen_building_blocks(SynthonSpace space , CachedDescriptorProvider cdh , StereoMolecule m, int threads) {

        List<SynthonSpace.ExpandedHit> hits = Collections.synchronizedList(new ArrayList<>() );

        //BitSet fp = this.getFP( getDescriptorHandler_PerThread(), m );
        BitSet fp = cdh.getFP_cached( m );

        //IDCodeParser icp = new IDCodeParser();

//        SSSearcher sss = new SSSearcher();
//        m.setFragment(true);
//        sss.setFragment(m);

        m.setFragment(true);
        SSSearcherProvider sss_provider = new SSSearcherProvider(m);

//        for(String ri : this.ffps_sorted_by_rxn_and_frag_BT.keySet()) {
//            for(int fi : this.ffps_sorted_by_rxn_and_frag_BT.get(ri).keySet()) {
        ExecutorService main_pool = Executors.newFixedThreadPool(threads);
        List<Future> tasks = new ArrayList<>();

        for(String ri : space.ffps_sorted_by_rxn_and_frag_BT.keySet()) {
            Runnable runnable_ri = new Runnable() {
                @Override
                public void run() {
                    SSSearcher sss = sss_provider.getSSS();
                    hits.addAll(screen_building_blocks_for_rxn(space,cdh,m,fp,sss,ri));
                }
            };
            Future fi = main_pool.submit(runnable_ri);
            tasks.add(fi);
        }
        // wait for pool:
        if(logLevel_findCandidates>0) {
            System.out.println("[SCB] Start screening building blocks:+\n");
        }

        int cnt_fc = 0;
        for(Future f_i : tasks) {
            try {
                f_i.get();
                if(logLevel_findCandidates>0) {
                    System.out.print(".");
                    if ( (++cnt_fc) % 60 == 0) {
                        System.out.print("\n");
                    }
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        if(logLevel_findCandidates>0) {
            System.out.println("\n[SCB] Start screening building blocks: Done!");
        }

        // required?
        main_pool.shutdown();


        return hits;
    }


    /**
     * screens building blocks of a single rxn.
     * fp must be fingerprint of m, SSS must be configured to search for m.
     *
     * @param space
     * @param cdh
     * @param m
     * @param fp
     * @param sss
     * @param rxn
     * @return
     */
    private static List<SynthonSpace.ExpandedHit> screen_building_blocks_for_rxn(SynthonSpace space , CachedDescriptorProvider cdh , StereoMolecule m, BitSet fp, SSSearcher sss, String rxn) {
        List<SynthonSpace.ExpandedHit> hits = new ArrayList<>();
        IDCodeParser icp = new IDCodeParser();

        String ri = rxn;
        for(int fi : space.ffps_sorted_by_rxn_and_frag_BT.get(ri).keySet()) {

            boolean found_hit_for_synthonset = false;
            FragType ft = new FragType(ri,fi);

            //BitSetTree bsti = this.ffps_sorted_by_rxn_and_frag_BT.get(ri).get(fi);
            BitSetTree bsti = space.ffps_sorted_by_rxn_and_frag_BT.get(ri).get(fi);

            BitSetTree.Node n_out[] = new BitSetTree.Node[1];
            bsti.testSubset(fp,n_out);
            List<BitSet> results = new ArrayList<>();
            bsti.root.collectSuperSets(fp,results);
            for(BitSet bsi : results) {
                //List<SynthonSpace.FragId> frags = this.fragments.get(bsi).stream().filter(si -> si.rxn_id.equals(ri) && si.frag==fi ).collect(Collectors.toList());
                List<SynthonSpace.FragId> frags = new ArrayList<>();
                try {
                    //frags = this.getAllFragments_byConnectors(getConnectorSetForFragType(ft)).get(bsi).stream().filter(si -> si.rxn_id.equals(ri) && si.frag == fi).collect(Collectors.toList());
                    frags = space.getAllFragments_byConnectors(space.getConnectorSetForFragType(ft)).get(bsi).stream().filter(si -> si.rxn_id.equals(ri) && si.frag == fi).collect(Collectors.toList());
                }
                catch(Exception ex) {
                    //ex.printStackTrace();
                    System.out.println("[Exception] Could not find frag for frag fp.. this is strange..");
                }
                for(SynthonSpace.FragId fid: frags) {
                    StereoMolecule sm = new StereoMolecule();
                    icp.parse(sm,fid.idcode);
                    // subset test..
                    sss.setMolecule(sm);
                    boolean test = sss.isFragmentInMolecule();
                    if(test) {
                        found_hit_for_synthonset = true;
                        System.out.println("Building Block Hit! -> Synthon Set "+ri+" : "+fi);
                        // create the full expanded hit..
                        Map<Integer, SynthonSpace.FragId> mapped_frag = new HashMap<>();
                        mapped_frag.put(fi,fid);
                        SynthonSpace.ExpandedHit hi = new SynthonSpace.ExpandedHit(ri,mapped_frag);
                        hits.add(hi);
                        break;
                    }
                }
                if(found_hit_for_synthonset) {
                    break;
                }
            }
        }

        return hits;
    }



    /**
     * All possible assigments of U,Np,Pu to the connectors.
     * num_splits different connectors are assigned
     *
     * @param m
     */
    private List<StereoMolecule> getAllConnectorAssignmentsForFragment(StereoMolecule m, int num_splits) {
        List<Integer> connector_pos = new ArrayList<>();
        for(int zi=0;zi<m.getAtoms();zi++) {
            if( m.getAtomicNo(zi)==92 || m.getAtomicNo(zi)==93 || m.getAtomicNo(zi)==94 ) {
                connector_pos.add(zi);
            }
        }

        // compute all assignments: 1. combinations, then all permutations for all combinations
        List<int[]> combis = CombinationGenerator.getAllOutOf(num_splits,connector_pos.size());
        List<List<Integer>> all_permutations = new ArrayList<>();
        for(int[] combi : combis) {
            List<Integer> clist = new ArrayList<>( Arrays.stream(combi).boxed().collect(Collectors.toList()) );
            List<List<Integer>> perms = all_permutations(clist);
            all_permutations.addAll(perms);
        }

        List<StereoMolecule> variants = new ArrayList<>();
        for(List<Integer> pi : all_permutations) {
            StereoMolecule mi = new StereoMolecule(m);
            for(int zi=0;zi<pi.size();zi++) {
                if(pi.get(zi)==0){mi.setAtomicNo(connector_pos.get(zi),92);}
                if(pi.get(zi)==1){mi.setAtomicNo(connector_pos.get(zi),93);}
                if(pi.get(zi)==2){mi.setAtomicNo(connector_pos.get(zi),94);}
            }
            mi.ensureHelperArrays(StereoMolecule.cHelperCIP);
            variants.add(mi);
        }
        return variants;
    }


    protected static List<List<Integer>> all_permutations(List<Integer> li) {
        if(li.size()==1) {
            ArrayList<List<Integer>> pi = new ArrayList<>();
            pi.add( new ArrayList<>(li) );
            return pi;
        }
        List<List<Integer>> all_permutations = new ArrayList<>();
        for(int zi=0;zi<li.size();zi++) {
            List<Integer> li2 = new ArrayList<>(li);
            Integer xi = li2.get(zi);
            li2.remove(zi);
            List<List<Integer>> perms = all_permutations(li2);
            for(List<Integer> pei : perms) {
                pei.add(xi);
                all_permutations.add(pei);
            }
        }
        return all_permutations;
    }


//    private FingerPrintGenerator fp_1         = new FingerPrintGenerator(false,true);
//    //private FragmentFingerPrintGenerator fp_2 = FragmentFingerPrintGenerator.getHyperFP(BITS);
//    private DescriptorHandlerLongFFP512  fp_3 = DescriptorHandlerLongFFP512.getDefaultInstance();
//    private DescriptorHandlerLongFFP1024_plus  fp_4 = DescriptorHandlerLongFFP1024_plus.getDefaultInstance();

    //private DescriptorHandler<long[],StereoMolecule> mFP = fp_4;
    //private DescriptorHandler<long[],StereoMolecule> mFP = DescriptorHandlerLongFFP1024_plus.getDefaultInstance();
    transient private DescriptorHandler<long[],StereoMolecule> mFP = null;
    private String mDescriptorHandlerShortName = "FFP1024_plus_ffp";

    transient Map<Thread,DescriptorHandler<long[],StereoMolecule>> mFPs_forThreads = new HashMap<>();

    public synchronized DescriptorHandler<long[],StereoMolecule> getDescriptorHandler_PerThread() {

        if(mFP==null) {
            // parse and restore the descriptor handler (based on mDescriptorHandlerShortName)
            this.mFP = resolveDescriptorHandlerFromName(mDescriptorHandlerShortName);
        }

        DescriptorHandler<long[],StereoMolecule> dhi = mFPs_forThreads.get( Thread.currentThread() );
        if(dhi==null) {
            dhi = mFP.getThreadSafeCopy();
            mFPs_forThreads.put(Thread.currentThread(),dhi);
        }

        return dhi;
    }

    /**
     * Now Threadsafe!
     *
     * @param mi
     * @return
     */
    public BitSet getFP(StereoMolecule mi) {
        return getFP( getDescriptorHandler_PerThread() , mi );
    }


    public static BitSet getFP( DescriptorHandler<long[],StereoMolecule> dh, StereoMolecule mi) {
        BitSet bs = BitSet.valueOf(dh.createDescriptor(mi));
        return bs;
    }

    public void setFP(DescriptorHandler<long[],StereoMolecule> dh, int bits) {
        this.mFP = dh;
        this.mDescriptorHandlerShortName = dh.getInfo().shortName;
        this.BITS = bits;
    }

//    /**
//     *
//     * @param m
//     * @return
//     */
//    public synchronized BitSet getFP(StereoMolecule m) {
//        //BitSet bs = fp_1.getFingerprint(m,BITS);
//        //BitSet bs =  BitSet.valueOf( fp_2.getFingerPrint(m) );
//        //BitSet bs = BitSet.valueOf(fp_4.createDescriptor(m));
//        BitSet bs = BitSet.valueOf(mFP.createDescriptor(m));
//        return bs;
//    }

    public DescriptorHandler<long[],StereoMolecule> getDescriptorHandler() {
        return mFP;
    }


//    Map<String,BitSet> cache_FPs = new HashMap<>();
//
//    /**
//     * TODO: improve multithreading. This way is pretty bad I guess.. Make FP computation ThreadPool / Service..
//     *
//     * @param m
//     * @return
//     */
//    public synchronized BitSet getFP_cached(StereoMolecule m) {
//        String hash = m.getIDCode();
//        if(cache_FPs.containsKey(hash)){
//            return cache_FPs.get(hash);
//        }
//        BitSet bs = getFP(m);
//        cache_FPs.put(hash,bs);
//        return bs;
//    }

    public String serializeToString() {
        StringBuilder sb = new StringBuilder();
        sb.append(""+BITS);
        sb.append("          ");
        sb.append(serialize_ffps_sorted_by_rxn_and_frag_BT());
        sb.append("          ");
        sb.append(serialize_fragments_to_string());
        sb.append("          ");
        sb.append(serialize_bsts_to_string());
        //sb.append( this.bst.serializeToString() );
        return sb.toString();
    }

    public void loadFromString(String data) {
        String splits[] = data.split("          ");
        int num_bits = Integer.parseInt(splits[0]);
        this.BITS = num_bits;
        init_ffps_sorted_by_rxn_and_frag_BT_fromString(splits[1]);
        init_fragments_from_string(splits[2]);
        init_bsts_from_string(splits[3]);
        //this.bst = new BitSetTree(splits[3]);
        //this.bsts_by_connectors = deserial

        // init rxn -> frag -> fragid map:
        this.fragment_map      = new HashMap<String,Map<Integer,List<FragId>>>();
        this.fragment_type_map = new HashMap<String,Map<Integer,FragType>>();

        //for( Set<SynthonSpace.FragId> fsi : this.fragments.values() ) {

        for(BitSet ki : this.fragments_by_connectors.keySet()) {
            for (Set<SynthonSpace.FragId> fsi : this.fragments_by_connectors.get(ki).values() ) {
                for (SynthonSpace.FragId fid : fsi) {
                    if (!this.fragment_map.containsKey(fid.rxn_id)) {
                        this.fragment_map.put(fid.rxn_id, new HashMap<>());
                    }
                    if (!this.fragment_map.get(fid.rxn_id).containsKey(fid.frag)) {
                        this.fragment_map.get(fid.rxn_id).put(fid.frag, new ArrayList());
                    }
                    this.fragment_map.get(fid.rxn_id).get(fid.frag).add(fid);

                    if (!this.fragment_type_map.containsKey(fid.rxn_id)) {
                        this.fragment_type_map.put(fid.rxn_id, new HashMap<>());
                    }
                    this.fragment_type_map.get(fid.rxn_id).put(fid.frag, new FragType(fid.rxn_id, fid.frag));
                }
            }
        }
        reinitHelperMaps();
    }


    private static DescriptorHandler<long[],StereoMolecule> resolveDescriptorHandlerFromName(String shortName) {
        if(shortName.equals("FFP1024_plus_ffp")) {
            return new DescriptorHandlerLongFFP1024_plus("ffp");
        }
        else if(shortName.equals("FFP1024_plus_pfp")) {
            return new DescriptorHandlerLongFFP1024_plus("pfp");
        }
        else if(shortName.equals("PPC2048")) {
            return new DescriptorHandlerPPCore();
        }
        else {
            return DescriptorHandlerStandard2DFactory.getFactory().create(shortName);
        }
    }

    /**
     * Initializes everything based on the fields:
     * fragments_by_connectors, ffps_sorted_by_rxn_and_frag_BT, bsts_by_connectors
     *
     * create the mFP based on the mDescriptorShortName
     *
     */
    public void initAfterJavaDeserialization() {
        mFPs_forThreads = new HashMap<>();

        mFP = resolveDescriptorHandlerFromName(mDescriptorHandlerShortName);

        this.fragment_map      = new HashMap<String,Map<Integer,List<FragId>>>();
        this.fragment_type_map = new HashMap<String,Map<Integer,FragType>>();

        //for( Set<SynthonSpace.FragId> fsi : this.fragments.values() ) {

        for(BitSet ki : this.fragments_by_connectors.keySet()) {
            for (Set<SynthonSpace.FragId> fsi : this.fragments_by_connectors.get(ki).values() ) {
                for (SynthonSpace.FragId fid : fsi) {
                    if (!this.fragment_map.containsKey(fid.rxn_id)) {
                        this.fragment_map.put(fid.rxn_id, new HashMap<>());
                    }
                    if (!this.fragment_map.get(fid.rxn_id).containsKey(fid.frag)) {
                        this.fragment_map.get(fid.rxn_id).put(fid.frag, new ArrayList());
                    }
                    this.fragment_map.get(fid.rxn_id).get(fid.frag).add(fid);

                    if (!this.fragment_type_map.containsKey(fid.rxn_id)) {
                        this.fragment_type_map.put(fid.rxn_id, new HashMap<>());
                    }
                    this.fragment_type_map.get(fid.rxn_id).put(fid.frag, new FragType(fid.rxn_id, fid.frag));
                }
            }
        }
        reinitHelperMaps();
    }

    private String serialize_bsts_to_string() {
        List<String> str_bsts = new ArrayList<>();
        for(BitSet ki : this.bsts_by_connectors.keySet()) {
            StringBuilder sbi = new StringBuilder();
            //sbi.append(serializeIntegerSet(ki));
            sbi.append(createBitSetString(ki,MAX_CONNECTORS));
            sbi.append("       ");
            sbi.append(this.bsts_by_connectors.get(ki).serializeToString());
            str_bsts.add(sbi.toString());
        }
        return String.join("        ",str_bsts);
    }

    private void init_bsts_from_string(String ds) {
        String[] bst_splits = ds.split("        ");
        for(String si : bst_splits) {
            String splits[] = si.split("       ");
            //Set<Integer> ki = deserializeIntegerSet(splits[0]);
            BitSet ki = parseBitSetString(splits[0]);
            BitSetTree bsti = new BitSetTree(splits[1]);
            this.bsts_by_connectors.put(ki,bsti);
        }
    }

    private String serialize_ffps_sorted_by_rxn_and_frag_BT() {
        List<String> strings = new ArrayList<>();
        for(String ki : this.ffps_sorted_by_rxn_and_frag_BT.keySet()) {
            for(Integer fi : this.ffps_sorted_by_rxn_and_frag_BT.get(ki).keySet() ) {
                StringBuilder sb = new StringBuilder();
                sb.append(ki);
                sb.append("     ");
                sb.append(fi);
                sb.append("     ");
                sb.append(this.ffps_sorted_by_rxn_and_frag_BT.get(ki).get(fi).serializeToString());
                strings.add(sb.toString());
            }
        }
        return String.join("      ",strings);
    }
    private void init_ffps_sorted_by_rxn_and_frag_BT_fromString(String s) {
        this.ffps_sorted_by_rxn_and_frag_BT = new HashMap<>();
        String entries[] = s.split("      ");
        for(String entry : entries) {
            String splits[] = entry.split("     ");

            String ki         = new String(splits[0]);
            int    fi         = Integer.parseInt(splits[1]);
            BitSetTree  bsti  = new BitSetTree(splits[2]);

            if(!this.ffps_sorted_by_rxn_and_frag_BT.containsKey(ki)){ this.ffps_sorted_by_rxn_and_frag_BT.put(ki,new HashMap<>()); }
            this.ffps_sorted_by_rxn_and_frag_BT.get(ki).put(fi,bsti);
        }
    }

    private String serialize_fragments_to_string() {
        List<String> x_si = new ArrayList<>();
        for( BitSet ki : this.fragments_by_connectors.keySet() ) {
            StringBuilder sb = new StringBuilder();
            Map<BitSet,Set<FragId>> fragments = this.getAllFragments_byConnectors(ki);
            //sb.append(serializeIntegerSet(ki));
            sb.append(createBitSetString(ki,MAX_CONNECTORS));
            sb.append("      ");
            List<String> si = new ArrayList<>();
            for (BitSet bsi : fragments.keySet()) {
                StringBuilder sb_i = new StringBuilder();
                sb_i.append(BitSetTree.bitSetToString(bsi));
                sb_i.append("    ");
                sb_i.append(serializeFragIdSetToString(fragments.get(bsi)));
                si.add(sb_i.toString());
            }
            sb.append( String.join("     ",si) );
            x_si.add(sb.toString());
        }
        return String.join("       ",x_si);
    }



    private void init_fragments_from_string(String si) {
        this.fragments_by_connectors = new HashMap<>();
        String frags[] = si.split("       ");
        for(String fi_a : frags) {
            String fa[] = fi_a.split("      ");
            //Set<Integer> ki = deserializeIntegerSet(fa[0]);
            BitSet ki = parseBitSetString(fa[0]);
            fragments_by_connectors.put(ki,new HashMap<>());

            String fa1_splits[] = fa[1].split("     ");
            for( String fa1_s : fa1_splits ) {
                String fi = fa1_s;
                String splits[] = fi.split("    ");
                BitSet bsi = BitSetTree.stringToBitSet(splits[0]);
                Set<FragId> fids = deserializeFragIdSetFromString(splits[1]);

                //this.fragments.put(bsi,fids);
                this.fragments_by_connectors.get(ki).put(bsi, fids);
            }
        }
    }

    private static String serializeIntegerSet(Set<Integer> s) {
        List<String> sli = s.stream().map(si -> si.toString()).collect(Collectors.toList());
        return String.join(";",sli);
    }
    public static Set<Integer> deserializeIntegerSet(String s) {
        String[] ss = s.split(";");
        Set<Integer> r = new HashSet<>();
        Arrays.stream(ss).forEach( si -> r.add(Integer.parseInt(si)) );
        return r;
    }


    public static String serializeFragIdToString(FragId id) {
        return null; // this has to be adapted, if you want to use the manual serialization
        //return id.rxn_id+" "+id.frag+" "+id.idcode+" "+id.fragment_id+" "+(id.has_blue?1:0)+" "+(id.has_red?1:0)+" "+(id.has_orange?1:0) + " "+ (id.has_green?1:0) + " "+ BitSetTree.bitSetToString(id.fp);
    }

    public static FragId deserializeFragIdFromString(String s) {
        return null; // this has to be adapted, if you want to use the manual serialization
//        String splits[] = s.split(" ");
//        return new FragId( new String(splits[0]),Integer.parseInt(splits[1]), new String(splits[2]) , new String(splits[3]) ,
//                Integer.parseInt(splits[4])==1 , Integer.parseInt(splits[5])==1 , Integer.parseInt(splits[6])==1 , Integer.parseInt(splits[7])==1 ,
//                BitSetTree.stringToBitSet(splits[8]+" "+splits[9]));
    }

    public static String serializeFragIdSetToString(Set<FragId> ids) {
        if(ids.isEmpty()) {
            return "<<EMPTY>>";
        }
        List<String> fs = new ArrayList<>();
        for(FragId fid : ids) {
            fs.add( serializeFragIdToString(fid) );
        }
        return String.join("  ",fs);
    }

    public static Set<FragId> deserializeFragIdSetFromString(String s) {
            if(s.equals("<<EMPTY>>")){
            return new HashSet<>();
        }
        Set<FragId> set = new HashSet<>();
        String splits[] = s.split("  ");
        for(String si : splits) {
            set.add( deserializeFragIdFromString(si) );
        }
        return set;
    }


    public BitSet getFPfromSmiles(String smiles) throws Exception {
        SmilesParser sp   = new SmilesParser();
        StereoMolecule sm = new StereoMolecule();
        sp.parse(sm,smiles);
        //FingerPrintGenerator fpgen = new FingerPrintGenerator();
        //BitSet bs_fp_i = fpgen.getFingerprint(sm,512);
        return getFP( getDescriptorHandler_PerThread() , sm );
    }

}
