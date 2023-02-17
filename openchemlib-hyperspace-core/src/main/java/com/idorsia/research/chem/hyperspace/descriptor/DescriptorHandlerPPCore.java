package com.idorsia.research.chem.hyperspace.descriptor;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.descriptor.AbstractDescriptorHandlerLongFP;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorInfo;
import com.actelion.research.chem.descriptor.ISimilarityCalculator;
import com.actelion.research.chem.phesa.pharmacophore.*;
import com.actelion.research.chem.phesa.pharmacophore.pp.AcceptorPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.ChargePoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.DonorPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import java.io.IOException;
import java.util.*;

public class DescriptorHandlerPPCore extends AbstractDescriptorHandlerLongFP<StereoMolecule> {



    //private static final int MAX_BITS = 512;
    private static final int MAX_DEPTH = 6;
    private static final boolean DEBUG = false;
    private static int debugCounter = 0;
    private Hashtable paths;

    boolean considerNothingButTriple = true;
    boolean considerDelocalized      = false;

    private static DescriptorHandlerPPCore defaultInstance = new DescriptorHandlerPPCore();

    public static DescriptorHandlerPPCore getDefaultInstance() {
        return defaultInstance;
    }


    public static class AtomInfo {
        public final int AtomicNo;
        public final int isConnector; // 0 means no, 1,2,3 means is connector
        public final boolean isAromatic;
        public final boolean isInRing;

        public final int pharmacophore; // -1 means no, 0,1,2,3,4,5,6 means yes..

        public AtomInfo(int atomic_no, int connector, boolean aromatic, boolean in_ring, int pp) {
            this.AtomicNo=atomic_no;
            this.isConnector = connector;
            this.isAromatic = aromatic;
            this.isInRing = in_ring;
            this.pharmacophore = pp;
        }
    }

    public static class BondInfo {
        public final int order;
        public final boolean isDelocalized;
        public final boolean isInRing;
        public final boolean isRotatable;
        public final int ringSize;

        public BondInfo(int order, boolean delocalized, boolean in_ring, int ring_size, boolean is_rotatable) {
            this.order = order;
            this.isDelocalized = delocalized;
            this.isInRing = in_ring;
            this.ringSize = ring_size;
            this.isRotatable = is_rotatable;
        }
    }


//    @Override
//    public String encode(long[] o) {
//        throw new Error("not implemented");
//    }
//
//    @Override
//    public long[] decode(String s) {
//        throw new Error("not implemented");
//    }
//
//    @Override
//    public long[] decode(byte[] bytes) {
//        throw new Error("not implemented");
//    }
//
//    @Override
//    public boolean calculationFailed(long[] o) {
//        throw new Error("not implemented");
//    }
//
//    @Override
//    public float getSimilarity(long[] d1, long[] d2) {
//        throw new Error("not implemented");
//    }

    @Override
    public DescriptorInfo getInfo() {
        return new DescriptorInfo("PPCore","PPC2048","-1",0,true,true,false,false);
    }

    @Override
    public String getVersion() {
        throw new Error("not implemented");    }

    @Override
    public long[] createDescriptor(StereoMolecule m_a) {

        StereoMolecule m = new StereoMolecule(m_a);
        m.ensureHelperArrays(Molecule.cHelperCIP);

        int bits = 512; int max_multiplicity = 4;

        // small ones get special treatment
        // they also get a bit pattern that is very far away from all regular ones
        // DEACTIVATED (may not necessary?)
        // probably if we activate it we should just set these bits in addition to the regular ones?
        if(m.getAtoms()<=-1) {
            BitSet bi = new BitSet(bits*max_multiplicity);
            for(int zi=0;zi<bits;zi+=6) {
                bi.set(zi*max_multiplicity+1,true);
            }
            List<Integer> all_atoms = new ArrayList<>();
            for(int zi=0;zi<m.getAtoms();zi++) { all_atoms.add(m.getAtomicNo(zi)); }
            // mark connectors
            if( all_atoms.contains(92) || all_atoms.contains(93) || all_atoms.contains(94) ) {
                for (int zi = 1; zi < bits; zi += 3) {
                    bi.set(zi * max_multiplicity + 2, true);
                }
            }
            // mark halogens
            if( all_atoms.contains(9) || all_atoms.contains(17) || all_atoms.contains(35) || all_atoms.contains(53) ) {
                for (int zi = 0; zi < bits/4; zi += 2) {
                    bi.set(zi * max_multiplicity + 1, true);
                }
            }
            // mark O
            if( all_atoms.contains(8) ) {
                for (int zi = 0; zi < bits/4; zi += 5) {
                    bi.set(zi * max_multiplicity + 1, true);
                }
            }
            // mark S
            if( all_atoms.contains(16) ) {
                for (int zi = 0; zi < bits/2; zi += 7) {
                    bi.set(zi * max_multiplicity + 2, true);
                }
            }

            long[] data = Arrays.copyOf( bi.toLongArray() , bits*max_multiplicity / 64 );
            //System.out.println("ok_small!!");
            return data;
        }


        List<IPharmacophorePoint> pp = new ArrayList<>();
        try {
            StereoMolecule m_2 = new StereoMolecule(m);
            m_2.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
            for(int zi=0;zi<m_2.getAtoms();zi++){
                if( m_2.getAtomicNo(zi)==92 || m_2.getAtomicNo(zi)==93 || m_2.getAtomicNo(zi)==94 ) {
                    m_2.setAtomicNo(zi,6);
                }
            }
            m_2.ensureHelperArrays(Molecule.cHelperCIP);
            pp = PharmacophoreCalculator.getPharmacophorePoints(m_2);
        }
        catch(Exception ex) {
            ex.printStackTrace();
            System.out.println("Exception in PharmacophoreCalculator for Molecule: "+ HyperspaceUtils.idcodeToSmiles(m.getIDCode()));
        }

        boolean[] rotatable_bonds = new boolean[m.getBonds()];
        TorsionDB.findRotatableBonds(m, true, rotatable_bonds);
        this.rotatable_bonds = rotatable_bonds;

        this.p_points = new int[m.getAtoms()];
        for(int zi=0;zi<m.getAtoms();zi++){ this.p_points[zi]=-1;}
        for(IPharmacophorePoint pi : pp) {
            if(pi instanceof AcceptorPoint) { this.p_points[((AcceptorPoint)pi).getCenterID()] = PharmacophoreCalculator.ACCEPTOR_ID; }
            if(pi instanceof DonorPoint) { this.p_points[((DonorPoint)pi).getCenterID()] = PharmacophoreCalculator.DONOR_ID; }
            if(pi instanceof ChargePoint) {
                if( ((ChargePoint)pi).getCharge()>0 ){ this.p_points[((DonorPoint)pi).getCenterID()] = PharmacophoreCalculator.CHARGE_POS_ID;}
                else if(((ChargePoint)pi).getCharge()<0 ){ this.p_points[ ((DonorPoint)pi).getCenterID()] = PharmacophoreCalculator.CHARGE_NEG_ID;}
            }
        }

        //System.out.println("ok");

        // now compute dfs paths for every p point as well as for every connector:
        Map<Integer,List<MolPath>> paths = new HashMap<>();

        for(int zi=0;zi<m.getAtoms();zi++) {
            int atomic_number = m.getAtomicNo(zi);
            if( atomic_number == 92 || atomic_number == 93 || atomic_number == 94 ) {
                // start connector paths:
                List<MolPath> paths_i = new ArrayList<>();
                boolean visited[] = new boolean[m.getAtoms()];
                MolPath mpi = new MolPath();
                //createPathFingerprintFeature_df_recursive(m,-1,zi,8,mpi,paths_i,visited);
                createPathFingerprintFeature_df_recursive(m,-1,zi,12,mpi,paths_i,visited);
                paths.put(zi,paths_i);
            }
            else if( this.p_points[zi] >= 0 ) {
                // start connector paths:
                List<MolPath> paths_i = new ArrayList<>();
                boolean visited[] = new boolean[m.getAtoms()];
                MolPath mpi = new MolPath();
                //createPathFingerprintFeature_df_recursive(m,-1,zi,8,mpi,paths_i,visited);
                createPathFingerprintFeature_df_recursive(m,-1,zi,12,mpi,paths_i,visited);
                paths.put(zi,paths_i);
            }
        }


        // now analyze paths:

        // 1. just raw info, i.e. length, flex,
        List<String> hashs_raw = new ArrayList<>();
        for(int ki : paths.keySet()) {
            for(MolPath mpi : paths.get(ki)) {
                mpi.computeProperties();

                // 1. determine all kinds of predicates
                String size = null;

                if(mpi.atoms.size()<4) {
                    size = "s";
                }
                if(mpi.atoms.size()>4 && mpi.atoms.size() < 7) {
                    size = "m";
                }
                if(mpi.atoms.size()>7) {
                    size="l";
                }

                String aromaticity = null;
                if( (1.0*mpi.delocalized_bonds)/mpi.atoms.size() < 0.1 ) {
                    aromaticity = "al";
                }
                else if((1.0*mpi.delocalized_bonds)/mpi.atoms.size() < 0.3 ) {
                    aromaticity = "am";
                }
                else {
                    aromaticity = "ah";
                }

                String flex = null;
                if( (1.0*mpi.flexibility) / mpi.atoms.size() < 0.15 ) {
                    flex = "fl";
                }
                else if((1.0*mpi.flexibility) / mpi.atoms.size() <= 0.25 ) {
                    flex = "fm";
                }
                else {
                    flex = "fh";
                }


                // special treatment when it actually contains at least two ppoints:
                String pp_class = "p<2";
                if(mpi.furthest_special_points_distance>0) {
                    pp_class = "p>2";
                }

                // special treatment for long paths:
                String flex_2 = "xs";
                if( mpi.furthest_special_points_distance > 6 && mpi.flexibility < 2 ) {
                    flex_2 = "xls";
                }
                else if (mpi.furthest_special_points_distance > 6 && mpi.flexibility >= 3) {
                    flex_2 = "xlf";
                }


                String med_dist_class = "md_0";
                // special treatment for medium dist:
                if( mpi.furthest_special_points_distance>3 && mpi.furthest_special_points_distance<6 && mpi.flexibility <2) {
                    med_dist_class = "md_1";
                }
                if( mpi.furthest_special_points_distance>3 && mpi.furthest_special_points_distance<6 && mpi.flexibility >=2) {
                    med_dist_class = "md_2";
                }


                // 2. some combinatorial stuff from predicates:
                List<String> predicates_point_types = new ArrayList<>();
                if(mpi.p_n_acceptor>0){predicates_point_types.add("ac");}
                if(mpi.p_n_donor>0){predicates_point_types.add("do");}
                if(mpi.p_n_charge_neg>0){predicates_point_types.add("ne");}
                if(mpi.p_n_charge_pos>0){predicates_point_types.add("po");}
                if(mpi.p_n_connector>0){predicates_point_types.add("cn");}

                List<String> predicates_size = new ArrayList<>();
                if(mpi.atoms.size()<4){predicates_size.add("0x<4");}
                if(mpi.atoms.size()>4){predicates_size.add("1x>4");}
                if(mpi.atoms.size()>6){predicates_size.add("2y>6");}
                if(mpi.atoms.size()<7){predicates_size.add("3z<7");}
                if(mpi.furthest_special_points_distance<4){predicates_size.add("4fd<4");}
                if(mpi.furthest_special_points_distance<6){predicates_size.add("5fd<6");}
                if(mpi.furthest_special_points_distance>6){predicates_size.add("5fd>6");}

                List<String> predicates_flex = new ArrayList<>();
                if(mpi.flexibility>0){ predicates_flex.add("0fx>0"); }
                if(mpi.flexibility<2){ predicates_flex.add("1fx<2"); }
                if(mpi.flexibility>2){ predicates_flex.add("2fx>2"); }

                List<String> predicates_ring = new ArrayList<>();
                if(mpi.ring_bonds>0){predicates_size.add("ar>0");}
                if(mpi.ring_bonds>2){predicates_size.add("br>2");}
                if(mpi.ring_bonds>4){predicates_size.add("cr>4");}
                if(mpi.delocalized_bonds>0){predicates_size.add("dcl>0");}
                if(mpi.delocalized_bonds<2){predicates_size.add("ecl<2");}
                if(mpi.delocalized_bonds>2){predicates_size.add("fcl>2");}

                List<String> predicates_ring_2 = new ArrayList<>(predicates_ring);
                if(mpi.has_ring_3||mpi.has_ring_4){predicates_ring_2.add("smrxx");}
                if(mpi.has_ring_3){predicates_ring_2.add("smr3");}
                if(mpi.has_ring_4){predicates_ring_2.add("smr4");}
                if(mpi.has_ring_5){predicates_ring_2.add("smr5");}
                if(mpi.has_ring_6){predicates_ring_2.add("smr6");}
                if(mpi.has_ring_7){predicates_ring_2.add("smr7");}
                if(mpi.has_ring_greater_7){predicates_ring_2.add("smrgx7");}



                String predicates_size_combined = String.join("x",predicates_size);
                String predicates_pt_combined   = String.join("x",predicates_point_types);
                String predicates_ring_combined = String.join("r",predicates_ring);
                String predicates_ring_2_combined = String.join("q",predicates_ring_2);
                String predicates_flex_combined = String.join("f",predicates_flex);

                // extra structural info for small paths..
                String extra_1 = (mpi.atoms.size()<=3)? mpi.bonds.stream().map(bi -> bi.isRotatable+"_"+bi.isDelocalized ).reduce((x,y)->x+";"+y).toString():"L";
                String extra_2 = (mpi.atoms.size()<=4)? mpi.bonds.stream().map(bi -> bi.isRotatable+"_"+bi.isDelocalized ).reduce((x,y)->x+";"+y).toString():"L";


                for(String psi : predicates_size) {
                    //hashs_raw.add("h1"+psi+":"+predicates_pt_combined);
                    //hashs_raw.add("h1"+psi+":"+predicates_ring_combined+":"+predicates_pt_combined);
                    //hashs_raw.add("h1"+psi+":"+predicates_flex_combined+":"+predicates_pt_combined);
                    hashs_raw.add("h1"+psi+":"+predicates_ring_2_combined+":"+predicates_flex_combined+":"+predicates_pt_combined+"_"+flex+"_"+aromaticity);
                    if(!pp_class.equals("p<2")) {
                        // special treatment for paths with pps at start and end
                        if(mpi.start_end_class[0]>-1 && mpi.start_end_class[0]>-1) {
                            hashs_raw.add("px_" + psi +"_"+ mpi.start_end_class[0] + "::" + mpi.start_end_class[1] + "_"+extra_1);
                            hashs_raw.add("px_" + psi +"_"+ flex + "_"+ flex_2 +"_"+mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ flex_2 +"_" + mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ flex_2 +"_" + predicates_ring_2_combined +"_"+  mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ flex_2 + "_ " +predicates_ring_combined+"_"+mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ aromaticity + mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ flex + mpi.start_end_class[0] + "::" + mpi.start_end_class[1]);
                            hashs_raw.add("px_" + psi +"_"+ flex + mpi.start_end_class[0] + "::" + mpi.start_end_class[1] + "->"+ mpi.atoms.get(0).AtomicNo + "_"+mpi.atoms.get(mpi.atoms.size()-1).AtomicNo);
                        }
                    }

                    for(String ri : predicates_ring_2) {
                        hashs_raw.add("xh1"+psi+":"+predicates_pt_combined+";"+ri+"_"+flex+"_"+predicates_ring_2_combined+"_"+aromaticity+"_"+extra_2);
                        //hashs_raw.add("xh1"+psi+":"+predicates_pt_combined+";"+ri+"_"+aromaticity);
                    }
                }

                // special treatment for double pp paths
                if(predicates_point_types.size()==2) {
                    hashs_raw.add("2p_y_"+predicates_size_combined+":"+predicates_pt_combined+"x"+flex+"_"+flex_2);
                    hashs_raw.add("2p_y_"+predicates_size_combined+":"+predicates_pt_combined+"_"+flex+"_"+aromaticity);
                    hashs_raw.add("2p_y_"+predicates_size_combined+":"+predicates_pt_combined+"_"+flex_2+"_"+aromaticity);
                    for(String pi : predicates_ring_2) {
                        hashs_raw.add("2p_y_"+pi+"_"+predicates_size_combined+":"+predicates_pt_combined);
                        hashs_raw.add("2p_y_"+pi+"x"+predicates_flex_combined+":"+predicates_pt_combined);
                    }
                }

                // special treatment for (at least) triple pp paths
                if(predicates_point_types.size()>=3) {
                    hashs_raw.add("trp__h1"+predicates_size_combined+":"+predicates_pt_combined);
                    for(String ppi : predicates_point_types) {
                        hashs_raw.add("trp_x1"+ppi+predicates_size_combined+":"+predicates_pt_combined);
                        hashs_raw.add("trp_x2"+ppi+predicates_size_combined+":"+predicates_pt_combined+":"+predicates_flex_combined);
                        hashs_raw.add("trp_x3"+ppi+predicates_size_combined+":"+predicates_pt_combined+":"+predicates_ring_combined);
                    }
                }

            }
        }

        //Set<String> all_hashs = new HashSet<>(hashs_raw);
        //System.out.println("all_hashs: "+all_hashs.size() + " in list: "+hashs_raw.size());
        //System.out.println("ok");

        // note, for all that exist more than once there exist symmetric ones, so we count (x<2)?x:(1+x/2)
        int counts[] = new int[bits];
        for(String si : hashs_raw) {
            int hash = (new Random(si.hashCode())).nextInt(bits);
            counts[hash] += 1;
//            if(counts[hash]==30){
//                System.out.println(si);
//            }
        }
        //System.out.println("ok!");
        BitSet bsi = new BitSet(bits*max_multiplicity);
        for(int zi=0;zi<bits;zi++) {
            int ci = (counts[zi]<2)?counts[zi]:(1+counts[zi]/2);
            ci     = Math.min(ci,4);
            // set bits:
            for(int bi=0;bi<ci;bi++) {
                bsi.set( zi*max_multiplicity + bi , true );
            }
        }

        //if(bsi.cardinality()>1100){ System.out.println("cardinality="+bsi.cardinality() + " num_atoms="+m.getAtoms());}

        //long[] data = new long[bits*max_multiplicity / 64];
        long[] data = Arrays.copyOf( bsi.toLongArray() , bits*max_multiplicity / 64 );
        //System.out.println("ok!!");

        return data;
    }

    @Override
    public DescriptorHandler<long[], StereoMolecule> getThreadSafeCopy() {
        return new DescriptorHandlerPPCore();
    }



    int p_points[] = null;
    boolean rotatable_bonds[]     = null;


    public final static class MolPath {

        public MolPath(){}

        /**
         * Copy constructor
         * @param mp
         */
        public MolPath(MolPath mp) {
            this.atoms.addAll(mp.atoms);
            this.bonds.addAll(mp.bonds);
        }

        public final List<AtomInfo> atoms = new ArrayList<>(16);
        public final List<BondInfo> bonds = new ArrayList<>(16);

        public void addAtomInfo(AtomInfo ai){this.atoms.add(ai);}
        public void addBondInfo(BondInfo bi){this.bonds.add(bi);}


        public int flexibility  = 0;
        public int ring_bonds   = 0;
        public int delocalized_bonds = 0;

        public boolean has_ring_3 = false;
        public boolean has_ring_4 = false;
        public boolean has_ring_5 = false;
        public boolean has_ring_6 = false;
        public boolean has_ring_7 = false;
        public boolean has_ring_greater_7 = false;

        public int p_n_connector   = 0;

        public int p_n_acceptor    = 0;
        public int p_n_donor       = 0;
        public int p_n_charge_pos  = 0;
        public int p_n_charge_neg  = 0;

        public int furthest_special_points_distance = 0;


        public int p_n_all = 0;


        // an identifier classifying the path based on the start and endpoint
        // NOTE: we canonize the order (nat. order of identifiers)
        public int[] start_end_class = new int[]{0,0};

        /**
         * Computes a number of properties of this MolPath
         */
        public void computeProperties() {

            this.flexibility = 0;
            this.ring_bonds  = 0;
            this.delocalized_bonds = 0;
            for(BondInfo bi : bonds){
                flexibility += (bi.isRotatable)?1:0;
                ring_bonds  += (bi.isInRing)?1:0;
                delocalized_bonds += (bi.isDelocalized)?1:0;
                if(bi.ringSize==3){has_ring_3=true;}
                if(bi.ringSize==4){has_ring_4=true;}
                if(bi.ringSize==5){has_ring_5=true;}
                if(bi.ringSize==6){has_ring_6=true;}
                if(bi.ringSize==7){has_ring_7=true;}
                if(bi.ringSize>7){has_ring_greater_7=true;}
            }

            this.p_n_connector = 0;
            this.p_n_acceptor = 0;
            this.p_n_donor    = 0;
            this.p_n_charge_neg = 0;
            this.p_n_charge_pos = 0;

            for(int zi=0;zi<this.atoms.size();zi++) {
                AtomInfo ai = atoms.get(zi);
                if ( (ai.isConnector>0) || ai.pharmacophore >= 0) {
                    furthest_special_points_distance = zi;
                    p_n_connector += (ai.isConnector > 0) ? 1 : 0;
                    p_n_acceptor += (ai.pharmacophore == PharmacophoreCalculator.ACCEPTOR_ID) ? 1 : 0;
                    p_n_donor += (ai.pharmacophore == PharmacophoreCalculator.DONOR_ID) ? 1 : 0;
                    p_n_charge_neg += (ai.pharmacophore == PharmacophoreCalculator.CHARGE_NEG_ID) ? 1 : 0;
                    p_n_charge_pos += (ai.pharmacophore == PharmacophoreCalculator.CHARGE_POS_ID) ? 1 : 0;
                    p_n_all += 1;
                }
            }

            // compute start end class:
            int class_first  = this.atoms.get(0).pharmacophore;
            int class_second = this.atoms.get(atoms.size()-1).pharmacophore;
            int se_class[] = new int[]{ class_first,class_second};
            Arrays.sort(se_class);
            this.start_end_class = se_class;

        }

    }




    private void createPathFingerprintFeature_df_recursive(StereoMolecule m,
                                                                  int last_pos, int pos,
                                                                  int remaining_steps, MolPath current_path,
                                                                  List<MolPath> paths, boolean visited[]) {
        if(remaining_steps==0){return;};

        MolPath next_path = new MolPath(current_path);
        if(last_pos>=0) {
            int bond_i = m.getBond(last_pos,pos);
            int bond_order = m.getBondOrder(bond_i);
            boolean bond_rotatable = this.rotatable_bonds[bond_i];
            int ring_size = m.getBondRingSize(bond_i);
            next_path.addBondInfo(new BondInfo(bond_order,m.isDelocalizedBond(bond_i),m.isRingBond(bond_i),ring_size,bond_rotatable));
        }

        int atomic_no = m.getAtomicNo(pos);
        int connector = 0;
        if(atomic_no >= 92 && atomic_no <= 94) { connector = atomic_no - 91; }
        AtomInfo atom_info = new AtomInfo(atomic_no,connector,m.isAromaticAtom(pos),m.isRingAtom(pos),this.p_points[pos]);
        next_path.addAtomInfo(atom_info);
        paths.add(next_path);
        visited[pos] = true;

        for( int zi=0;zi<m.getConnAtoms(pos);zi++) {
            boolean visited_next[] = Arrays.copyOf(visited,visited.length);
            int next_pos = m.getConnAtom(pos,zi);
            if(visited[next_pos]){continue;}
            createPathFingerprintFeature_df_recursive(m,pos,next_pos,remaining_steps-1,next_path,paths,visited_next);
        }
    }

    public static void main(String args[]) {
        test_01();
        test_02();
    }

    public static void test_02() {
        try {
            List<StereoMolecule> list = HyperspaceUtils.parseInputFile("C:\\datasets\\Osiris_Project_Profile_CCR8_c.txt");
            List<Integer> cardinalities = new ArrayList<>();
            long ts_a = System.currentTimeMillis();
            int cnt = 0;;
            DescriptorHandlerPPCore ppc = new DescriptorHandlerPPCore();
            for(StereoMolecule mi : list) {
                cardinalities.add(BitSet.valueOf(ppc.createDescriptor(mi)).cardinality());
                cnt++;
            }
            long ts_b = System.currentTimeMillis();
            System.out.println("Structures: "+ cnt);
            System.out.println("Time: "+(ts_b-ts_a));
            System.out.println("Cardinalities avg: "+cardinalities.stream().mapToInt(vi -> vi.intValue()).summaryStatistics().getAverage() );

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void test_01() {

        String smiles_tiny        = "CCc1ccccc1";
        String smiles_test        = "";

        String smiles_01          = "";
        String smiles_01b         = "";
        String smiles_01_complete = "";
        String smiles_02          = "";
        String smiles_02_complete = "";


        SmilesParser sp = new SmilesParser();

        StereoMolecule mtiny  = new StereoMolecule();

        StereoMolecule mtest  = new StereoMolecule();

        StereoMolecule m01    = new StereoMolecule();
        StereoMolecule m01b   = new StereoMolecule();
        StereoMolecule m01_c  = new StereoMolecule();

        StereoMolecule m02   = new StereoMolecule();
        StereoMolecule m02_c = new StereoMolecule();

        try{
            sp.parse(mtiny,smiles_tiny);
            sp.parse(mtest,smiles_test);
            sp.parse(m01,smiles_01);
            sp.parse(m01b,smiles_01b);
            sp.parse(m01_c,smiles_01_complete);
            sp.parse(m02,smiles_02);
            sp.parse(m02_c,smiles_02_complete);
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }



        DescriptorHandlerPPCore ppc = new DescriptorHandlerPPCore();

        BitSet bstest = BitSet.valueOf(ppc.createDescriptor(mtest));
        BitSet bs01   = BitSet.valueOf(ppc.createDescriptor(m01));
        BitSet bs02   = BitSet.valueOf(ppc.createDescriptor(m02));


        System.out.println("ok");
    }

}
