package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.*;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

public class SynthonAssembler {
    static int       logLevel_assemble = 0;

    public static StereoMolecule assembleSynthons_faster(List<StereoMolecule> m_parts) {
        return assembleSynthons_faster(m_parts,null);
    }

    /**
     *
     *
     * @param m_parts
     * @param atom_maps maps atoms in m_parts to atoms in assembled fragment
     * @return
     */
    public static StereoMolecule assembleSynthons_faster(List<StereoMolecule> m_parts, List<int[]> atom_maps) {
        if(logLevel_assemble>0) {
            System.out.println("Fragments to assemble:");
            for(int zi=0;zi<m_parts.size();zi++) {
                System.out.println( idcodeToSmiles(m_parts.get(zi).getIDCode()));
            }
        }

        if(logLevel_assemble>0) {
            System.out.println( HyperspaceUtils.idcodeToSmiles(m_parts.get(0).getIDCode()) );
        }

        // copy all into assembled_pre
        StereoMolecule assembled_pre = new StereoMolecule();//StereoMolecule(m_parts.get(0));


        //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
        //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours);

        for(int zi=0;zi<m_parts.size();zi++) {
            // assemble ..

            m_parts.get(zi).ensureHelperArrays(Molecule.cHelperCIP);

            // for the addFragment function to work correctly for query features, we have to
            // set isFragment in assembled_pre!
            if(m_parts.get(zi).isFragment()) { assembled_pre.setFragment(true); }


            int atom_map[] = new int[m_parts.get(zi).getAtoms()];
            assembled_pre.addFragment(m_parts.get(zi), 0, atom_map);
            //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
            //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours);

            if(atom_maps!=null) {
                atom_maps.add(atom_map);
            }

            if(logLevel_assemble>0) {
                System.out.println( HyperspaceUtils.idcodeToSmiles(assembled_pre.getIDCode()) );
            }
        }

        //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
        assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours);


        // now in one go: add bonds to all connector neighbors..
        for(int zc=92;zc<=95;zc++) {
            int pos_connector_a = -1;
            int pos_connector_b = -1;
            for (int zi = 0; zi < assembled_pre.getAtoms(); zi++) {
                if(assembled_pre.getAtomicNo(zi)==zc) {
                    if(pos_connector_a<0){pos_connector_a=zi;}
                    else{pos_connector_b=zi;}
                }
            }
            int num_connectors_found = ((pos_connector_a>=0)?1:0) + ((pos_connector_b>=0)?1:0);
            if(num_connectors_found==0) {
                // not found, thats ok!
            }
            else if(num_connectors_found==1) {
                // single connector found, this might also be ok (e.g. in sss, when assembling combinatorial hits)
            }
            else if(num_connectors_found==2) {
                // asssemble!
                int conn_atom_a = assembled_pre.getConnAtom(pos_connector_a, 0);
                int conn_atom_b = assembled_pre.getConnAtom(pos_connector_b, 0);

                int bond_a = assembled_pre.getConnBond(pos_connector_a, 0);
                int bond_b = assembled_pre.getConnBond(pos_connector_a, 0);
                if (assembled_pre.getBondType(bond_a) != assembled_pre.getBondType(bond_b)) {
                    System.out.println("[ERROR] error while assembling fragments -> different bond types?");
                }
                assembled_pre.addBond(conn_atom_a, conn_atom_b, assembled_pre.getBondType(bond_a));

                //m_parts.get(0).ensureHelperArrays(StereoMolecule.cHelperCIP); // new
                //assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours); // new

                // remove Us
                assembled_pre.markAtomForDeletion(pos_connector_a);
                assembled_pre.markAtomForDeletion(pos_connector_b);
            }
        }

        int atom_map_deletion[] = assembled_pre.deleteMarkedAtomsAndBonds();
        assembled_pre.ensureHelperArrays(Molecule.cHelperCIP);

        if(atom_maps!=null) {
            // in this case we have to adjust atom positions due to the
            // call to deleteMarkedAtomsAndBonds()
            for(int[] atom_map_i : atom_maps) {
                for(int zi=0;zi<atom_map_i.length;zi++) {
                    // atom_map_i is (frag->assembled).
                    // atom_map_deletion is (assembled->assembled_after_deletion), i.e.
                    atom_map_i[zi] = atom_map_deletion[ atom_map_i[zi] ];
                }
            }
        }

        return assembled_pre;
    }

    public static StereoMolecule assembleSynthons(List<StereoMolecule> m_parts) {

        if(logLevel_assemble>0) {
            System.out.println("Fragments to assemble:");
            for(int zi=0;zi<m_parts.size();zi++) {
                System.out.println( idcodeToSmiles(m_parts.get(zi).getIDCode()));
            }
        }

        if(logLevel_assemble>0) {
            System.out.println( HyperspaceUtils.idcodeToSmiles(m_parts.get(0).getIDCode()) );
        }

        if(m_parts.size()==0) {
            return new StereoMolecule();
        }

        // copy all into the m_parts.get(0)
        m_parts.get(0).ensureHelperArrays(Molecule.cHelperCIP);
        StereoMolecule assembled_pre = new StereoMolecule(m_parts.get(0));
        assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
        assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours);

        //for(int zi=1;zi<m_parts.size();zi++) {
        for(int zi=1;zi<m_parts.size();zi++) {
            // assemble ..
            int atom_map[] = new int[200];

            m_parts.get(zi).ensureHelperArrays(Molecule.cHelperCIP);

//            m_parts.get(0).addFragment(m_parts.get(zi), 0, atom_map);
//            m_parts.get(0).ensureHelperArrays(StereoMolecule.cHelperCIP);
//            m_parts.get(0).ensureHelperArrays(StereoMolecule.cHelperNeighbours);
            assembled_pre.addFragment(m_parts.get(zi), 0, atom_map);
            assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
            assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours);

            if(logLevel_assemble>0) {
                System.out.println( HyperspaceUtils.idcodeToSmiles(assembled_pre.getIDCode()) );
            }

        }

//        // check if we do not have to connect anything. Then we break early
//        if(connectors_u<2 && connectors_np<2 && connectors_pu<2 ) {
//            if(logLevel_assemble>0){ System.out.println("Nothing to assemble -> return");}
//            StereoMolecule m_assembled = new StereoMolecule(m_parts.get(0));
//            return m_assembled;
//        }



        // now connect uranium atoms with same color:
        boolean created_connection = false;
        boolean still_connecting = true;

        while(still_connecting) {
            still_connecting = false;
            created_connection = false;
            for (int zi = 0; zi < assembled_pre.getAtoms() - 1; zi++) {
                if ( assembled_pre.getAtomicNo(zi) != 92 && assembled_pre.getAtomicNo(zi) != 93 && assembled_pre.getAtomicNo(zi) != 94 && assembled_pre.getAtomicNo(zi) != 95 ) {
                    continue;
                }
                else {
                    for (int zj = zi + 1; zj < assembled_pre.getAtoms(); zj++) {
                        //if (  assembled_pre.getAtomicNo(zi)==assembled_pre.getAtomicNo(zj) && ((assembled_pre.getAtomicNo(zj) == 92 || assembled_pre.getAtomicNo(zj) == 93 || assembled_pre.getAtomicNo(zj) == 94) || assembled_pre.getAtomicNo(zj) == 95) ) {
                        if (  assembled_pre.getAtomicNo(zi)==assembled_pre.getAtomicNo(zj) ) {
                            if(logLevel_assemble>0) {
                                System.out.println("Assemble connector "+assembled_pre.getAtomicNo(zi));
                            }

                            // connect uranium neighbor atoms with the right bond type, and then remove uranium..
                            if (assembled_pre.getConnAtoms(zi) != 1 || assembled_pre.getConnAtoms(zj) != 1) {
                                System.out.println("error while assembling fragments -> too many neighbors?");
                            }
                            int conn_atom_a = assembled_pre.getConnAtom(zi, 0);
                            int conn_atom_b = assembled_pre.getConnAtom(zj, 0);

                            int bond_a = assembled_pre.getConnBond(zi, 0);
                            int bond_b = assembled_pre.getConnBond(zj, 0);
                            if (assembled_pre.getBondType(bond_a) != assembled_pre.getBondType(bond_b)) {
                                System.out.println("error while assembling fragments -> different bond types?");
                            }
                            assembled_pre.addBond(conn_atom_a, conn_atom_b, assembled_pre.getBondType(bond_a));

                            //m_parts.get(0).ensureHelperArrays(StereoMolecule.cHelperCIP); // new
                            assembled_pre.ensureHelperArrays(StereoMolecule.cHelperNeighbours); // new

                            // remove Us
                            assembled_pre.markAtomForDeletion(zi);
                            assembled_pre.markAtomForDeletion(zj);

                            assembled_pre.deleteMarkedAtomsAndBonds();
                            assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);
                            if(logLevel_assemble>0) {
                                System.out.println("After connecting:");
                                System.out.println( HyperspaceUtils.idcodeToSmiles(assembled_pre.getIDCode()) );
                            }

                            //if(logLevel>0){System.out.println( "Assembly status: " + Hyperspace.stereoMoleculeToSmiles(m_parts.get(0)) );}

                            still_connecting = true;
                            created_connection = true;
                            break;

                        }
                        if(created_connection) {
                            break;
                        }
                    }
                }
//                if(!still_connecting) {
//                    break;
//                }
            }
        }
        assembled_pre.deleteMarkedAtomsAndBonds();
        assembled_pre.ensureHelperArrays(StereoMolecule.cHelperCIP);

        if(logLevel_assemble>0) {
            System.out.println("Final assembly:");
            System.out.println( idcodeToSmiles(assembled_pre.getIDCode() ) );
        }

        // ! actually it is also ok if we just remove the MATCHED uraniums.. (then we see for incomplete matches where the additional bbs go.
        StereoMolecule assembled = assembled_pre;
        if(false) {
            for (int zi = 0; zi < assembled.getAtoms(); zi++) {
                if (assembled.getAtomicNo(zi) == 92) {
                    assembled.markAtomForDeletion(zi);
                }
            }
            assembled.deleteMarkedAtomsAndBonds();
        }
        assembled.ensureHelperArrays(StereoMolecule.cHelperCIP);

        return assembled;
    }

    public static String idcodeToSmiles(String idcode) {
        StereoMolecule m = new StereoMolecule();
        IDCodeParser p = new IDCodeParser();
        p.parse(m,idcode);
        IsomericSmilesCreator smiles_creator = new IsomericSmilesCreator(m);
        return smiles_creator.getSmiles();
    }





    public static class ExpandedCombinatorialHit {
        public final String assembled_idcode;
        public final List<SynthonSpace.FragId> fragments;
        public ExpandedCombinatorialHit(String assembled_idcode, List<SynthonSpace.FragId> frags) {
            this.assembled_idcode = assembled_idcode;
            this.fragments = frags;
        }
    }



    /**
     *
     * @param hi
     * @return first element is idcode of assembled hit. in the second element, the first
     */
    public static List<ExpandedCombinatorialHit> expandCombinatorialHit(SynthonSpace.CombinatorialHit hi) {
        return expandCombinatorialHit(hi,1024);
    }

    /**
     *
     * @param hi
     * @return
     */
    public static List<ExpandedCombinatorialHit> expandCombinatorialHit(SynthonSpace.CombinatorialHit hi, int max_expanded) {
        List<List<SynthonSpace.FragId>> frag_sets =  new ArrayList<>( hi.hit_fragments.values() );

        List<List<SynthonSpace.FragId>> current_assemblies = new ArrayList<>();
        current_assemblies.add(new ArrayList<>());

        for(int zi=0;zi<frag_sets.size();zi++) {
            if(Thread.currentThread().isInterrupted()) {return new ArrayList<>();}
            List<List<SynthonSpace.FragId>> next_assemblies = new ArrayList<>();
            for( List<SynthonSpace.FragId> ca : current_assemblies) {
                if(next_assemblies.size()>=max_expanded) {break;}
                for(SynthonSpace.FragId fin : frag_sets.get(zi)) {
                    if(next_assemblies.size()>=max_expanded) {break;}
                    List<SynthonSpace.FragId> new_list = new ArrayList<>(ca);
                    new_list.add(fin);
                    next_assemblies.add(new_list);
                }
            }
            current_assemblies = next_assemblies;
        }

        // and then assemble..
        //List<Pair<StereoMolecule,Pair<StereoMolecule[],String[]>>> results = new ArrayList<>();
        List<ExpandedCombinatorialHit> results = new ArrayList<>();
        IDCodeParser icp = new IDCodeParser();
        for(List<SynthonSpace.FragId> api : current_assemblies) {
            if(Thread.currentThread().isInterrupted()) {return results;}

            List<StereoMolecule> parts     = new ArrayList<>();
            List<String>         parts_ids = new ArrayList<>();
            for( SynthonSpace.FragId fid : api ) {
                parts_ids.add(fid.fragment_id);
                StereoMolecule mpi = new StereoMolecule();
                icp.parse(mpi,fid.idcode);
                parts.add(mpi);
            }

            StereoMolecule parts_as_array[] = parts.toArray(new StereoMolecule[parts.size()]);
            StereoMolecule assembly_i = SynthonAssembler.assembleSynthons(parts);

            //results.add( Pair.of( assembly_i , parts_as_array ) );
            results.add( new ExpandedCombinatorialHit( assembly_i.getIDCode() , api ) );
            //System.out.println("AssembledMol: Smiles: " + HyperspaceUtils.idcodeToSmiles(assembly_i.getIDCode()) );
        }

        return results;
    }

    public static ExpandedCombinatorialHit expandAssembly(List<SynthonSpace.FragId> api) {
        List<StereoMolecule> parts     = new ArrayList<>();
        List<String>         parts_ids = new ArrayList<>();
        IDCodeParser icp = new IDCodeParser();
        for( SynthonSpace.FragId fid : api ) {
            parts_ids.add(fid.fragment_id);
            StereoMolecule mpi = new StereoMolecule();
            icp.parse(mpi,fid.idcode);
            parts.add(mpi);
        }

        StereoMolecule parts_as_array[] = parts.toArray(new StereoMolecule[parts.size()]);
        StereoMolecule assembly_i = SynthonAssembler.assembleSynthons(parts);

        //results.add( Pair.of( assembly_i , parts_as_array ) );
        return new ExpandedCombinatorialHit( assembly_i.getIDCode() , api );
        //System.out.println("AssembledMol: Smiles: " + HyperspaceUtils.idcodeToSmiles(assembly_i.getIDCode()) );
    }

    public static void expandCombinatorialHit(SynthonSpace.CombinatorialHit hi, Consumer<ExpandedCombinatorialHit> hitConsumer, long max_expanded) {
        List<List<SynthonSpace.FragId>> frag_sets =  new ArrayList<>( hi.hit_fragments.values() );

        long cnt = 0;

        if(frag_sets.size()==2) {
            for(SynthonSpace.FragId fa : frag_sets.get(0)) {
                for(SynthonSpace.FragId fb : frag_sets.get(1)) {
                    List<SynthonSpace.FragId> frags = new ArrayList<>();
                    frags.add(fa); frags.add(fb);
                    hitConsumer.accept(expandAssembly(frags));
                    if( cnt++ >= max_expanded) {break;}
                }
            }
        }
        else if(frag_sets.size()==3) {
            for(SynthonSpace.FragId fa : frag_sets.get(0)) {
                for(SynthonSpace.FragId fb : frag_sets.get(1)) {
                    for(SynthonSpace.FragId fc : frag_sets.get(2)) {
                        List<SynthonSpace.FragId> frags = new ArrayList<>();
                        frags.add(fa);
                        frags.add(fb);
                        frags.add(fc);
                        hitConsumer.accept(expandAssembly(frags));
                        if( cnt++ >= max_expanded) {break;}
                    }
                }
            }
        }
        else {
            throw new RuntimeException("Not supported");
        }

    }






    public static void assembly_test_01() {


        System.out.print("\nLoad Enamine REAL Space SynthonSpace: ");


        SynthonSpace space = new SynthonSpace();

        //String input_file = "C:\\dev\\git2\\openchemlib\\enamine_pfp_space.data";
        //String input_file = "C:\\dev\\git2\\openchemlib\\enamine_pfp_space_2.data";
        //String input_file = "C:\\dev\\git2\\openchemlib\\enamine_pfp_test_01.data";
        //String input_file = "C:\\dev\\git2\\openchemlib\\enamine_pfp_test_02.data";
        String input_file = "C:\\dev\\git2\\openchemlib\\enamine_pfp_space_3.data";
        String data       = HyperspaceUtils.readFileIntoString(input_file);
        space.loadFromString(data);

        StereoMolecule a = new StereoMolecule();
        StereoMolecule b = new StereoMolecule();
        StereoMolecule c = new StereoMolecule();

        SmilesParser sp = new SmilesParser();

    }
}
