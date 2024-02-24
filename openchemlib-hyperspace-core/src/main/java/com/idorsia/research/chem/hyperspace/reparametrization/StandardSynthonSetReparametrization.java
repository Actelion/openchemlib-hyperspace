package com.idorsia.research.chem.hyperspace.reparametrization;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class StandardSynthonSetReparametrization implements SynthonSetReparametrization {

    private boolean completeRingSystem = true;
    private int distance = 4;

    public StandardSynthonSetReparametrization() {

    }

    public StandardSynthonSetReparametrization(boolean completeRingSystem, int distance) {
        this.completeRingSystem = completeRingSystem;
        this.distance = distance;
    }

    @Override
    public StereoMolecule processSynthon(StereoMolecule mi) {
        StereoMolecule mi2 = new StereoMolecule(mi);
        mi2.ensureHelperArrays(Molecule.cHelperCIP);
        return process_synthon(mi2,this.distance,this.completeRingSystem);
    }



    /**
     * This method performs the crucial step: starting at the connectors, it
     * considers atoms up to distance away from the connector. Then, it visits
     * all ring systems that are already partially visited.
     * Then, the method cuts off all rests, and adds for the remaining parts
     * exit vector atoms (atomic number 90, Thorium)
     *
     *
     * @param distance
     * @return
     */
    public static StereoMolecule process_synthon(StereoMolecule mi, int distance, boolean complete_ring_systems) {
        mi.ensureHelperArrays(Molecule.cHelperCIP);

        // visit nodes from every connector atom
        boolean visited[] = new boolean[mi.getAtoms()];
        for(int zi=0;zi<mi.getAtoms();zi++) {
            if( mi.getAtomicNo(zi)==92 || mi.getAtomicNo(zi)==93 || mi.getAtomicNo(zi)==94 || mi.getAtomicNo(zi)==95 ) {
                boolean visited_i []= new boolean[mi.getAtoms()];
                ArrayList<String> features = new ArrayList<>();
                visit_df_recursive( mi,distance,-1,zi,"<"+mi.getAtomicNo(zi)+":"+zi+">",features,visited_i);
                // add to "visited":
                for(int zj=0;zj<visited.length;zj++) {
                    visited[zj] |= visited_i[zj];
                }
            }
        }

        // complete ring systems:
        if(complete_ring_systems) {
            for (int zi = 0; zi < mi.getAtoms(); zi++) {
                if (visited[zi]) {
                    complete_ring_systems_recursive(mi, zi, visited);
                }
            }
        }

        // add exit vectors:
        Map<Integer,Integer> exit_vectors = new HashMap<>();
        //compute_exit_vectors(mi,exit_vectors,visited);
        compute_exit_vectors_2(mi,exit_vectors,visited);

        // cut :)
        for(int zi=0;zi<mi.getAtoms();zi++) {
            if(!visited[zi]) {
                mi.markAtomForDeletion(zi);
            }
        }
        int old_to_new_pos[] = mi.deleteMarkedAtomsAndBonds();
        if(old_to_new_pos==null){old_to_new_pos = new int[mi.getAtoms()];for(int zi=0;zi<mi.getAtoms();zi++){ old_to_new_pos[zi]=zi; }}

        // update exit vector positions
        Map<Integer,Integer> exit_vectors_new = new HashMap<>();
        for( int poi :  exit_vectors.keySet()) {
            exit_vectors_new.put( old_to_new_pos[poi] , exit_vectors.get(poi) );
        }

        //add exit vectors:
        for(int pi : exit_vectors_new.keySet()) {
            for(int zi=0;zi<exit_vectors_new.get(pi);zi++){
                //int pni = mi.addAtom( 92);
                int pni = mi.addAtom( 90); // !! DONT take 92,93,94,95 .. because we will have to assemble these..
                mi.addBond( pi,pni,Molecule.cBondTypeSingle);
            }
        }
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        return mi;
    }

    public static void visit_df_recursive(StereoMolecule m, int remaining_steps, int last_pos, int pos, String current_path, List<String> features, boolean visited[]) {
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
            //boolean visited_next[] = Arrays.copyOf(visited,visited.length);
            boolean visited_next[] = visited;
            int next_pos = m.getConnAtom(pos,zi);
            if(visited[next_pos]){continue;}
            visit_df_recursive(m,remaining_steps-1,pos,next_pos,next_path,features,visited_next);
        }
    }

    public static void complete_ring_systems_recursive( StereoMolecule m, int pos, boolean visited[]) {
        if(m.isRingAtom(pos)) {
            for(int zi =0;zi<m.getConnAtoms(pos);zi++) {
                int ni = m.getConnAtom(pos,zi);
                if(visited[ni]){

                }
                else{
                    if(m.isRingAtom(ni)) {
                        visited[ni] = true;
                        complete_ring_systems_recursive( m,  pos, visited);
                    }
                }
            }
        }
    }

    public static void compute_exit_vectors( StereoMolecule m, Map<Integer,Integer> exit_vectors, boolean visited[] ) {
        for(int pos=0;pos<m.getAtoms();pos++) {
            if(visited[pos]) {
                int nonvisited_neighbors = 0;
                for(int zi =0;zi<m.getConnAtoms(pos);zi++) {
                    int ni = m.getConnAtom(pos, zi);
                    if(!visited[ni]) {
                        nonvisited_neighbors++;
                    }
                }
                if(nonvisited_neighbors>0) {
                    exit_vectors.put(pos, nonvisited_neighbors);
                }
            }
        }
    }

    /**
     * Requires that the atom that is identified as exit vector is not
     * only a single heavy atom, but has at least two heavy atom
     * neighbors
     *
     * @param m
     * @param exit_vectors
     * @param visited
     */
    public static void compute_exit_vectors_2( StereoMolecule m, Map<Integer,Integer> exit_vectors, boolean visited[] ) {
        for(int pos=0;pos<m.getAtoms();pos++) {
            if(visited[pos]) {
                int nonvisited_neighbors = 0;
                for(int zi =0;zi<m.getConnAtoms(pos);zi++) {
                    int ni = m.getConnAtom(pos, zi);
                    if(!visited[ni]) {
                        if(m.getConnAtoms(ni)>=2) {
                            nonvisited_neighbors++;
                        }
                    }
                }
                if(nonvisited_neighbors>0) {
                    exit_vectors.put(pos, nonvisited_neighbors);
                }
            }
        }
    }

}
