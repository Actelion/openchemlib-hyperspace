package com.idorsia.research.chem.hyperspace.gui;

import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.SynthonSpace;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.*;
import java.util.List;

import static com.idorsia.research.chem.hyperspace.SynthonAssembler.expandCombinatorialHit;

public class JHyperspaceSubstructureSearchResultsTable extends JPanel {

    public static class SubstructureSearchResultsTableModel_A  extends AbstractTableModel
    {
        private int max_synthons = 3;

        private Map<String,String> fragid_to_bbid = new HashMap<>();

        public SubstructureSearchResultsTableModel_A(Map<String,String> fragid_to_bbid, int max_synthons) {
            this.fragid_to_bbid = fragid_to_bbid;
            this.max_synthons = max_synthons;
        }

        List<SynthonSpace.CombinatorialHit> combi_hits = new ArrayList<>();

        public Map<String,String> getFragId2BBId() {
            return this.fragid_to_bbid;
        }

        public SynthonSpace.CombinatorialHit getHit(int idx) {
            return this.combi_hits.get(idx);
        }

        public void setResults(List<SynthonSpace.CombinatorialHit> combi_hits) {
            this.combi_hits = new ArrayList<>(combi_hits);
            fireTableDataChanged();
        }

        @Override
        public int getRowCount() {
            return combi_hits.size();
        }

        @Override
        public int getColumnCount() {
            return 2+2*max_synthons;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SynthonSpace.CombinatorialHit hit = combi_hits.get(rowIndex);
            if(columnIndex==0) {
                return hit.rxn;
            }
            if(columnIndex==1){
                int results = hit.hit_fragments.values().stream().mapToInt( fi -> fi.size() ).reduce( (x,y)->x*y ).getAsInt();
                return results;
            }
            if(columnIndex>=2 && columnIndex<2+max_synthons){
                if(hit.sri==null) {
                    return "";
                }
                if( hit.sri.fragments.length > (columnIndex-2) ) {
                    return hit.sri.fragments[columnIndex-2].getIDCode();
                }
                return "";
            }
            // TODO: this is actually not correct. we somehow have to add information to the combinatorial hit, which
            //       split result fragment belongs to which frag type.
            if(columnIndex>=2+max_synthons && columnIndex<2+2*max_synthons){
                int frag_idx = (columnIndex-2)-max_synthons;
                if( hit.sri.fragments.length > frag_idx ) {
                    //List<Integer> num_hits = hit.hit_fragments.values().stream().map(fi -> fi.size()).collect(Collectors.toList());
                    //return ""+num_hits.get(frag_idx);
                    int num_hits = hit.hit_fragments.get( hit.mapping.get(frag_idx).getLeft() ).size();
                    return ""+num_hits;
                }
                return "";
            }
//            if(columnIndex > 2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    List<SynthonSpace.FragType> fti = new ArrayList<>(hit.hit_fragments.keySet());
//                    fti.sort((x,y)->x.);
//                }
//            }

            return null;
        }
    }

    public static class SubstructureSearchResultsTableModel_B  extends AbstractTableModel {
        private int max_synthons = 3;

        private Map<String,String> fragid_to_bbid = new HashMap<>();


        public SubstructureSearchResultsTableModel_B(Map<String,String> fragid_to_bbid, int max_synthons) {
            this.fragid_to_bbid = fragid_to_bbid;
            this.max_synthons = max_synthons;
        }

        SynthonSpace.CombinatorialHit hit = null;
        List<SynthonAssembler.ExpandedCombinatorialHit> expanded_hit = new ArrayList<>();

        public void setExpandedHit(List<SynthonAssembler.ExpandedCombinatorialHit> expanded_hit) {
            this.expanded_hit = expanded_hit;
            this.fireTableDataChanged();
        }

        public void setHit(SynthonSpace.CombinatorialHit combi_hit) {

            //List<Pair<StereoMolecule, StereoMolecule[]>> expanded = new ArrayList<>();
            Thread compute_expansion = new Thread() {
                @Override
                public void run() {
                    List<SynthonAssembler.ExpandedCombinatorialHit> expanded = expandCombinatorialHit(combi_hit);
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            setExpandedHit(expanded);
                        }
                    });
                }
            };
            compute_expansion.start();
        }

        @Override
        public int getRowCount() {
            return expanded_hit.size();
        }

        @Override
        public int getColumnCount() {
            return 1 + 2 * max_synthons;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            SynthonAssembler.ExpandedCombinatorialHit hit = this.expanded_hit.get(rowIndex);
            if (columnIndex == 0) {
                return hit.assembled_idcode;
            }
            if (columnIndex >= 1 && columnIndex < 1 + max_synthons) {
                int frag_idx = columnIndex - 1;
                if (hit.fragments.size() > frag_idx) {
                    return hit.fragments.get(frag_idx).idcode;
                } else {
                    return "";
                }
            }
            if (columnIndex >= 1 + max_synthons && columnIndex < 1 + 2 * max_synthons) {
                int frag_idx = (columnIndex - 1)-max_synthons;
                if (hit.fragments.size() > frag_idx) {
                    //String bb_id = fragid_to_bbid.get( hit.fragments.get(frag_idx).fragment_id );
                    String bb_id = hit.fragments.get(frag_idx).fragment_id;
                    if(bb_id != null) { return bb_id; }
                    return "<ERROR>";
                } else {
                    return "N/A";
                }
            }
//            if(columnIndex==1){
//                int results = hit.hit_fragments.values().stream().mapToInt( fi -> fi.size() ).reduce( (x,y)->x*y ).getAsInt();
//                return results;
//            }
//            if(columnIndex>=2 && columnIndex<=2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    hit.sri.fragments[columnIndex-2].getIDCode();
//                }
//            }
//            if(columnIndex > 2+max_synthons){
//                if( hit.sri.fragments.length > (columnIndex-2) ) {
//                    List<SynthonSpace.FragType> fti = new ArrayList<>(hit.hit_fragments.keySet());
//                    fti.sort((x,y)->x.);
//                }
//            }
            return null;
        }
    }


    JSplitPane sp_a = new JSplitPane();
    JPanel jp_left;
    JPanel jp_right;
    JPanel jp_left_top;
    JPanel jp_right_top;
    JSlider js_left_size;
    JSlider js_right_size;
    JScrollPane jscp_left  = new JScrollPane();
    JScrollPane jscp_right = new JScrollPane();
    JTable jt_left = new JTable();
    JTable jt_right = new JTable();

    SubstructureSearchResultsTableModel_A tablemodel_a = null;
    SubstructureSearchResultsTableModel_B tablemodel_b = null;

    public JHyperspaceSubstructureSearchResultsTable(SubstructureSearchResultsTableModel_A model) {
        this.tablemodel_a = model;
        this.tablemodel_b = new SubstructureSearchResultsTableModel_B(this.tablemodel_a.getFragId2BBId(),3);
        initGUI();
        initTableGraphics();
        initListeners();
    }

    private void initGUI() {
        this.removeAll();
        this.setLayout(new BorderLayout());

        this.sp_a = new JSplitPane();
        jt_left  = new JTable(this.tablemodel_a);
        jt_right = new JTable(this.tablemodel_b);

        jp_left = new JPanel(); jp_left.setLayout(new BorderLayout());
        jp_right = new JPanel(); jp_right.setLayout(new BorderLayout());

        jp_left_top = new JPanel();
        jp_right_top = new JPanel();

        jp_left_top.setLayout(new FlowLayout());
        jp_right_top.setLayout(new FlowLayout());

        js_left_size = new JSlider(8,128);
        js_right_size = new JSlider(8,128);

        jp_left_top.add(js_left_size);
        jp_right_top.add(js_right_size);

        jp_left.add(jp_left_top,BorderLayout.NORTH);
        jp_right.add(jp_right_top,BorderLayout.NORTH);


        jscp_left  = new JScrollPane(jt_left);
        jscp_right = new JScrollPane(jt_right);

        jp_left.add(jscp_left,BorderLayout.CENTER);
        jp_right.add(jscp_right,BorderLayout.CENTER);

        this.sp_a.setLeftComponent(jp_left);
        this.sp_a.setRightComponent(jp_right);

        this.add(this.sp_a,BorderLayout.CENTER);
    }

    private void initTableGraphics() {
        // TODO: use something different..

//        this.jt_left.getColumnModel().getColumn(2).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_left.getColumnModel().getColumn(3).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_left.getColumnModel().getColumn(4).setCellRenderer(new SubstanceFragmentRenderer());
//
//        this.jt_right.getColumnModel().getColumn(0).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_right.getColumnModel().getColumn(1).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_right.getColumnModel().getColumn(2).setCellRenderer(new SubstanceFragmentRenderer());
//        this.jt_right.getColumnModel().getColumn(3).setCellRenderer(new SubstanceFragmentRenderer());
    }

    private void initListeners() {

        this.js_left_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_left.setRowHeight(js_left_size.getValue());
            }
        });
        this.js_right_size.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                jt_right.setRowHeight(js_right_size.getValue());
            }
        });

        this.jt_left.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        this.jt_left.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                int model_idx = jt_left.convertRowIndexToModel( jt_left.getSelectedRow() );
                if(model_idx<0) {
                    System.out.println("invalid hit selected -> skip"); return;
                }
                SynthonSpace.CombinatorialHit hit = tablemodel_a.getHit(model_idx);
                tablemodel_b.setHit(hit);
            }
        });

    }

    public void setHits(List<SynthonSpace.CombinatorialHit> hits) {
        this.tablemodel_a.setResults(hits);
    }


//    public static class ExpandedCombinatorialHit {
//        public final String assembled_idcode;
//        public final List<SynthonSpace.FragId> fragments;
//        public ExpandedCombinatorialHit(String assembled_idcode, List<SynthonSpace.FragId> frags) {
//            this.assembled_idcode = assembled_idcode;
//            this.fragments = frags;
//        }
//    }
//
//    /**
//     *
//     * @param hi
//     * @return first element is idcode of assembled hit. in the second element, the first
//     */
//    public static List<ExpandedCombinatorialHit> expandCombinatorialHit(SynthonSpace.CombinatorialHit hi) {
//        List<List<SynthonSpace.FragId>> frag_sets =  new ArrayList<>( hi.hit_fragments.values() );
//
//        List<List<SynthonSpace.FragId>> current_assemblies = new ArrayList<>();
//        current_assemblies.add(new ArrayList<>());
//
//        for(int zi=0;zi<frag_sets.size();zi++) {
//            List<List<SynthonSpace.FragId>> next_assemblies = new ArrayList<>();
//            for( List<SynthonSpace.FragId> ca : current_assemblies) {
//                for(SynthonSpace.FragId fin : frag_sets.get(zi)) {
//                    List<SynthonSpace.FragId> new_list = new ArrayList<>(ca);
//                    new_list.add(fin);
//                    next_assemblies.add(new_list);
//                }
//            }
//            current_assemblies = next_assemblies;
//        }
//
//        // and then assemble..
//        //List<Pair<StereoMolecule,Pair<StereoMolecule[],String[]>>> results = new ArrayList<>();
//        List<ExpandedCombinatorialHit> results = new ArrayList<>();
//        IDCodeParser icp = new IDCodeParser();
//        for(List<SynthonSpace.FragId> api : current_assemblies) {
//            List<StereoMolecule> parts     = new ArrayList<>();
//            List<String>         parts_ids = new ArrayList<>();
//            for( SynthonSpace.FragId fid : api ) {
//                parts_ids.add(fid.fragment_id);
//                StereoMolecule mpi = new StereoMolecule();
//                icp.parse(mpi,fid.idcode);
//                parts.add(mpi);
//            }
//
//            StereoMolecule parts_as_array[] = parts.toArray(new StereoMolecule[parts.size()]);
//            StereoMolecule assembly_i = SynthonAssembler.assembleSynthons(parts);
//
//            //results.add( Pair.of( assembly_i , parts_as_array ) );
//            results.add( new ExpandedCombinatorialHit( assembly_i.getIDCode() , api ) );
//            System.out.println("AssembledMol: Smiles: " + HyperspaceUtils.idcodeToSmiles(assembly_i.getIDCode()) );
//        }
//
//        return results;
//    }

}
