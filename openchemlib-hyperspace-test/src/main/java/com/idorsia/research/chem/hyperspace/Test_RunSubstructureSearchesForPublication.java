package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.DWARFileParser;
import org.apache.commons.io.input.CountingInputStream;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class Test_RunSubstructureSearchesForPublication {

    /**
     *
     * Has three parameter:
     * 1. the path to the hyperspace data file containing the combinatorial library to search in.
     * 2. the path to dwarfile containing the queries
     * 3. the column name of the column containing the queries in the dwarfile
     *
     * @param args
     */
    public static void main(String args[]) {

        String path_to_datafile = args[0];
        String path_to_query_dwar_file = args[1];
        String column_name = args[2];


        List<StereoMolecule> mols_a = loadMolecules_A(path_to_query_dwar_file,column_name);



        SynthonSpace space = null;

        try {
            //space = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
            // todo: update HyperspaceIOUtils.loadSynthonSpace function to include zip..
            //space = loadSpace("/home/liphath1/dev5/git/hyperspace/REAL_Space_latest_FragFp.data");
            //space = loadSpace("/home/liphath1/dev5/git/hyperspace/IVCS3_3S_FragFp.data");
            space = loadSpace(path_to_datafile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        runSearchesAndCreateReport(space,mols_a);
    }


    public static void runSearchesAndCreateReport(SynthonSpace space, List<StereoMolecule> mols_a) {
        int rand_i = (int) (Math.abs( Math.random() ) * 10000);
        BufferedWriter out = null;
        try {
            out = new BufferedWriter(new FileWriter("benchmark_run_out_"+rand_i+".csv"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        //List<String> header = new ArrayList<>();
        try {
            out.write("Structure[idcode],ffp_bits,chbb,ch1,ch2,ch3,chtot,thbb,th1,th2,th3,thtot,tbb,t1,t2,t3,ttot,tbbw,t1w,t2w,t3w,ttotw,splits_all,splits_valid,rxn_mappings,processed_labeled_decomps,total_sss,exacthit,t_exact,tw_exact\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

        // sort from small to large:
        mols_a.sort( (x,y) -> Integer.compare( x.getAtoms() , y.getAtoms() ) );

        for(int zi=0;zi<mols_a.size();zi++) {
            StereoMolecule mi = mols_a.get(zi);
            StereoMolecule mi_fixed = new StereoMolecule(mi);
            mi_fixed.setFragment(true);
            mi_fixed.ensureHelperArrays(Molecule.cHelperCIP);

            for(int za=0;za<mi_fixed.getAtoms();za++) {
                mi_fixed.setAtomQueryFeature(za,StereoMolecule.cAtomQFNoMoreNeighbours, true);
            }
            mi_fixed.ensureHelperArrays(Molecule.cHelperCIP);

            // run fragment search:
            SearchReport report_i = runSearch(space,mi,null);

            // run exact search:
            SearchReport report_i_exact = runSearch(space,mi_fixed,null);
            boolean exact_hit = (report_i_exact.thits_s1+report_i_exact.thits_s2+report_i_exact.thits_s3) > 0;
            long    t_exact = report_i_exact.t1+report_i_exact.t2+report_i_exact.t3;
            long    t_exact_w = report_i_exact.t1w+report_i_exact.t2w+report_i_exact.t3w;

            try {
                out.write(report_i.toCSV() + "," + exact_hit+","+t_exact+","+t_exact_w+"\n");
                out.flush();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        try {
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * NOTE: if cdp is zero, then a new cached dp is created for the run.
     *
     * @param space
     * @param mi_a
     * @param cdp
     */
    public static SearchReport runSearch(SynthonSpace space, StereoMolecule mi_a, CachedDescriptorProvider cdp) {

        System.out.println("Smiles : "+HyperspaceUtils.idcodeToSmiles(mi_a.getIDCode()));
        System.out.println("Atoms  : "+mi_a.getAtoms());

        List<SynthonSpace.CombinatorialHit> hits_1s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_2s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_3s = new ArrayList<>();
        List<SynthonSpace.CombinatorialHit> hits_4s = new ArrayList<>();

        Set<String> empty_set = new HashSet<>();

        boolean empty_cache = false;
        if(cdp==null) {
            cdp = new CachedDescriptorProvider(space.getDescriptorHandler().getInfo().shortName);
            empty_cache = true;
        }

        StereoMolecule mi = new StereoMolecule(mi_a);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        //fi.setFragment(true);

        long ts_bb_a = System.currentTimeMillis();
        List<SynthonSpace.ExpandedHit> hits_bbs = SynthonSpace.screen_building_blocks(space,cdp,mi,8);
        long ts_bb_b = System.currentTimeMillis();

        List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics_1 = Collections.synchronizedList(new ArrayList<>());
        //List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics_2 = Collections.synchronizedList(new ArrayList<>());
        //List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> output_statistics_3 = Collections.synchronizedList(new ArrayList<>());


        long ts_a = System.currentTimeMillis();
        hits_1s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 1, 2, empty_set , 1000, 8, output_statistics_1);

        long ts_b = System.currentTimeMillis();
        hits_2s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 2, 3, empty_set , 1000, 8, output_statistics_1);

        long ts_c = System.currentTimeMillis();
        hits_3s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 3, 3, empty_set , 1000, 8, output_statistics_1);

        long ts_d = System.currentTimeMillis();

        long tbb = ts_bb_b - ts_bb_a;
        long t1 = ts_b-ts_a;
        long t2 = ts_c-ts_b;
        long t3 = ts_d-ts_c;

        System.out.println("Hits (bb): " + hits_bbs.size());
        System.out.println("Hits (1s): " + hits_1s.size());
        System.out.println("Hits (2s): " + hits_2s.size());
        System.out.println("Hits (3s): " + hits_3s.size());


        // analzye output statistics:
        //System.out.println("Search complete!"+" [Splits="+num_splits+"]"+"Query= "+mol.getIDCode());
        System.out.println("Statistics:");
        int  splits_total = output_statistics_1.size();
        int  splits_valid = (int) output_statistics_1.stream().filter( si -> si.valid_split ).count();
        System.out.println("Splits: "+splits_total+" Valid: "+splits_valid);
        int total_possible_rxn_mappings    = output_statistics_1.stream().flatMapToInt(si -> si.possible_rxn_mappings.values().stream().mapToInt( xmi -> xmi.size() )).sum();
        int total_processed_labeled_splits = output_statistics_1.stream().mapToInt(si -> si.total_labeled_connector_splits_processed).sum();
        int total_enumerated_sss           = output_statistics_1.stream().mapToInt(si -> si.total_sss_performed).sum();
        System.out.println("Possible Rxn Maps: "+total_possible_rxn_mappings);
        System.out.println("Processed labeled splits: "+total_processed_labeled_splits);
        System.out.println("Synthon substructure searches: "+total_enumerated_sss);


        //

        long t_bb_warm = -1;
        long t1_warm = -1;
        long t2_warm = -1;
        long t3_warm = -1;

        if(true) {
            // now run it a second time with warm cache:
            List<SynthonSpace.CombinatorialHit> hits2_1s = new ArrayList<>();
            List<SynthonSpace.CombinatorialHit> hits2_2s = new ArrayList<>();
            List<SynthonSpace.CombinatorialHit> hits2_3s = new ArrayList<>();

            long ts_warm_bb_a = System.currentTimeMillis();
            List<SynthonSpace.ExpandedHit> hits2_bbs = SynthonSpace.screen_building_blocks(space,cdp,mi,8);
            long ts_warm_bb_b = System.currentTimeMillis();

            long ts_warm_a = System.currentTimeMillis();
            hits2_1s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 1, 2, empty_set, 1000, 8);
            long ts_warm_b = System.currentTimeMillis();
            hits2_2s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 2, 3, empty_set, 1000, 8);
            long ts_warm_c = System.currentTimeMillis();
            hits2_3s = space.findExpandedHits_withConnProximityMatching(space, cdp, mi, 3, 3, empty_set, 1000, 8);
            long ts_warm_d = System.currentTimeMillis();
            t_bb_warm = ts_warm_bb_b - ts_warm_bb_a;
            t1_warm = ts_warm_b-ts_warm_a;
            t2_warm = ts_warm_c-ts_warm_b;
            t3_warm = ts_warm_d-ts_warm_c;
        }

        //System.out.println("Hits (4s): " + hits_4s.size());

        // add cardinality of fingerprint:
        int bits_in_ffp = cdp.getFP_cached(mi).cardinality();

        if(hits_1s.size()>0 || hits_2s.size()>0 || hits_3s.size()>0 ) {
            System.out.println("HITS!!");

            // compute thits:
            long thits_bb = hits_bbs.stream().mapToLong( hi -> computeNumberOfHitsOfExpandedHit(space,hi) ).sum();
            long thits1 = hits_1s.stream().mapToLong( hi -> hi.hit_fragments.values().stream().mapToLong( li -> li.size() ).reduce((x,y)->x*y).getAsLong() ).sum();
            long thits2 = hits_2s.stream().mapToLong( hi -> hi.hit_fragments.values().stream().mapToLong( li -> li.size() ).reduce((x,y)->x*y).getAsLong() ).sum();
            long thits3 = hits_3s.stream().mapToLong( hi -> hi.hit_fragments.values().stream().mapToLong( li -> li.size() ).reduce((x,y)->x*y).getAsLong() ).sum();

            return new SearchReport(mi.getIDCode(),bits_in_ffp,hits_bbs.size(),hits_1s.size(),hits_2s.size(),hits_3s.size(),thits_bb,thits1,thits2,thits3,tbb,t1,t2,t3,t_bb_warm,t1_warm,t2_warm,t3_warm,
                    splits_total,splits_valid,total_possible_rxn_mappings,total_processed_labeled_splits,total_enumerated_sss);
        }
        else {
            System.out.println("NO HITS!!");
            return new SearchReport(mi.getIDCode(),bits_in_ffp,0,0,0,0,0,0,0,0,tbb,t1,t2,t3,t_bb_warm,t1_warm,t2_warm,t3_warm,
                    splits_total,splits_valid,total_possible_rxn_mappings,total_processed_labeled_splits,total_enumerated_sss);
        }

    }

    public static long computeNumberOfHitsOfExpandedHit(SynthonSpace space, SynthonSpace.ExpandedHit hi) {
        // determine hits in synthon rxn:
        Map<Integer,Integer> map_frag_to_num_fragments = new HashMap<>();
        space.getFragTypes(hi.rxn).values().stream().forEach( fi ->map_frag_to_num_fragments.put( fi.frag, space.getSynthonSet(fi.rxn_id,fi.frag).size() ) );
        long exp_hits = 1;
        for(Integer fi : map_frag_to_num_fragments.keySet()) {
            if(fi==hi.hit_fragments.values().iterator().next().frag) {
                exp_hits = exp_hits * 1; // dont do anything
            }
            else {
                exp_hits = exp_hits * map_frag_to_num_fragments.get(fi);
            }
        }
        return exp_hits;
    }

    public static class SearchReport {
        public final String idc;
        public final int bits_in_ffp;

        public final int chits_bb;
        public final int chits_s1;
        public final int chits_s2;
        public final int chits_s3;

        public final long thits_bb;
        public final long thits_s1;
        public final long thits_s2;
        public final long thits_s3;

        public final long tbb;
        public final long t1;
        public final long t2;
        public final long t3;

        public final long tbbw;
        public final long t1w;
        public final long t2w;
        public final long t3w;

        public final int splits_total;
        public final int splits_valid;
        public final int rxn_mappings;
        public final int processed_labeled_splits;
        public final int performed_sss;



        public SearchReport(String idc, int bits_in_ffp,
                            int chits_bb, int chits_s1, int chits_s2, int chits_s3,
                            long thits_bb, long thits_s1, long thits_s2, long thits_s3,
                            long tbb, long t1, long t2, long t3,
                            long tbbw, long t1w, long t2w, long t3w,
                            int splits_total, int splits_valid, int rxn_mappings, int processed_labeled_splits, int performed_sss) {
            this.idc = idc;
            this.bits_in_ffp = bits_in_ffp;
            this.chits_bb = chits_bb;
            this.chits_s1 = chits_s1;
            this.chits_s2 = chits_s2;
            this.chits_s3 = chits_s3;
            this.thits_bb = thits_bb;
            this.thits_s1 = thits_s1;
            this.thits_s2 = thits_s2;
            this.thits_s3 = thits_s3;
            this.tbb = tbb;
            this.t1 = t1;
            this.t2 = t2;
            this.t3 = t3;
            this.tbbw = tbbw;
            this.t1w = t1w;
            this.t2w = t2w;
            this.t3w = t3w;
            this.splits_total = splits_total;
            this.splits_valid = splits_valid;
            this.rxn_mappings = rxn_mappings;
            this.processed_labeled_splits = processed_labeled_splits;
            this.performed_sss = performed_sss;
        }

        public String toCSV() {
            List<String> parts = new ArrayList<>();
            parts.add( idc);
            parts.add(""+bits_in_ffp);
            parts.add(""+chits_bb);
            parts.add(""+chits_s1);
            parts.add(""+chits_s2);
            parts.add(""+chits_s3);
            parts.add(""+(chits_bb+chits_s1+chits_s2+chits_s3));
            parts.add(""+thits_bb);
            parts.add(""+thits_s1);
            parts.add(""+thits_s2);
            parts.add(""+thits_s3);
            parts.add(""+(thits_bb+thits_s1+thits_s2+thits_s3));

            parts.add(""+tbb);
            parts.add(""+t1);
            parts.add(""+t2);
            parts.add(""+t3);
            parts.add(""+(t1+t2+t3));

            parts.add(""+tbbw);
            parts.add(""+t1w);
            parts.add(""+t2w);
            parts.add(""+t3w);
            parts.add(""+(t1w+t2w+t3w));

            parts.add(""+splits_total);
            parts.add(""+splits_valid);
            parts.add(""+rxn_mappings);
            parts.add(""+processed_labeled_splits);
            parts.add(""+performed_sss);


            return String.join(",",parts);
        }
    }




    public static List<StereoMolecule> loadMolecules_A(String dwarfile, String column_name) {

        //DWARFileParser in = new DWARFileParser("/home/liphath1/datasets/ListOfSmallMoleculeDrugsFrlomDrugCentral2021.dwar");
        //int si_mol = in.getSpecialFieldIndex("Structure of SMILES");
        DWARFileParser in = new DWARFileParser(dwarfile);
        int si_mol        = in.getSpecialFieldIndex(column_name);

        List<StereoMolecule> mols = new ArrayList<>();
        while(in.next()) {
            StereoMolecule mi = HyperspaceUtils.parseIDCode(in.getSpecialFieldData(si_mol));

            if(mi.getAtoms()<=32 && mi.getAtoms()>=12) {
                mols.add(mi);
            }
        }

        System.out.println("total molecules: "+mols.size());
        return mols;
    }



    public static SynthonSpace loadSpace(String file) throws IOException{
        File fi_file = new File(file);
        long file_size_in_bytes = fi_file.length();

        FileInputStream file_in = new FileInputStream(fi_file);
        //counting_in_stream = new CountingInputStream(new BufferedInputStream(file_in));
        //counting_in_stream = new CountingInputStream( new GZIPInputStream( new BufferedInputStream(file_in)));
        CountingInputStream counting_in_stream = new CountingInputStream( new BufferedInputStream(file_in));

        ObjectInputStream in = new ObjectInputStream( new GZIPInputStream( counting_in_stream ) );

        SynthonSpace space_a = null;
        try {
            space_a = (SynthonSpace) in.readObject();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        space_a.initAfterJavaDeserialization();
        return space_a;
    }

}
