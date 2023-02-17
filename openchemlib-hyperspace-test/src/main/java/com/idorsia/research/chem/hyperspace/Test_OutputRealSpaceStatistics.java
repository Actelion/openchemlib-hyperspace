package com.idorsia.research.chem.hyperspace;

import java.io.IOException;

public class Test_OutputRealSpaceStatistics {

    /**
     *
     * Has one parameter: the path to the hyperspace data file containing
     * the combinatorial library to search in.
     *
     * @param args
     */
    public static void main(String args[]) {
        String path_to_datafile = args[0];

        SynthonSpace space = null;

        try {
            //space = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
            // todo: update HyperspaceIOUtils.loadSynthonSpace function to include zip..
            //space = Test_RunSubstructureSearchesForPublication.loadSpace("/home/liphath1/dev5/git/hyperspace/REAL_Space_latest_FragFp.data");
            space = Test_RunSubstructureSearchesForPublication.loadSpace( path_to_datafile );
        } catch (IOException e) {
            e.printStackTrace();
        }

        for( String ccs : space.rxns_by_connector_config.keySet()) {
            System.out.println(ccs + " -> "+space.rxns_by_connector_config.get(ccs).size());
        }

    }

}
