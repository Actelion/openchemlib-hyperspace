package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.StereoMolecule;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Test_RunMoreSubstructureSearchesForPublication_Scaffolds {


    /**
     *
     * Has one parameter: the path to the hyperspace data file containing
     * the combinatorial library to search in.
     *
     * @param args
     */
    public static void main(String args[]) {
        String path_to_datafile = args[0];

        //String idc_a = "ffcqB@BR\\HD`drnwJzrmFV\\Vl@DuUT@A@@P"; // O=C(c(nc1)ccc1F)N1CCN(Cc2ccccc2)CC1
        List<StereoMolecule> molecules = new ArrayList<>();

        if(true) {
            String idc_0 = "fakpRDEV[HLZ`BKSdbbbRfJTTTQbppVb}hIAb`ZAb@@LO@Hs`iJ|FNRu[t";
            String idc_a = "fakpRDEV[HLZ`BKSdbbbRfJTTTQbppVb}hIAb`ZAb@@H";
            String idc_b = "fj}Pb@EYdd@QdbbbRfNbbTVFBuDmAHLPPLP@A`hCIJ\\FJ\\e@";
            String idc_c = "f`qQbDEVd{p@btyHhhdibmLLEj`dFJZ@@pd@cNBdiw`";
            String idc_d = "fakpPDHNZAEjpeInUmvvYYxIVb}ZjjjZjijh@MoXIQgARUkpXyKvh\\NHzcfNxxSiNT{s`nR{Kfnzx[in@";
            String idc_e = "fakpPDHNZAEjpeInUmvvYYxIVb}ZjjjZjijh@LO@PRgIJlzNRu[vwlDhs`iJuxL\\e{TNGD]QsG\\\\ItgJ]ypWI]esW]\\Mtw@";

            String idc_f = "fdeaB@ER@DYHhhdichdlLEhLZBPX`@`@B";
            String idc_g = "fdeaB@ER@DYHhhdichdlLEhLZBPX`@`@CqFk@";
            //String idc_f = "";
            //String idc_g = "";
            molecules.add(HyperspaceUtils.parseIDCode(idc_0));
            molecules.add(HyperspaceUtils.parseIDCode(idc_a));
            molecules.add(HyperspaceUtils.parseIDCode(idc_b));
            molecules.add(HyperspaceUtils.parseIDCode(idc_c));
            //molecules.add(HyperspaceUtils.parseIDCode(idc_d));
            //molecules.add(HyperspaceUtils.parseIDCode(idc_e));
            molecules.add(HyperspaceUtils.parseIDCode(idc_f));
            molecules.add(HyperspaceUtils.parseIDCode(idc_g));
        }

        if(false) {
            String idc_0 = "eo`RJH@LBMJaEOBf`@aJYuf{Y]UYylJHljYifB`bjjjh@@`"; // CC(C)N(CCNC(c1csc(NC(c(c(O)c2)cc(OC)c2OC)=O)n1)=O)C(C)C
            String idc_a = "dg}H`LK`BDig{e\\h{ifX@H@H";
            String idc_b = "fdypb@LDtIP@aJYuf{WqPLZYifB@`@H";
            String idc_c = "flePR@LDykpPABTs}sNjxiYuLsA@UT@AyB`";
            String idc_d = "fbeHb@LDUXT`AFRJKIJEIQQSEJYiZjjj@B";
            String idc_e = "fbeHb@LDUXT`AFRJKIJEIQQSEJYiZjjj@CrC@";

            String idc_f = "dg}H`LK`BDig{e\\h{ifX@H@OBjSl";
            String idc_g = "dg}H`LK`BDig{e\\h{ifX@H@LNALKewaUIv";

            String idc_p3 = "fn}qP@DIZ@UoQrJIQJIYIIQGIQfhJBJijh@C[ANZjaceN@";
            String idc_p4 = "facqP@DEz@UoYrJIQJIYIIKPiJLuAPQUTuP@FqB\\tznFNTzs`";
            String idc_p5 = "fisqP@DUF@UoERYYY]UW{UQfhJBJjij`@MjDyitZdXySkNBx";

            String idc_x_p3 = "fkwqP@DIZ@UoQrJIQJIYIIQGIQQQQfhJBJijjjh@C[ANZjaceN@";
            String idc_x_p4 = "fkwqP@DIZ@UoQrJIQJIYIIQGIQQQQfhJBJijjjh@C[ANZjaceN@";
            String idc_x_p5 = "fkwqP@DIZ@UoQrJIQJIYIIQGIQQQQfhJBJijjjh@C[ANZjaceN@";

            molecules.add(HyperspaceUtils.parseIDCode(idc_0));
            molecules.add(HyperspaceUtils.parseIDCode(idc_a));
            molecules.add(HyperspaceUtils.parseIDCode(idc_b));
            molecules.add(HyperspaceUtils.parseIDCode(idc_c));
            molecules.add(HyperspaceUtils.parseIDCode(idc_d));
            molecules.add(HyperspaceUtils.parseIDCode(idc_e));
            molecules.add(HyperspaceUtils.parseIDCode(idc_f));
            molecules.add(HyperspaceUtils.parseIDCode(idc_g));

            molecules.add(HyperspaceUtils.parseIDCode(idc_x_p3));
            molecules.add(HyperspaceUtils.parseIDCode(idc_p4));

        }

        if(false) {
            String idc_0 = "fg\u007Fa`@F`QgHheECDeLhddmihihu^mU@PUAPQUTDEPD";
            String idc_a = "ff}`@@@YIEEDuEEDdehesVjjj`@@@@@`";
            String idc_b = "ff}`@@@YIEEDuEEDdehesVjjj`@@@HRPODFb";
        }




        SynthonSpace space = null;

        try {
            //space = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data");
            // todo: update HyperspaceIOUtils.loadSynthonSpace function to include zip..
            //space = HyperspaceIOUtils.loadSynthonSpace("/home/liphath1/dev5/git/hyperspace/REAL_Space_latest_FragFp.data");
            space = HyperspaceIOUtils.loadSynthonSpace(path_to_datafile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Test_RunSubstructureSearchesForPublication.runSearchesAndCreateReport(space,molecules);


    }

}
