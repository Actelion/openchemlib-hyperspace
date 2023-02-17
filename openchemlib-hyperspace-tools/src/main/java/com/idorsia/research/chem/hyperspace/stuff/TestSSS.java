package com.idorsia.research.chem.hyperspace.stuff;

import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser2;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser;
import com.idorsia.research.chem.hyperspace.util.ConnectedPairSampler;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class TestSSS {

    private static String template = "{\n" +
            "  \"ServiceProviders\": [\n" +
            "    {\"ServiceName\":  \"TEST Space SSS\", \"ServiceProvider\":  \"HyperspaceSSS\",\n" +
            "      \"Config\": { \"SpaceName\": \"TEST Space\" , \"File\": \"test_space.data\" , \"MaxNumberOfThreads\": 1 }\n" +
            "    }]}";



    /**
     * Takes two parameters: 1: synthon space 2: MODE  [3: file containing idcodes with test queries, OR number of test queries per synthon reaction]
     * MODE can be either "runProvidedQueries" or createFragmentQueries
     *
     *
     * The program does:
     * 1. parse the synthon space
     * 2. initialize hyperspace
     * 3. create all assembled synthons and initialize normal substructure search
     * 4. runs the test queries and checks if the results are the same.
     *
     *
     *
     * @param args
     */
    public static void main(String args[]) {
        String synthonspace = args[0];
        String modus        = args[1];



        File f_synthonspace = new File(synthonspace);
        SynthonSpace sspace = null;
        try {
            sspace = SynthonSpaceParser2.parseBuildingBlocks(f_synthonspace, "FragFp", 2, 3);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        if(modus.equals("GUI")) {
            //Start gui, but first store space on disk, then load in gui..
            HyperspaceIOUtils.saveSynthonSpace(sspace,"test_space.data");

            try {
                FileWriter out = new FileWriter(new File("hyperspace_test_config.json"));
                out.write(template);
                out.flush();
                out.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            String arguments[] = new String[] {"hyperspace_test_config.json"};
            HyperspaceSearchGUI.main(arguments);

            return;
        }

        AssembledMolecules allMols_a = assembleAllMolecules(sspace);
        List<String> allMols = allMols_a.all_molecules;
        Map<String,String> allMols_IDs = allMols_a.all_mols_synthons;

        if(true) {
            for(String mi : allMols) {
                System.out.println(mi);
            }
        }

        List<String> queries = new ArrayList<>();
        if(modus.equals("runProvidedQueries")) {
            String testqueries = args[2];
            try( BufferedReader in = new BufferedReader(new FileReader(new File(testqueries))) ) {
                String li = null;
                while( (li=in.readLine()) != null ) {
                    queries.add(li);
                }
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        else if(modus.equals("createFragmentQueries")) {
            int num_queries_per_rxn = 8;
            Random random = new Random();
            for(String rxn_i : allMols_a.rxn_to_molecules.keySet()) {
                List<String> mols_ri = new ArrayList<>(allMols_a.rxn_to_molecules.get(rxn_i));
                //Collections.shuffle(mols_ri,random);
                for(int zi=0;zi<num_queries_per_rxn;zi++) {
                    String full_mi = mols_ri.get(random.nextInt(mols_ri.size()));
                    StereoMolecule smi = HyperspaceUtils.parseIDCode(full_mi);
                    ConnectedPairSampler sampler = new ConnectedPairSampler(smi);

                    if(true) {
                        List<ConnectedPairSampler.ConnectedPairOfSubgraphs> rsi = sampler.computePairs(6, 8, 1, 3);
                        //for(ConnectedPairSampler.ConnectedPairOfSubgraphs xi : rsi) {
                        int bridge_length = 1 + random.nextInt(2);
                        StereoMolecule qfi = rsi.get(random.nextInt(rsi.size())).getPairAsMoleculeWithBridgeQF(bridge_length, bridge_length);
                        qfi.setFragment(true);
                        queries.add(qfi.getIDCode());
                        System.out.println("QF: " + qfi.getIDCode());
                        //}
                    }
                    if(false) {
                        StereoMolecule ffi = sampleRandomFragment(smi,12+random.nextInt(6),random);
                        ffi.setFragment(true);
                        if(ffi != null) {
                            queries.add(ffi.getIDCode());
                            System.out.println("QF: " + ffi.getIDCode());
                        }
                    }
                }
            }
        }

        // run test
        runTests(sspace,allMols,allMols_IDs,queries);
    }

    public static StereoMolecule sampleRandomFragment(StereoMolecule m, int size, Random r) {
        m.ensureHelperArrays(Molecule.cHelperCIP);
        if(m.getAtoms()<size) {return null;}
        if(m.getAtoms()==size) {return m;}

        BitSet bsi = new BitSet();
        for(int zi=0;zi<size;zi++) {
            if(zi==0) {
                bsi.set(r.nextInt(m.getAtoms()));
                continue;
            }
            int extended_atom = -1;
            BitSet candidates = new BitSet();
            bsi.stream().forEach( xi -> HyperspaceUtilsOCL.findNeighbors(m,xi).stream().forEach( ni -> candidates.set(ni) ) );
            candidates.andNot(bsi);
            if(candidates.isEmpty()) {
                break;
            }
            List<Integer> ci = candidates.stream().boxed().collect(Collectors.toList());
            int selected = ci.get(r.nextInt(ci.size()));
            bsi.set(selected);
        }
        StereoMolecule frag = new StereoMolecule();
        boolean selection[] = new boolean[m.getAtoms()];
        bsi.stream().forEach( xi -> selection[xi] = true);
        m.copyMoleculeByAtoms(frag,selection,true,null);
        frag.ensureHelperArrays(Molecule.cHelperCIP);
        return frag;
    }

    public static void runTests(SynthonSpace sspace, List<String> assembled, Map<String,String> assembled_ids, List<String> queries) {
        FastSubstructureSearcher fss = new FastSubstructureSearcher();
        Map<String,String> ids = assembled_ids;//new HashMap<>(); queries.stream().forEach(xi -> ids.put(xi,xi));
        fss.initWithIDCodes(assembled,ids,2,false);

        List<ResultComparison> results = new ArrayList<>();
        boolean all_ok = true;
        int number_not_ok = 0;
        for(String qi : queries) {
            ResultComparison rci = runTest(sspace,fss,qi);
            results.add(rci);
            all_ok &= rci.OK;
            if(!rci.OK) {number_not_ok++;}
        }
        System.out.println("Query results summarized:");
        for(int zi=0;zi<queries.size();zi++) {
            ResultComparison rci = results.get(zi);
            System.out.println(rci.query+ " : HS-Hits="+rci.results_Hyperspace.size()+" Std-Hits="+rci.results_Normal.size() );
        }

        System.out.println("All OK: "+all_ok);
        if(!all_ok) {
            System.out.println("Number not OK: "+number_not_ok+"/"+queries.size());
            System.out.println("Look at disagreeing results:");
            for(int zi=0;zi<queries.size();zi++) {
                ResultComparison rci = results.get(zi);
                if(!rci.OK) {
                    System.out.println(rci.query + " : HS-Hits=" + rci.results_Hyperspace.size() + " Std-Hits=" + rci.results_Normal.size());
                    System.out.println("Results HS:");
                    rci.results_Hyperspace.stream().forEach( ri -> System.out.println(ri) );
                    System.out.println("Results Normal:");
                    rci.results_Normal.stream().forEach( ri -> System.out.println(ri) );
                    System.out.println("Missing in HS:");
                    List<String> list_a = new ArrayList<>(rci.results_Normal);
                    list_a.removeAll(rci.results_Hyperspace);
                    list_a.stream().forEach(ri -> System.out.println(ri));
                }
            }

        }
    }

    public static ResultComparison runTest(SynthonSpace sspace, FastSubstructureSearcher fss, String qi) {
        CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
        StereoMolecule mi = HyperspaceUtils.parseIDCode(qi);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        List<SynthonSpace.CombinatorialHit> combihits = new ArrayList<>();
        List<StereoMolecule> expanded_searches = null;
        try {
            expanded_searches = SubstructureSearchHelper.expandBridgedSearches(mi);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        for(StereoMolecule qxi : expanded_searches) {
            combihits.addAll( SubstructureSearchHelper.run_substructure_search_01(sspace, cdp, qxi, 2, true, false) );
        }

        List<StereoMolecule> mol_hs = SubstructureSearchHelper.enumerateHits(combihits);
        List<String> mols_idc_hs = mol_hs.stream().map(xi -> xi.getIDCode()).distinct().sorted().collect(Collectors.toList());
        List<String> mols_idc_fss = fss.findSubstructure(cdp,mi,100000).stream().distinct().sorted().collect(Collectors.toList());

        //System.out.println("mkay, now compare for query: "+qi);
        if(!mols_idc_hs.equals(mols_idc_fss)) {
            System.out.println("Hits_SH[idcode]");
            for(String mhsi : mols_idc_hs) {
                System.out.println(mhsi);
            }
            System.out.println("Hits_FSS[idcode]\tSynthons[idcode]");
            for(String mfsi : mols_idc_fss) {
                try {
                    System.out.println(mfsi + "\t" + HyperspaceUtils.parseSmiles(fss.resolveIdcodeToIDs(Collections.singletonList(mfsi)).get(mfsi).get(0)).getIDCode());
                }
                catch(Exception ex){}
            }
        }
        return new ResultComparison(qi,mols_idc_hs,mols_idc_fss);
    }

    public static class ResultComparison {
        public final String query;
        public final List<String> results_Hyperspace;
        public final List<String> results_Normal;
        public final boolean OK;

        public ResultComparison(String query, List<String> results_Hyperspace, List<String> results_Normal) {
            this.query = query;
            this.results_Hyperspace = results_Hyperspace.stream().distinct().sorted().collect(Collectors.toList());
            this.results_Normal = results_Normal.stream().distinct().sorted().collect(Collectors.toList());
            this.OK = results_Hyperspace.equals(results_Normal);
        }
    }


    public static class AssembledMolecules {
        public final List<String> all_molecules;
        public final Map<String,String> all_mols_synthons;
        public final Map<String,List<String>> rxn_to_molecules;

        public AssembledMolecules(List<String> all_molecules, Map<String, String> all_mols_synthons, Map<String, List<String>> rxn_to_molecules) {
            this.all_molecules = all_molecules;
            this.all_mols_synthons = all_mols_synthons;
            this.rxn_to_molecules = rxn_to_molecules;
        }
    }

    /**
     *
     * @param sspace
     * @return
     */
    //public static Pair<List<String>,Map<String,String>> assembleAllMolecules(SynthonSpace sspace) {
    public static AssembledMolecules assembleAllMolecules(SynthonSpace sspace) {
        List<String> all_molecules = new ArrayList<>();
        Map<String,String> all_mols_synthons = new HashMap<>();
        Map<String,List<String>> rxn_to_molecules = new HashMap<>();
        for( String ri : sspace.getRxnIds()) {
            List<String> mols_rxn_i = new ArrayList<>();
            List<SynthonSpace.FragType> fts = new ArrayList<>(sspace.getFragTypes(ri).values());
            List<List<SynthonSpace.FragId>> rfrags = new ArrayList<>();
            for(SynthonSpace.FragType fti : fts) {
                List<SynthonSpace.FragId> frags = sspace.getSynthonSet(ri,fti.frag);
                rfrags.add(frags);
            }
            // assemble..
            if(rfrags.size()==2) {
                for(int zi=0;zi<rfrags.get(0).size();zi++) {
                    for(int zj=0;zj<rfrags.get(1).size();zj++) {
                        List<StereoMolecule> ssi = new ArrayList<>();
                        ssi.add(SynthonSpaceParser.parseIDCode(rfrags.get(0).get(zi).idcode));
                        ssi.add(SynthonSpaceParser.parseIDCode(rfrags.get(1).get(zj).idcode));
                        String as_idc = SynthonAssembler.assembleSynthons(ssi).getIDCode();
                        all_molecules.add(as_idc);
                        List<String> smiles_synthons = new ArrayList<>();
                        for(StereoMolecule xi : ssi) {
                            try{ smiles_synthons.add( new IsomericSmilesCreator(xi).getSmiles());}
                            catch (Exception ex){ex.printStackTrace();}
                        }
                        all_mols_synthons.put(as_idc,String.join(".",smiles_synthons));
                        mols_rxn_i.add(as_idc);
                    }
                }
            }
            if(rfrags.size()==3) {
                for(int zi=0;zi<rfrags.get(0).size();zi++) {
                    for(int zj=0;zj<rfrags.get(1).size();zj++) {
                        for(int zk=0;zk<rfrags.get(2).size();zk++) {
                            List<StereoMolecule> ssi = new ArrayList<>();
                            ssi.add(SynthonSpaceParser.parseIDCode(rfrags.get(0).get(zi).idcode));
                            ssi.add(SynthonSpaceParser.parseIDCode(rfrags.get(1).get(zj).idcode));
                            ssi.add(SynthonSpaceParser.parseIDCode(rfrags.get(2).get(zk).idcode));
                            String as_idc = SynthonAssembler.assembleSynthons(ssi).getIDCode();
                            all_molecules.add(as_idc);
                            List<String> smiles_synthons = new ArrayList<>();
                            for(StereoMolecule xi : ssi) {
                                try{ smiles_synthons.add( new IsomericSmilesCreator(xi).getSmiles());}
                                catch (Exception ex){ex.printStackTrace();}
                            }
                            all_mols_synthons.put(as_idc,String.join(".",smiles_synthons));
                            mols_rxn_i.add(as_idc);
                        }
                    }
                }
            }
            rxn_to_molecules.put(ri,mols_rxn_i);
        }
        return new AssembledMolecules(all_molecules,all_mols_synthons,rxn_to_molecules);//Pair.of(all_molecules,all_mols_synthons);
    }

}
