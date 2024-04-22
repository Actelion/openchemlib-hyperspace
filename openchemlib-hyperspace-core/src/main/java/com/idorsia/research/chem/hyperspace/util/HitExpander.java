package com.idorsia.research.chem.hyperspace.util;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.hyperspace.SimpleCombinatorialHit;
import com.actelion.research.chem.hyperspace.SimpleSynthon;
import com.actelion.research.chem.io.DWARFileCreator;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class HitExpander {

    public static class SimpleExpandedHit {
        private String idcAssembled;
        private String idcAssembledCoords;
        private SimpleSynthon[] synthons;
        public SimpleExpandedHit(String idcAssembled,String idcAssembledCoords, SimpleSynthon[] synthons) {
            this.idcAssembled = idcAssembled;
            this.idcAssembledCoords = idcAssembledCoords;
            this.synthons = synthons;
        }
        public String getSynthonRxn() { return synthons[0].rxnId; }
        public String getIDCodeAssembled() {return this.idcAssembled;}
    }

    public static List<SimpleExpandedHit> expandSimpleHit(SimpleCombinatorialHit hi, int max_expanded) {
        String synthon_rxn = hi.rxnId;
        List<SimpleExpandedHit> exp_hits_i = new ArrayList<>();
        if( hi.synthons.length==2) {
            exp_hits_i = expandCombiHits2(hi.synthons[0],hi.synthons[1],max_expanded);
        }
        else if( hi.synthons.length==3) {
            exp_hits_i = expandCombiHits3(hi.synthons[0],hi.synthons[1],hi.synthons[2],max_expanded);
        }

        return exp_hits_i;
    }

    /**
     *
     * @param hi
     * @param max_expanded
     * @param numThreads if -1 it uses Runtime.availableProcessors threads
     * @return
     */
    public static List<SimpleExpandedHit> expandSimpleHit_Parallel(SimpleCombinatorialHit hi, int max_expanded, int numThreads) {
        String synthon_rxn = hi.rxnId;
        List<SimpleExpandedHit> exp_hits_i = new ArrayList<>();
        if( hi.synthons.length==2) {
            exp_hits_i = expandCombiHits2_Parallel(hi.synthons[0],hi.synthons[1],max_expanded,numThreads);
        }
        else if( hi.synthons.length==3) {
            exp_hits_i = expandCombiHits3_Parallel(hi.synthons[0],hi.synthons[1],hi.synthons[2],max_expanded,numThreads);
        }

        return exp_hits_i;
    }


    public static List<SimpleExpandedHit> expandCombiHits2(SimpleSynthon[] ma, SimpleSynthon[] mb, int max_expanded) {
        List<SimpleExpandedHit> expanded = new ArrayList<>();
        //IDCodeParser icp = new IDCodeParser();
        List<StereoMolecule> mol_a = Arrays.stream(ma).map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_b = Arrays.stream(mb).map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        for(int zi=0;zi<ma.length;zi++) {
            SimpleSynthon mai = ma[zi];
            for(int zj=0;zj<mb.length;zj++) {
                if(expanded.size()>=max_expanded) {
                    break;
                }

                SimpleSynthon mbi = mb[zj];
                List<StereoMolecule> assembly = new ArrayList<>();
                assembly.add(mol_a.get(zi)); assembly.add(mol_b.get(zj));
                StereoMolecule assembled = SynthonReactor.react(assembly);
                SimpleSynthon[] synthons_i = new SimpleSynthon[]{mai , mbi };
                expanded.add( new SimpleExpandedHit(assembled.getIDCode(),assembled.getIDCoordinates(), synthons_i));
            }
        }
        return expanded;
    }


    public static List<SimpleExpandedHit> expandCombiHits3(SimpleSynthon[] ma, SimpleSynthon[] mb, SimpleSynthon[] mc, int max_expanded) {
        List<SimpleExpandedHit> expanded = new ArrayList<>();
        //IDCodeParser icp = new IDCodeParser();
        List<StereoMolecule> mol_a = Arrays.stream(ma).map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_b = Arrays.stream(mb).map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_c = Arrays.stream(mc).map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        for(int zi=0;zi<ma.length;zi++) {
            SimpleSynthon mai = ma[zi];
            for(int zj=0;zj<mb.length;zj++) {
                SimpleSynthon mbi = mb[zj];
                for(int zk=0;zk<mc.length;zk++) {
                    if(expanded.size()>=max_expanded) {
                        break;
                    }

                    SimpleSynthon mci = mc[zk];
                    List<StereoMolecule> assembly = new ArrayList<>();
                    assembly.add(mol_a.get(zi)); assembly.add(mol_b.get(zj)); assembly.add(mol_c.get(zk));
                    StereoMolecule assembled = SynthonReactor.react(assembly);
                    SimpleSynthon[] synthons_i = new SimpleSynthon[]{mai , mbi , mci };
                    expanded.add( new SimpleExpandedHit(assembled.getIDCode(),assembled.getIDCoordinates(), synthons_i));
                }
            }
        }
        return expanded;
    }


    public static List<SimpleExpandedHit> expandCombiHits2_Parallel(SimpleSynthon[] ma, SimpleSynthon[] mb, int max_expanded, int threads) {
        // Initialize the ExecutorService with a fixed thread pool. Adjust the pool size as needed.
        int numThreads = (threads > 0) ? threads : Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        List<Future<SimpleExpandedHit>> futureList = new ArrayList<>();
        List<StereoMolecule> mol_a = Arrays.stream(ma).parallel().map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_b = Arrays.stream(mb).parallel().map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());

        for (int zi = 0; zi < ma.length; zi++) {
            for (int zj = 0; zj < mb.length; zj++) {
                final int finalZi = zi;
                final int finalZj = zj;
                // Submit callable tasks to the executor
                Callable<SimpleExpandedHit> task = () -> {
                    List<StereoMolecule> assembly = Arrays.asList(mol_a.get(finalZi), mol_b.get(finalZj));
                    StereoMolecule assembled = SynthonReactor.react(assembly);
                    SimpleSynthon[] synthons_i = new SimpleSynthon[]{ma[finalZi], mb[finalZj]};
                    return new SimpleExpandedHit(assembled.getIDCode(), assembled.getIDCoordinates(), synthons_i);
                };
                // Early exit check before task submission to avoid unnecessary task submissions.
                if (futureList.size() >= max_expanded) {
                    break;
                }
                futureList.add(executor.submit(task));
            }
        }

        List<SimpleExpandedHit> expanded = new ArrayList<>();
        // Collect and process results
        for (Future<SimpleExpandedHit> future : futureList) {
            try {
                SimpleExpandedHit hit = future.get(); // This call blocks until the task is complete.
                expanded.add(hit);
                if (expanded.size() >= max_expanded) {
                    break;
                }
            } catch (Exception e) {
                // Handle exceptions appropriately
                e.printStackTrace();
            }
        }
        executor.shutdown(); // Shutdown the executor service
        return expanded;
    }

    public static List<SimpleExpandedHit> expandCombiHits3_Parallel(SimpleSynthon[] ma, SimpleSynthon[] mb, SimpleSynthon[] mc, int max_expanded, int threads) {
        // Initialize the ExecutorService with a fixed thread pool. Adjust the pool size as needed.
        int numThreads = (threads > 0) ? threads : Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        List<Future<SimpleExpandedHit>> futureList = new ArrayList<>();
        List<StereoMolecule> mol_a = Arrays.stream(ma).parallel().map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_b = Arrays.stream(mb).parallel().map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());
        List<StereoMolecule> mol_c = Arrays.stream(mc).parallel().map(si -> HyperspaceUtils.parseIDCode(si.idcode)).collect(Collectors.toList());

        for (int zi = 0; zi < ma.length; zi++) {
            for (int zj = 0; zj < mb.length; zj++) {
                for (int zk = 0; zj < mb.length; zk++) {
                    final int finalZi = zi;
                    final int finalZj = zj;
                    final int finalZk = zk;
                    // Submit callable tasks to the executor
                    Callable<SimpleExpandedHit> task = () -> {
                        List<StereoMolecule> assembly = Arrays.asList(mol_a.get(finalZi), mol_b.get(finalZj), mol_c.get(finalZk));
                        StereoMolecule assembled = SynthonReactor.react(assembly);
                        SimpleSynthon[] synthons_i = new SimpleSynthon[]{ma[finalZi], mb[finalZj], mc[finalZk]};
                        return new SimpleExpandedHit(assembled.getIDCode(), assembled.getIDCoordinates(), synthons_i);
                    };
                    // Early exit check before task submission to avoid unnecessary task submissions.
                    if (futureList.size() >= max_expanded) {
                        break;
                    }
                    futureList.add(executor.submit(task));
                }
            }
        }

        List<SimpleExpandedHit> expanded = new ArrayList<>();
        // Collect and process results
        for (Future<SimpleExpandedHit> future : futureList) {
            try {
                SimpleExpandedHit hit = future.get(); // This call blocks until the task is complete.
                expanded.add(hit);
                if (expanded.size() >= max_expanded) {
                    break;
                }
            } catch (Exception e) {
                // Handle exceptions appropriately
                e.printStackTrace();
            }
        }
        executor.shutdown(); // Shutdown the executor service
        return expanded;
    }


    public static void exportExpandedHitsWithQueriesToDWAR(String file, List<Pair<String,SimpleExpandedHit>> hitsWithQueries) {
        exportExpandedHitsWithQueriesToDWAR(file,hitsWithQueries, new ArrayList<>());
    }

    /**
     *
     * @param file
     * @param hitsWithQueries
     * @param additional_values triple: left value: name, middle value: is_idcode, right value: values
     */
    public static void exportExpandedHitsWithQueriesToDWAR(String file, List<Pair<String,SimpleExpandedHit>> hitsWithQueries,
                                                           List<Triple<String,Boolean,List<String>>> additional_values ) {
        try(FileWriter f_out = new FileWriter(file)) {
            BufferedWriter bw_out = new BufferedWriter(f_out);
            DWARFileCreator out = new DWARFileCreator(bw_out);
            int idx_s_query = out.addStructureColumn("Query",null);

            // put additional values in front(?)
            int indices_additional_v[]       = new int[additional_values.size()];
            for(int zi=0;zi<additional_values.size();zi++) {
                if (additional_values.get(zi).getMiddle()) {
                    indices_additional_v[zi] = out.addStructureColumn(additional_values.get(zi).getLeft(),null);
                } else {
                    indices_additional_v[zi] = out.addAlphanumericalColumn(additional_values.get(zi).getLeft());
                }
            }


            int idx_s_hit   = out.addStructureColumn("Structure",null);
            int idx_srxn    = out.addAlphanumericalColumn("SynthonRxn");
            int idx_s_a     = out.addStructureColumn("Synthon_A",null);
            int idx_s_b     = out.addStructureColumn("Synthon_B",null);
            int idx_s_c     = out.addStructureColumn("Synthon_C",null);
            int idx_id_a    = out.addAlphanumericalColumn("SynthonId_A");
            int idx_id_b    = out.addAlphanumericalColumn("SynthonID_B");
            int idx_id_c    = out.addAlphanumericalColumn("SynthonID_C");

            out.writeHeader(hitsWithQueries.size());

            String idc_empty = (new StereoMolecule()).getIDCode();

            for(int zi=0;zi<hitsWithQueries.size();zi++) {
                Pair<String,SimpleExpandedHit> hi = hitsWithQueries.get(zi);
                out.setRowStructure(hi.getLeft(),idx_s_query);
                out.setRowStructure(hi.getRight().idcAssembled,idx_s_hit);
                out.setRowValue(hi.getRight().getSynthonRxn(),idx_srxn);
                out.setRowStructure(hi.getRight().synthons[0].idcode,idx_s_a);
                out.setRowStructure( ((hi.getRight().synthons.length>=2)?hi.getRight().synthons[1].idcode:idc_empty) ,idx_s_b );
                out.setRowStructure( ((hi.getRight().synthons.length>=3)?hi.getRight().synthons[2].idcode:idc_empty) ,idx_s_c );
                out.setRowValue(hi.getRight().synthons[0].synthonId,idx_id_a);
                out.setRowStructure( ((hi.getRight().synthons.length>=2)?hi.getRight().synthons[1].synthonId:"") ,idx_id_b );
                out.setRowStructure( ((hi.getRight().synthons.length>=3)?hi.getRight().synthons[2].synthonId:"") ,idx_id_c );

                // additional values:
                for(int zj=0;zj<indices_additional_v.length;zj++) {
                    if (additional_values.get(zj).getMiddle()) {
                        out.setRowStructure(additional_values.get(zj).getRight().get(zi), indices_additional_v[zj]);
                    } else {
                        out.setRowValue(additional_values.get(zj).getRight().get(zi), indices_additional_v[zj]);
                    }
                }

                out.writeCurrentRow();
            }
            out.writeEnd();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
