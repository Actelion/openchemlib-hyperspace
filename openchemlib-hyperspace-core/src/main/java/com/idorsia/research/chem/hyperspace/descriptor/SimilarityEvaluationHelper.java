package com.idorsia.research.chem.hyperspace.descriptor;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;

public class SimilarityEvaluationHelper {

    public static class EvaluatedSimilarities {
        private List<String> mols;

        private Map<String,Map<String,Object>> descriptors;
        private Map<String,Map <String, Map<String,Double>>> similarities;

        public EvaluatedSimilarities(List<String> mols, Map<String,Map<String,Object>> descriptors, Map<String, Map<String, Map<String, Double>>> similarities) {
            this.mols = mols;
            this.descriptors = descriptors;
            this.similarities = similarities;
        }


        public List<String> getMols() {
            return mols;
        }

        /**
         * desciprotrs.get(i).get(j) contains j th descriptor for i'th lib molecule
         */
        public Map<String,Map<String,Object>> getDescriptors() {
            return descriptors;
        }


        /**
         * Map: descriptor-shortname -> ( refmol-idc -> ( libmol-idcode -> similarity )
         */
        public Map<String, Map<String, Map<String, Double>>> getSimilarities() {
            return similarities;
        }

    }

    public static EvaluatedSimilarities evaluateSimilarities(List<StereoMolecule> refmols , List<String> libmols, List<DescriptorHandler> dhs, int threads, int batchsize_descriptor, int batchsize_similarity) {

        if(threads<=0) {
            threads = Runtime.getRuntime().availableProcessors();
        }


        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        List<Future<String>> jobs = new ArrayList<>();

        Map<String,Map<String,Object>> all_descriptors = new HashMap<>();

        // map: fp-shortname -> ( refmol-idc -> ( libmol-idc -> similarity ) )
        Map<String,Map <String, Map<String,Double>>> all_similarities = new HashMap<>();

        for(DescriptorHandler dhi : dhs) {

            Map<String,Object> descriptors_i = new ConcurrentHashMap<>();

            List<Future> jobs_descr = new ArrayList<>();

            int batch_start = 0;
            int zbi = 0;
            //for(int zbi=0;zbi<libmols.size();zbi+=batchsize_descriptor) {
            do {
                int final_batch_start = batch_start;
                zbi += batchsize_descriptor;
                int batch_end = Math.min(libmols.size(),zbi);

                Runnable ri = new Runnable() {
                    @Override
                    public void run() {
                        IDCodeParser icp = new IDCodeParser();
                        for(int zi=final_batch_start;zi<batch_end;zi++) {

                            StereoMolecule mi = new StereoMolecule();
                            icp.parse(mi, libmols.get(zi));
                            mi.ensureHelperArrays(Molecule.cHelperCIP);

                            Object fp_mi = dhi.createDescriptor(mi);
                            descriptors_i.put(libmols.get(zi),fp_mi);
                        }
                    }
                };

                jobs_descr.add(executorService.submit(ri));
                batch_start = batch_end;
            } while(zbi<libmols.size());

            // wait for jobs:
            for(Future fi : jobs_descr) {
                try {
                    fi.get();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                }
            }

            all_descriptors.put(dhi.getInfo().shortName,descriptors_i);

            Map<String,Map<String,Double>> similarities_dh_i = new HashMap<>();

            for(StereoMolecule rmi : refmols) {
                rmi.ensureHelperArrays(Molecule.cHelperCIP);

                Object fp_rmi = dhi.createDescriptor(rmi);
                Map<String,Double> similarities_i = new ConcurrentHashMap<>();

                List<Future> jobs_sim = new ArrayList<>();

                batch_start = 0;
                zbi = 0;
                //for(int zbi=0;zbi<libmols.size();zbi+=batchsize_similarity) {

                do{
                    int final_batch_start = batch_start;
                    zbi+= batchsize_similarity;

                    int batch_end = Math.min(libmols.size(),zbi);

                    Runnable ri = new Runnable() {
                        @Override
                        public void run() {
                            for(int zi=final_batch_start;zi<batch_end;zi++) {
                                Object fp_mi = descriptors_i.get(libmols.get(zi));

                                if(fp_mi!=null) {
                                    double sim_i = dhi.getSimilarity(fp_rmi, fp_mi);
                                    similarities_i.put(libmols.get(zi), sim_i);
                                }
                                else {
                                    double sim_i = Double.NaN;
                                    similarities_i.put(libmols.get(zi), sim_i);
                                }
                            }
                        }
                    };

                    jobs_sim.add( executorService.submit(ri) );

                    batch_start = batch_end;
                } while (zbi<libmols.size());

                // wait for jobs:
                for(Future fi : jobs_sim) {
                    try {
                        fi.get();
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    } catch (ExecutionException e) {
                        e.printStackTrace();
                    }
                }
                //all_descriptors.add(descriptors_i);
                similarities_dh_i.put(rmi.getIDCode(),similarities_i);
            }
            all_similarities.put(dhi.getInfo().shortName,similarities_dh_i);
        }

        executorService.shutdown();

        return new EvaluatedSimilarities(libmols,all_descriptors,all_similarities);
    }

}
