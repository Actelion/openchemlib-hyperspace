package com.idorsia.research.chem.hyperspace3d.benchmark;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace3d.batch.*;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.model.*;
import com.idorsia.research.chem.hyperspace3d.screening.ScreeningObjective;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public final class Hyperspace3DMicrobenchmark {
    private Hyperspace3DMicrobenchmark() {}

    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("usage: benchmark MODEL_BUNDLE");
        StereoMolecule molecule = new StereoMolecule();
        new SmilesParser().parse(molecule, "CCOc1ccc2nc(N)sc2c1");
        OCLDeepSpaceFeaturizer featurizer = new OCLDeepSpaceFeaturizer();
        DeepSpaceTensorBatchBuilder builder = new DeepSpaceTensorBatchBuilder(featurizer);
        int warmup = 20, iterations = 100;
        for (int i = 0; i < warmup; i++) featurizer.featurize(molecule);

        double assembly = measure(iterations,
                () -> SynthonAssembler.assembleSynthons_faster(List.of(molecule)));
        double features = measure(iterations, () -> featurizer.featurize(molecule));
        List<StereoMolecule> batchMolecules = new ArrayList<>();
        List<FeaturizationResult> batchFeatures = new ArrayList<>();
        for (int i = 0; i < 64; i++) {
            batchMolecules.add(molecule);
            batchFeatures.add(featurizer.featurize(molecule));
        }
        double packing = measure(iterations, () -> builder.buildFromFeatures(batchFeatures));

        DeepSpaceModelBundle bundle = DeepSpaceModelBundle.load(Path.of(args[0]));
        DeepSpaceOnnxEnvironment runtime =
                new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.CPU);
        try (DeepSpaceV1Encoder encoder = new DeepSpaceV1Encoder(runtime, bundle);
             DeepSpaceV1Comparator comparator = new DeepSpaceV1Comparator(runtime, bundle)) {
            DeepSpaceTensorBatch tensors = builder.build(batchMolecules);
            float[][] encoded = encoder.encode(tensors);
            float[] query = encoded[0];
            double encoding = measure(20, () -> encoder.encode(tensors));
            double comparison = measure(iterations, () -> comparator.compare(query, encoded));
            int comparisonBatch = 32768;
            float[] flatCandidates = new float[comparisonBatch * 128];
            for (int row = 0; row < comparisonBatch; row++) {
                System.arraycopy(query, 0, flatCandidates, row * 128, 128);
            }
            comparator.compareFlat(query, flatCandidates, comparisonBatch);
            double largeComparison = measure(5, () -> comparator.compareFlat(
                    query, flatCandidates, comparisonBatch));
            ProductTuple tuple = new ProductTuple("benchmark", List.of("s0"), List.of(0));
            DeepSpaceBatchAssemblyScorer scorer = new DeepSpaceBatchAssemblyScorer(
                    ignored -> List.of(molecule), featurizer, encoder, comparator,
                    bundle.manifest(), ScreeningObjective.direct("phesa_total"), query);
            List<AssemblyCandidate> candidates = new ArrayList<>();
            for (int i = 0; i < 64; i++) {
                candidates.add(new AssemblyCandidate(i,
                        new ProductTuple("benchmark", List.of("s" + i), List.of(i)),
                        "bench", 0, Map.of()));
            }
            double complete = measure(10, () -> scorer.scoreBatch(candidates));
            System.out.printf("assembly_cpu_us_per_molecule=%.3f%n", assembly / iterations / 1e3);
            System.out.printf("featurization_cpu_us_per_molecule=%.3f%n", features / iterations / 1e3);
            System.out.printf("tensor_packing_cpu_us_per_64=%.3f%n", packing / iterations / 1e3);
            System.out.printf("onnx_encoding_cpu_ms_per_64=%.3f%n", encoding / 20 / 1e6);
            System.out.printf("comparator_cpu_us_per_64=%.3f%n", comparison / iterations / 1e3);
            System.out.printf("comparator_cpu_molecules_per_second_batch_32768=%.0f%n",
                    comparisonBatch * 5.0 / (largeComparison / 1e9));
            System.out.printf("score_batch_warm_cache_cpu_ms_per_64=%.3f%n", complete / 10 / 1e6);
        }
    }

    private static double measure(int iterations, Runnable operation) {
        long start = System.nanoTime();
        for (int i = 0; i < iterations; i++) operation.run();
        return System.nanoTime() - start;
    }
}
