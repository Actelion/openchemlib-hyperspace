package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Traverses a SynthonSpace and applies a SynthonDownsampler set-by-set.
 */
public class SynthonDownsamplingOrchestrator {

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        return downsample(space, downsampler, request, 1);
    }

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request,
                                                int threads) {
        List<SynthonSpace.FragType> allTypes = new ArrayList<>();
        for (String rxnId : space.getRxnIds()) {
            Map<Integer, SynthonSpace.FragType> ft = space.getFragTypes(rxnId);
            if (ft != null) {
                allTypes.addAll(ft.values());
            }
        }
        return downsample(space, allTypes, downsampler, request, threads);
    }

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                Collection<SynthonSpace.FragType> fragTypes,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        return downsample(space, fragTypes, downsampler, request, 1);
    }

    public SynthonDownsamplingResult downsample(SynthonSpace space,
                                                Collection<SynthonSpace.FragType> fragTypes,
                                                SynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request,
                                                int threads) {
        int workerCount = Math.max(1, threads);
        ExecutorService executor = Executors.newFixedThreadPool(workerCount);
        List<Future<SynthonSetDownsamplingResult>> futures = new ArrayList<>();
        try {
            for (SynthonSpace.FragType fragType : fragTypes) {
                List<SynthonSpace.FragId> synthons = space.getSynthonSet(fragType.rxn_id, fragType.frag);
                if (synthons == null || synthons.isEmpty()) {
                    continue;
                }
                Callable<SynthonSetDownsamplingResult> task = () -> {
                    List<SynthonSpace.FragId> working = new ArrayList<>(synthons);
                    return downsampler.downsample(space, fragType, working, request);
                };
                futures.add(executor.submit(task));
            }
            int total = futures.size();
            int reportEvery = Math.max(1, total / 20);
            int completed = 0;
            List<SynthonSetDownsamplingResult> perSetResults = new ArrayList<>(total);
            for (Future<SynthonSetDownsamplingResult> future : futures) {
                SynthonSetDownsamplingResult result = future.get();
                completed++;
                if (total > 0 && (completed == total || completed % reportEvery == 0)) {
                    int pct = (int) Math.round(completed * 100.0 / total);
                    System.out.println("Downsampling progress: " + completed + " / " + total + " (" + pct + "%)");
                }
                if (result != null) {
                    perSetResults.add(result);
                }
            }
            return new SynthonDownsamplingResult(space, downsampler.getName(), request, perSetResults);
        } catch (Exception e) {
            throw new IllegalStateException("Failed to downsample synthon space", e);
        } finally {
            executor.shutdownNow();
        }
    }
}
