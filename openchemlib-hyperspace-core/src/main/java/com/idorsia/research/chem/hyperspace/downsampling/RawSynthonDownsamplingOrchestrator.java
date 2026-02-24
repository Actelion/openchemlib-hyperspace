package com.idorsia.research.chem.hyperspace.downsampling;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSet;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Applies a {@link RawSynthonDownsampler} to every synthon set in a {@link RawSynthonSpace}.
 */
public class RawSynthonDownsamplingOrchestrator {

    public SynthonDownsamplingResult downsample(RawSynthonSpace rawSpace,
                                                RawSynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request) {
        return downsample(rawSpace, downsampler, request, 1);
    }

    public SynthonDownsamplingResult downsample(RawSynthonSpace rawSpace,
                                                RawSynthonDownsampler downsampler,
                                                SynthonDownsamplingRequest request,
                                                int threads) {
        int workerCount = Math.max(1, threads);
        ExecutorService executor = Executors.newFixedThreadPool(workerCount);
        List<Future<SynthonSetDownsamplingResult>> futures = new ArrayList<>();
        try {
            for (Map.Entry<String, RawSynthonSpace.ReactionData> entry : rawSpace.getReactions().entrySet()) {
                String reactionId = entry.getKey();
                entry.getValue().getRawFragmentSets().forEach((fragIdx, synthons) -> {
                    if (synthons == null || synthons.isEmpty()) {
                        return;
                    }
                    Callable<SynthonSetDownsamplingResult> task = () -> {
                        List<RawSynthon> working = new ArrayList<>(synthons);
                        RawSynthonSet set = new RawSynthonSet(reactionId, fragIdx, working);
                        return downsampler.downsample(set, request);
                    };
                    futures.add(executor.submit(task));
                });
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
            return new SynthonDownsamplingResult(downsampler.getName(), request, perSetResults);
        } catch (Exception e) {
            throw new IllegalStateException("Failed to downsample raw synthon space", e);
        } finally {
            executor.shutdownNow();
        }
    }
}
