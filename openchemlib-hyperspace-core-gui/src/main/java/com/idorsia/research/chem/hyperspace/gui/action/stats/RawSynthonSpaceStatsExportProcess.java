package com.idorsia.research.chem.hyperspace.gui.action.stats;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.gui.process.AbstractHyperspaceProcess;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsExporter;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsOptions;
import com.idorsia.research.chem.hyperspace.stats.RawSynthonSpaceStatsReport;

import java.nio.file.Path;

public final class RawSynthonSpaceStatsExportProcess extends AbstractHyperspaceProcess implements AbstractHyperspaceProcess.HasProgress {
    private final Path rawInput;
    private final Path outputDirectory;
    private final int examplesPerReaction;
    private final int productSamplesPerReaction;
    private final long seed;
    private final int threads;
    private volatile double progress;

    public RawSynthonSpaceStatsExportProcess(Path rawInput,
                                             Path outputDirectory,
                                             int examplesPerReaction,
                                             int productSamplesPerReaction,
                                             long seed,
                                             int threads) {
        this.rawInput = rawInput;
        this.outputDirectory = outputDirectory;
        this.examplesPerReaction = examplesPerReaction;
        this.productSamplesPerReaction = productSamplesPerReaction;
        this.seed = seed;
        this.threads = threads;
    }

    public void startAsync() {
        setProcessStatus(ProcessStatus.COMPUTING);
        Thread thread = new Thread(this::runExport, "RawspaceStatsExport");
        thread.setDaemon(true);
        thread.start();
    }

    private void runExport() {
        try {
            setProcessStatusMessage("Loading " + rawInput.getFileName());
            RawSynthonSpace rawSpace = HyperspaceIOUtils.loadRawSynthonSpace(rawInput.toString());
            RawSynthonSpaceStatsOptions options = RawSynthonSpaceStatsOptions.builder()
                    .examplesPerReaction(examplesPerReaction)
                    .productSamplesPerReaction(productSamplesPerReaction)
                    .seed(seed)
                    .threads(threads)
                    .progressListener((completed, total, reactionId) -> {
                        progress = total <= 0 ? 1.0 : (double) completed / total;
                        setProcessStatusMessage("Processed " + completed + " / " + total + " reactions");
                    })
                    .build();
            RawSynthonSpaceStatsReport report = RawSynthonSpaceStatsExporter.analyze(rawSpace, options);
            setProcessStatusMessage("Writing report");
            report.writeToDirectory(outputDirectory);
            progress = 1.0;
            setProcessStatusMessage("Wrote " + outputDirectory);
            setProcessStatus(ProcessStatus.DONE);
        } catch (Exception ex) {
            progress = 1.0;
            setProcessStatusMessage(ex.getMessage() == null ? ex.getClass().getSimpleName() : ex.getMessage());
            setProcessStatus(ProcessStatus.FAILED);
        }
    }

    @Override
    public String getName() {
        return "Export Rawspace Stats";
    }

    @Override
    public double getProgress() {
        return progress;
    }
}
