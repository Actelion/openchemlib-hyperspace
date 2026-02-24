package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.screening.CandidateSampler;
import com.idorsia.research.chem.hyperspace.screening.ContinuousScreeningOrchestrator;
import com.idorsia.research.chem.hyperspace.screening.ScreeningMetrics;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.nio.file.Path;

public final class ContinuousScreeningCLI {

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("ContinuousScreeningCLI", options, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        RawSynthonSpace fullRaw = HyperspaceIOUtils.loadRawSynthonSpace(cmd.getOptionValue("rawFull"));
        RawSynthonSpace downsampledRaw = HyperspaceIOUtils.loadRawSynthonSpace(cmd.getOptionValue("rawDownsampled"));
        PheSAMolecule queryDescriptor = buildQueryDescriptor(cmd);

        CandidateSampler.Config samplerConfig = new CandidateSampler.Config(
                Long.parseLong(cmd.getOptionValue("candidateAttempts", "1000")),
                Integer.parseInt(cmd.getOptionValue("minAtoms", "1")),
                Integer.parseInt(cmd.getOptionValue("maxAtoms", "0")),
                Integer.parseInt(cmd.getOptionValue("maxRotatable", "15")),
                Double.parseDouble(cmd.getOptionValue("candidateThreshold", "0.6"))
        );

        LocalOptimizationRequest fullRequest = LocalOptimizationRequest.builder()
                .beamSize(Integer.parseInt(cmd.getOptionValue("fullBeam", "10")))
                .neighborPoolSize(Integer.parseInt(cmd.getOptionValue("fullTopL", "100")))
                .sampledNeighbors(Integer.parseInt(cmd.getOptionValue("fullSampleNeighbors", "10")))
                .perPositionCap(Integer.parseInt(cmd.getOptionValue("fullCap", "2")))
                .maxRounds(Integer.parseInt(cmd.getOptionValue("fullRounds", "6")))
                .patience(Integer.parseInt(cmd.getOptionValue("fullPatience", "3")))
                .minPhesaSimilarity(Double.parseDouble(cmd.getOptionValue("fullMinSimilarity", "0.6")))
                .minScoreThreshold(0.0)
                .reportAllCandidates(true)
                .randomSeed(Long.parseLong(cmd.getOptionValue("randomSeed", "13")))
                .build();

        boolean microEnabled = cmd.hasOption("microEnabled");
        LocalOptimizationRequest microRequest = null;
        if (microEnabled) {
            microRequest = LocalOptimizationRequest.builder()
                    .beamSize(Integer.parseInt(cmd.getOptionValue("microBeam", "5")))
                    .neighborPoolSize(Integer.parseInt(cmd.getOptionValue("microTopL", "50")))
                    .sampledNeighbors(Integer.parseInt(cmd.getOptionValue("microSampleNeighbors", "5")))
                    .perPositionCap(Integer.parseInt(cmd.getOptionValue("microCap", "2")))
                    .maxRounds(Integer.parseInt(cmd.getOptionValue("microRounds", "2")))
                    .patience(Integer.parseInt(cmd.getOptionValue("microPatience", "1")))
                    .minPhesaSimilarity(Double.parseDouble(cmd.getOptionValue("microMinSimilarity", "0.6")))
                    .randomSeed(Long.parseLong(cmd.getOptionValue("randomSeed", "13")))
                    .build();
        }

        ContinuousScreeningOrchestrator.Config config = new ContinuousScreeningOrchestrator.Config()
                .withFullRaw(fullRaw)
                .withDownsampledRaw(downsampledRaw)
                .withQueryDescriptor(queryDescriptor)
                .withSamplerConfig(samplerConfig)
                .withFullOptimizationRequest(fullRequest)
                .withFullOptimizerThreads(Integer.parseInt(cmd.getOptionValue("fullThreads", "4")))
                .withQueueCapacity(Integer.parseInt(cmd.getOptionValue("queueCapacity", "1000")))
                .withProgressIntervalSeconds(Integer.parseInt(cmd.getOptionValue("progressSeconds", "60")))
                .withReactionWeightExponent(Double.parseDouble(cmd.getOptionValue("reactionWeightExponent", "1.0")))
                .withReactionMinWeight(Double.parseDouble(cmd.getOptionValue("reactionMinWeight", "0.01")))
                .withHitOutput(Path.of(cmd.getOptionValue("outputHits")))
                .withDuplicateCacheSize(Integer.parseInt(cmd.getOptionValue("dedupeMax", "200000")))
                .withRandomSeed(Long.parseLong(cmd.getOptionValue("randomSeed", "13")));
        if (microEnabled) {
            config.withMicroEnabled(true).withMicroOptimizationRequest(microRequest);
        }
        config.validate();

        long iterations = Long.parseLong(cmd.getOptionValue("iterations", "-1"));
        try (ContinuousScreeningOrchestrator orchestrator = new ContinuousScreeningOrchestrator(config)) {
            orchestrator.run(iterations);
            ScreeningMetrics metrics = orchestrator.getMetrics();
            System.out.println("Sampled candidates: " + metrics.getSampled());
            System.out.println("Duplicate seeds skipped: " + metrics.getDuplicateSeeds());
            System.out.println("Submitted to optimizer: " + metrics.getSubmitted());
            System.out.println("Total hits: " + metrics.getHits());
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawFull").hasArg().required(true)
                .desc("Path to full RawSynthonSpace (.rawspace/.gz)").build());
        options.addOption(Option.builder().longOpt("rawDownsampled").hasArg().required(true)
                .desc("Path to downsampled RawSynthonSpace").build());
        options.addOption(Option.builder().longOpt("querySmiles").hasArg()
                .desc("Query ligand as SMILES").build());
        options.addOption(Option.builder().longOpt("queryIdcode").hasArg()
                .desc("Query ligand as IDCode").build());
        options.addOption(Option.builder().longOpt("outputHits").hasArg().required(true)
                .desc("Output TSV for optimized hits").build());
        options.addOption(Option.builder().longOpt("iterations").hasArg()
                .desc("Number of sampling iterations (-1 = infinite)").build());
        options.addOption(Option.builder().longOpt("candidateAttempts").hasArg()
                .desc("Attempts per reaction when sampling candidates").build());
        options.addOption(Option.builder().longOpt("minAtoms").hasArg()
                .desc("Minimum atoms for sampled assembly").build());
        options.addOption(Option.builder().longOpt("maxAtoms").hasArg()
                .desc("Maximum atoms (0 = unlimited)").build());
        options.addOption(Option.builder().longOpt("maxRotatable").hasArg()
                .desc("Maximum rotatable bonds").build());
        options.addOption(Option.builder().longOpt("candidateThreshold").hasArg()
                .desc("Minimum similarity for sampled candidate").build());
        options.addOption(Option.builder().longOpt("fullBeam").hasArg().desc("Full optimizer beam size").build());
        options.addOption(Option.builder().longOpt("fullTopL").hasArg().desc("Full optimizer neighbor pool size").build());
        options.addOption(Option.builder().longOpt("fullSampleNeighbors").hasArg().desc("Sampled neighbors per expansion").build());
        options.addOption(Option.builder().longOpt("fullCap").hasArg().desc("Full optimizer per-position cap").build());
        options.addOption(Option.builder().longOpt("fullRounds").hasArg().desc("Full optimizer max rounds").build());
        options.addOption(Option.builder().longOpt("fullPatience").hasArg().desc("Full optimizer patience").build());
        options.addOption(Option.builder().longOpt("fullMinSimilarity").hasArg().desc("Full optimizer min similarity").build());
        options.addOption(Option.builder().longOpt("fullThreads").hasArg().desc("Worker threads for full optimizer").build());
        options.addOption(Option.builder().longOpt("queueCapacity").hasArg()
                .desc("Max queued sample+opt jobs before applying backpressure").build());
        options.addOption(Option.builder().longOpt("progressSeconds").hasArg()
                .desc("Interval in seconds for progress logging (0 to disable)").build());
        options.addOption(Option.builder().longOpt("reactionWeightExponent").hasArg()
                .desc("Exponent applied to reaction size weighting (default 1.0)").build());
        options.addOption(Option.builder().longOpt("reactionMinWeight").hasArg()
                .desc("Minimum weight as a fraction of the max reaction weight (default 0.01)").build());
        options.addOption(Option.builder().longOpt("microEnabled").desc("Enable micro optimization on downsampled space").build());
        options.addOption(Option.builder().longOpt("microBeam").hasArg().desc("Micro optimizer beam size").build());
        options.addOption(Option.builder().longOpt("microTopL").hasArg().desc("Micro optimizer neighbor pool size").build());
        options.addOption(Option.builder().longOpt("microSampleNeighbors").hasArg().desc("Micro optimizer sampled neighbors").build());
        options.addOption(Option.builder().longOpt("microCap").hasArg().desc("Micro optimizer per-position cap").build());
        options.addOption(Option.builder().longOpt("microRounds").hasArg().desc("Micro optimizer rounds").build());
        options.addOption(Option.builder().longOpt("microPatience").hasArg().desc("Micro optimizer patience").build());
        options.addOption(Option.builder().longOpt("microMinSimilarity").hasArg().desc("Micro optimizer min similarity").build());
        options.addOption(Option.builder().longOpt("dedupeMax").hasArg()
                .desc("Maximum entries stored in duplicate filter").build());
        options.addOption(Option.builder().longOpt("randomSeed").hasArg()
                .desc("Random seed").build());
        return options;
    }

    private static PheSAMolecule buildQueryDescriptor(CommandLine cmd) throws Exception {
        StereoMolecule molecule = parseQuery(cmd);
        ConformerSetGenerator generator = new ConformerSetGenerator(1);
        ConformerSet conformers = generator.generateConformerSet(molecule);
        if (conformers.isEmpty()) {
            throw new IllegalStateException("Failed to generate conformer for query molecule");
        }
        DescriptorHandlerShape handler = new DescriptorHandlerShape();
        PheSAMolecule descriptor = handler.createDescriptor(conformers);
        if (handler.calculationFailed(descriptor)) {
            throw new IllegalStateException("Unable to create PheSA descriptor for query molecule");
        }
        return descriptor;
    }

    private static StereoMolecule parseQuery(CommandLine cmd) throws Exception {
        String smiles = cmd.getOptionValue("querySmiles");
        String idcode = cmd.getOptionValue("queryIdcode");
        if (smiles == null && idcode == null) {
            throw new IllegalArgumentException("Provide --querySmiles or --queryIdcode");
        }
        if (smiles != null && idcode != null) {
            throw new IllegalArgumentException("Specify only one of --querySmiles or --queryIdcode");
        }
        StereoMolecule molecule = new StereoMolecule();
        if (smiles != null) {
            new SmilesParser().parse(molecule, smiles);
        } else {
            new IDCodeParser().parse(molecule, idcode);
        }
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule;
    }
}
