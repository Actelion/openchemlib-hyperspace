package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.localopt.LocalBeamOptimizer;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationLogLevel;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.localopt.PheSAAssemblyScorer;
import com.idorsia.research.chem.hyperspace.localopt.SeedAssembly;
import com.idorsia.research.chem.hyperspace.localopt.SeedTsvParser;
import com.idorsia.research.chem.hyperspace.localopt.SynthonSetAccessor;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public class SynthonLocalOptimizerCLI {

    public static void main(String[] args) throws Exception {
        Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("SynthonLocalOptimizerCLI", options, true);
            throw new IllegalArgumentException("Unable to parse CLI parameters", e);
        }

        Path seedsPath = Path.of(cmd.getOptionValue("seeds"));
        Path outputPath = Path.of(cmd.getOptionValue("output"));

        SynthonSetAccessor provider = createProvider(cmd);
        PheSAMolecule queryDescriptor = buildQueryDescriptor(cmd);

        LocalOptimizationRequest request = LocalOptimizationRequest.builder()
                .beamSize(Integer.parseInt(cmd.getOptionValue("beam", "10")))
                .neighborPoolSize(Integer.parseInt(cmd.getOptionValue("topL", "100")))
                .sampledNeighbors(Integer.parseInt(cmd.getOptionValue("sampleNeighbors", "10")))
                .perPositionCap(Integer.parseInt(cmd.getOptionValue("cap", "2")))
                .maxRounds(Integer.parseInt(cmd.getOptionValue("rounds", "6")))
                .patience(Integer.parseInt(cmd.getOptionValue("patience", "3")))
                .minPhesaSimilarity(Double.parseDouble(cmd.getOptionValue("minPhesa", "0.6")))
                .randomSeed(Long.parseLong(cmd.getOptionValue("seed", "13")))
                .logLevel(parseLogLevel(cmd.getOptionValue("log", "NONE")))
                .build();

        List<SeedAssembly> seeds = SeedTsvParser.parse(seedsPath);
        PheSAAssemblyScorer scorer = new PheSAAssemblyScorer(queryDescriptor, request.getMinPhesaSimilarity());
        LocalBeamOptimizer optimizer = new LocalBeamOptimizer(provider, scorer);
        int threads = Integer.parseInt(cmd.getOptionValue("seedThreads", "1"));
        LocalOptimizerScheduler scheduler = new LocalOptimizerScheduler(threads, optimizer, request, outputPath);
        scheduler.run(seeds);
    }

    private static SynthonSetAccessor createProvider(CommandLine cmd) throws IOException {
        String spacePath = cmd.getOptionValue("space");
        String downsampledPath = cmd.getOptionValue("downsampled");
        String rawPath = cmd.getOptionValue("raw");
        if (rawPath != null) {
            if (spacePath != null || downsampledPath != null) {
                throw new IllegalArgumentException("Use only one of --space, --downsampled, or --raw");
            }
            RawSynthonSpace rawSpace = HyperspaceIOUtils.loadRawSynthonSpace(rawPath);
            String mode = cmd.getOptionValue("rawMode", "FULL").toUpperCase();
            switch (mode) {
                case "FULL":
                    return new SynthonSetAccessor(rawSpace);
                case "DOWNSAMPLED":
                    return new SynthonSetAccessor(rawSpace.toDownsampledSpace());
                default:
                    throw new IllegalArgumentException("Unknown --rawMode value: " + mode);
            }
        }
        if (spacePath != null) {
            SynthonSpace space = HyperspaceIOUtils.loadSynthonSpace(spacePath);
            return new SynthonSetAccessor(space);
        }
        if (downsampledPath != null) {
            DownsampledSynthonSpace dspace = HyperspaceIOUtils.loadDownsampledSynthonSpace(downsampledPath);
            return new SynthonSetAccessor(dspace);
        }
        throw new IllegalArgumentException("Provide --space, --downsampled, or --raw");
    }

    private static PheSAMolecule buildQueryDescriptor(CommandLine cmd) {
        StereoMolecule molecule = parseQuery(cmd);
        ConformerSetGenerator generator = new ConformerSetGenerator(1);
        ConformerSet conformers = generator.generateConformerSet(molecule);
        if (conformers.isEmpty()) {
            throw new IllegalStateException("Failed to build conformer for query molecule");
        }
        DescriptorHandlerShape handler = new DescriptorHandlerShape();
        PheSAMolecule descriptor = handler.createDescriptor(conformers);
        if (handler.calculationFailed(descriptor)) {
            throw new IllegalStateException("Unable to create PheSA descriptor for query molecule");
        }
        return descriptor;
    }

    private static StereoMolecule parseQuery(CommandLine cmd) {
        String smiles = cmd.getOptionValue("querySmiles");
        String idcode = cmd.getOptionValue("queryIdcode");
        if (smiles == null && idcode == null) {
            throw new IllegalArgumentException("Provide --querySmiles or --queryIdcode");
        }
        if (smiles != null && idcode != null) {
            throw new IllegalArgumentException("Specify only one of --querySmiles or --queryIdcode");
        }
        StereoMolecule molecule = new StereoMolecule();
        try {
            if (smiles != null) {
                new SmilesParser().parse(molecule, smiles);
            } else {
                new IDCodeParser().parse(molecule, idcode);
            }
        } catch (Exception e) {
            throw new IllegalArgumentException("Unable to parse query molecule", e);
        }
        molecule.ensureHelperArrays(StereoMolecule.cHelperCIP);
        return molecule;
    }

    private static LocalOptimizationLogLevel parseLogLevel(String value) {
        try {
            return LocalOptimizationLogLevel.valueOf(value.toUpperCase());
        } catch (Exception e) {
            return LocalOptimizationLogLevel.NONE;
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("seeds").required(true).hasArg()
                .desc("Input TSV produced by seed finder").build());
        options.addOption(Option.builder().longOpt("output").required(true).hasArg()
                .desc("Output TSV path for optimized assemblies").build());
        options.addOption(Option.builder().longOpt("space").hasArg()
                .desc("Path to full SynthonSpace .data file").build());
        options.addOption(Option.builder().longOpt("downsampled").hasArg()
                .desc("Path to DownsampledSynthonSpace .gz file").build());
        options.addOption(Option.builder().longOpt("raw").hasArg()
                .desc("Path to RawSynthonSpace JSON file").build());
        options.addOption(Option.builder().longOpt("rawMode").hasArg()
                .desc("When using --raw choose FULL or DOWNSAMPLED view (default FULL)").build());
        options.addOption(Option.builder().longOpt("querySmiles").hasArg()
                .desc("Query ligand as SMILES").build());
        options.addOption(Option.builder().longOpt("queryIdcode").hasArg()
                .desc("Query ligand as IDCode").build());
        options.addOption(Option.builder().longOpt("beam").hasArg()
                .desc("Beam size").build());
        options.addOption(Option.builder().longOpt("topL").hasArg()
                .desc("Neighbor pool size L").build());
        options.addOption(Option.builder().longOpt("sampleNeighbors").hasArg()
                .desc("Sampled neighbors per expansion M").build());
        options.addOption(Option.builder().longOpt("cap").hasArg()
                .desc("Per-position diversity cap r").build());
        options.addOption(Option.builder().longOpt("rounds").hasArg()
                .desc("Maximum outer loops T").build());
        options.addOption(Option.builder().longOpt("patience").hasArg()
                .desc("Early stopping patience S").build());
        options.addOption(Option.builder().longOpt("minPhesa").hasArg()
                .desc("Minimum PheSA similarity").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Random seed").build());
        options.addOption(Option.builder().longOpt("log").hasArg()
                .desc("Log level: NONE, IMPROVEMENTS, VERBOSE").build());
        options.addOption(Option.builder().longOpt("seedThreads").hasArg()
                .desc("Number of parallel seed optimizations").build());
        return options;
    }
}
