package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.downsampling.DownsampledSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.seed.QueryAwareSeedFinder;
import com.idorsia.research.chem.hyperspace.seed.QueryAwareSeedFinderRequest;
import com.idorsia.research.chem.hyperspace.seed.QueryAwareSeedFinderResult;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.io.IOException;
import java.nio.file.Path;

/**
 * CLI entry point that executes the query-aware seed finder.
 */
public class SynthonSeedFinderCLI {

    public static void main(String[] args) throws IOException {
        Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("SynthonSeedFinderCLI", options, true);
            throw new IllegalArgumentException("Unable to parse CLI parameters", e);
        }

        String downsampledArg = cmd.getOptionValue("downsampled");
        String rawArg = cmd.getOptionValue("raw");
        Path outputPath = Path.of(cmd.getOptionValue("output"));
        long attempts = Long.parseLong(cmd.getOptionValue("attempts", "1000"));
        int minAtoms = Integer.parseInt(cmd.getOptionValue("minAtoms", "1"));
        int maxAtoms = Integer.parseInt(cmd.getOptionValue("maxAtoms", "0"));
        int maxRotBonds = Integer.parseInt(cmd.getOptionValue("maxRotBonds", Integer.toString(Integer.MAX_VALUE)));
        double minPhesa = Double.parseDouble(cmd.getOptionValue("minPhesa", "0.7"));
        long seed = Long.parseLong(cmd.getOptionValue("seed", "13"));
        int maxHitsPerFile = Integer.parseInt(cmd.getOptionValue("maxHitsPerFile", "0"));
        int threads = Integer.parseInt(cmd.getOptionValue("threads", Integer.toString(Runtime.getRuntime().availableProcessors())));

        StereoMolecule queryMolecule = buildQueryMolecule(cmd);
        DescriptorHandlerShape descriptorHandlerShape = new DescriptorHandlerShape();
        ConformerSetGenerator generator = new ConformerSetGenerator(1);
        ConformerSet queryConformers = generator.generateConformerSet(queryMolecule);
        if (queryConformers.isEmpty()) {
            throw new IllegalStateException("Failed to generate conformer for query molecule");
        }
        PheSAMolecule queryDescriptor = descriptorHandlerShape.createDescriptor(queryConformers);
        if (descriptorHandlerShape.calculationFailed(queryDescriptor)) {
            throw new IllegalStateException("Unable to generate PheSA descriptor for query molecule");
        }

        DownsampledSynthonSpace downsampledSpace;
        if (rawArg != null) {
            if (downsampledArg != null) {
                throw new IllegalArgumentException("Specify only one of --downsampled or --raw");
            }
            RawSynthonSpace rawSpace = HyperspaceIOUtils.loadRawSynthonSpace(rawArg);
            downsampledSpace = rawSpace.toDownsampledSpace();
        } else {
            if (downsampledArg == null) {
                throw new IllegalArgumentException("Provide either --downsampled or --raw");
            }
            downsampledSpace = HyperspaceIOUtils.loadDownsampledSynthonSpace(downsampledArg);
        }
        QueryAwareSeedFinder finder = new QueryAwareSeedFinder(downsampledSpace);
        QueryAwareSeedFinderRequest request = QueryAwareSeedFinderRequest.builder(outputPath)
                .attemptsPerReaction(attempts)
                .minAtoms(minAtoms)
                .maxAtoms(maxAtoms)
                .maxRotatableBonds(maxRotBonds)
                .minPhesaSimilarity(minPhesa)
                .randomSeed(seed)
                .maxHitsPerFile(maxHitsPerFile)
                .numThreads(threads)
                .build();

        QueryAwareSeedFinderResult result = finder.run(queryDescriptor, request);
        System.out.println("Total attempts: " + result.getTotalAttempts());
        System.out.println("Assemblies: " + result.getTotalAssemblies());
        System.out.println("Threshold passes: " + result.getTotalThresholdPasses());
        System.out.println("PheSA evaluations: " + result.getTotalPhesaEvaluations());
        System.out.println("Hits: " + result.getTotalHits());
    }

    private static StereoMolecule buildQueryMolecule(CommandLine cmd) {
        String smiles = cmd.getOptionValue("querySmiles");
        String idcode = cmd.getOptionValue("queryIdcode");
        if (smiles == null && idcode == null) {
            throw new IllegalArgumentException("Provide either --querySmiles or --queryIdcode");
        }
        if (smiles != null && idcode != null) {
            throw new IllegalArgumentException("Please specify only one of --querySmiles or --queryIdcode");
        }
        StereoMolecule mol = new StereoMolecule();
        try {
            if (smiles != null) {
                new SmilesParser().parse(mol, smiles);
            } else {
                new IDCodeParser().parse(mol, idcode);
            }
        } catch (Exception e) {
            throw new IllegalArgumentException("Unable to parse query molecule", e);
        }
        mol.ensureHelperArrays(StereoMolecule.cHelperCIP);
        return mol;
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("downsampled").hasArg()
                .desc("Path to the downsampled synthon space file").build());
        options.addOption(Option.builder().longOpt("raw").hasArg()
                .desc("Path to a raw synthon space JSON that includes downsampled sets").build());
        options.addOption(Option.builder().longOpt("output").hasArg().required(true)
                .desc("Output TSV path for hits").build());
        options.addOption(Option.builder().longOpt("querySmiles").hasArg()
                .desc("Query molecule as SMILES").build());
        options.addOption(Option.builder().longOpt("queryIdcode").hasArg()
                .desc("Query molecule as OpenChemLib IDCode").build());
        options.addOption(Option.builder().longOpt("attempts").hasArg()
                .desc("Attempts per reaction").build());
        options.addOption(Option.builder().longOpt("minAtoms").hasArg()
                .desc("Minimum atoms allowed in assembled molecule").build());
        options.addOption(Option.builder().longOpt("maxAtoms").hasArg()
                .desc("Maximum atoms allowed (0 = unlimited)").build());
        options.addOption(Option.builder().longOpt("maxRotBonds").hasArg()
                .desc("Maximum allowed rotatable bonds").build());
        options.addOption(Option.builder().longOpt("minPhesa").hasArg()
                .desc("Minimum PheSA similarity to keep a hit").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Random seed").build());
        options.addOption(Option.builder().longOpt("maxHitsPerFile").hasArg()
                .desc("Maximum hits per TSV before rotating files (0 = unlimited)").build());
        options.addOption(Option.builder().longOpt("threads").hasArg()
                .desc("Number of worker threads").build());
        return options;
    }
}
