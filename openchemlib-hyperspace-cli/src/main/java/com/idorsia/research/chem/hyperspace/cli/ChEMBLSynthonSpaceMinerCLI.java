package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import com.idorsia.research.chem.hyperspace.tools.chembl.ObservedSynthonSpaceMiner;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.nio.file.Path;

/**
 * CLI for mining observed-pattern pseudo-reaction synthon spaces from ChEMBL molecule tables.
 */
public final class ChEMBLSynthonSpaceMinerCLI {

    private ChEMBLSynthonSpaceMinerCLI() {
    }

    public static void main(String[] args) throws Exception {
        Options cliOptions = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(cliOptions, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("ChEMBLSynthonSpaceMinerCLI", cliOptions, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        ObservedSynthonSpaceMiner.Options options = ObservedSynthonSpaceMiner.Options.builder()
                .input(Path.of(cmd.getOptionValue("input")))
                .spaceName(cmd.getOptionValue("spaceName", "chembl_observed_synthon_space"))
                .idColumn(cmd.getOptionValue("idColumn", "chembl_id"))
                .structureColumn(cmd.getOptionValue("structureColumn", "idcode"))
                .maxMolecules(Integer.parseInt(cmd.getOptionValue("maxMolecules", "0")))
                .maxCuts(Integer.parseInt(cmd.getOptionValue("maxCuts", "2")))
                .topReactions(Integer.parseInt(cmd.getOptionValue("topReactions", "20")))
                .minSetSize(Integer.parseInt(cmd.getOptionValue("minSetSize", "100")))
                .minSourceMolecules(Integer.parseInt(cmd.getOptionValue("minSourceMolecules", "100")))
                .maxFragmentHeavyAtoms(Integer.parseInt(cmd.getOptionValue("maxFragmentHeavyAtoms", "40")))
                .maxCutsetsPerMolecule(Integer.parseInt(cmd.getOptionValue("maxCutsetsPerMolecule", "500")))
                .sampleProducts(Integer.parseInt(cmd.getOptionValue("sampleProducts", "25")))
                .seed(Long.parseLong(cmd.getOptionValue("seed", "7")))
                .allowAdjacentTwoCuts(cmd.hasOption("allowAdjacentTwoCuts"))
                .splitFilter(ObservedSynthonSpaceMiner.SplitFilter.parse(cmd.getOptionValue("splitFilter", "both")))
                .descriptorShortName(cmd.getOptionValue("descriptor", "FragFp"))
                .descriptorBits(Integer.parseInt(cmd.getOptionValue("descriptorBits", "1024")))
                .build();

        ObservedSynthonSpaceMiner.Result result = new ObservedSynthonSpaceMiner().mine(options);
        Path output = Path.of(cmd.getOptionValue("output"));
        RawSynthonSpaceIO.write(result.rawSpace(), output);

        System.out.println("Raw observed synthon space written to: " + output);
        System.out.println("Input rows: " + result.stats().inputRows());
        System.out.println("Processed molecules: " + result.stats().processedMolecules());
        System.out.println("Accepted splits: " + result.stats().acceptedSplits());
        System.out.println("Rejected splits: " + result.stats().rejectedSplits());
        System.out.println("Selected reactions: " + result.reactions().size());
        for (ObservedSynthonSpaceMiner.Summary summary : result.reactions()) {
            System.out.println(summary.splitCount() + "-split "
                    + summary.synthonCounts()
                    + " products~" + summary.estimatedProductCount()
                    + " support=" + summary.sourceMoleculeCount());
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("input").hasArg().required(true)
                .desc("Input ChEMBL-style TSV file; .gz is supported").build());
        options.addOption(Option.builder().longOpt("output").hasArg().required(true)
                .desc("Output RawSynthonSpace JSON file; .gz enables gzip").build());
        options.addOption(Option.builder().longOpt("spaceName").hasArg()
                .desc("Logical name for the generated raw space").build());
        options.addOption(Option.builder().longOpt("idColumn").hasArg()
                .desc("Input column containing molecule IDs (default chembl_id)").build());
        options.addOption(Option.builder().longOpt("structureColumn").hasArg()
                .desc("Input column containing OCL idcodes (default idcode)").build());
        options.addOption(Option.builder().longOpt("maxMolecules").hasArg()
                .desc("Maximum molecules to process; 0 means all (default 0)").build());
        options.addOption(Option.builder().longOpt("maxCuts").hasArg()
                .desc("Maximum cut count: 1 or 2 (default 2)").build());
        options.addOption(Option.builder().longOpt("topReactions").hasArg()
                .desc("Number of ranked pseudo-reactions to export (default 20)").build());
        options.addOption(Option.builder().longOpt("minSetSize").hasArg()
                .desc("Minimum unique synthons per exported fragment set (default 100)").build());
        options.addOption(Option.builder().longOpt("minSourceMolecules").hasArg()
                .desc("Minimum source molecule support per exported reaction (default 100)").build());
        options.addOption(Option.builder().longOpt("maxFragmentHeavyAtoms").hasArg()
                .desc("Maximum non-connector heavy atoms per synthon (default 40)").build());
        options.addOption(Option.builder().longOpt("maxCutsetsPerMolecule").hasArg()
                .desc("Cap accepted cutsets per source molecule (default 500)").build());
        options.addOption(Option.builder().longOpt("sampleProducts").hasArg()
                .desc("Number of random assembled products to store per reaction (default 25)").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Sampling seed (default 7)").build());
        options.addOption(Option.builder().longOpt("allowAdjacentTwoCuts")
                .desc("Allow two cuts sharing an atom; disabled by default").build());
        options.addOption(Option.builder().longOpt("splitFilter").hasArg()
                .desc("Restrict exported splits: 1, 2, or both (default both)").build());
        options.addOption(Option.builder().longOpt("descriptor").hasArg()
                .desc("Descriptor short name metadata for later build (default FragFp)").build());
        options.addOption(Option.builder().longOpt("descriptorBits").hasArg()
                .desc("Descriptor bit count metadata for later build (default 1024)").build());
        return options;
    }
}
