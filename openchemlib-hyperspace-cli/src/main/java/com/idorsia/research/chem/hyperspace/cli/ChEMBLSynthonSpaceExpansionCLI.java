package com.idorsia.research.chem.hyperspace.cli;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import com.idorsia.research.chem.hyperspace.tools.chembl.ChemblSynthonSpaceExpander;
import com.idorsia.research.chem.hyperspace.tools.chembl.ObservedSynthonSpaceMiner;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.nio.file.Path;

/**
 * CLI for expanding observed ChEMBL pseudo-reaction rawspaces with compatible mined synthons.
 */
public final class ChEMBLSynthonSpaceExpansionCLI {

    private ChEMBLSynthonSpaceExpansionCLI() {
    }

    public static void main(String[] args) throws Exception {
        Options cliOptions = buildOptions();
        CommandLine cmd;
        try {
            cmd = new DefaultParser().parse(cliOptions, args);
        } catch (ParseException e) {
            new HelpFormatter().printHelp("ChEMBLSynthonSpaceExpansionCLI", cliOptions, true);
            throw new IllegalArgumentException("Unable to parse arguments", e);
        }

        RawSynthonSpace seed = RawSynthonSpaceIO.read(Path.of(cmd.getOptionValue("rawIn")));
        ChemblSynthonSpaceExpander.Options options = ChemblSynthonSpaceExpander.Options.builder()
                .chemblInput(Path.of(cmd.getOptionValue("chemblInput")))
                .idColumn(cmd.getOptionValue("idColumn", "chembl_id"))
                .structureColumn(cmd.getOptionValue("structureColumn", "idcode"))
                .catalogMaxMolecules(Integer.parseInt(cmd.getOptionValue("catalogMaxMolecules", "10000")))
                .catalogSplitFilter(ObservedSynthonSpaceMiner.SplitFilter.parse(
                        cmd.getOptionValue("catalogSplitFilter", "both")))
                .targetConnectorCount(ChemblSynthonSpaceExpander.TargetConnectorCount.parse(
                        cmd.getOptionValue("targetConnectorCount", "both")))
                .maxCutsetsPerMolecule(Integer.parseInt(cmd.getOptionValue("maxCutsetsPerMolecule", "500")))
                .maxAddedPerSet(Integer.parseInt(cmd.getOptionValue("maxAddedPerSet", "1000")))
                .maxCandidateHeavyAtoms(Integer.parseInt(cmd.getOptionValue("maxCandidateHeavyAtoms", "24")))
                .minFragFpSimilarity(Double.parseDouble(cmd.getOptionValue("minFragFpSimilarity", "0.25")))
                .connectorRegionSize(Integer.parseInt(cmd.getOptionValue("connectorRegionSize", "3")))
                .seed(Long.parseLong(cmd.getOptionValue("seed", "7")))
                .allowAdjacentTwoCuts(cmd.hasOption("allowAdjacentTwoCuts"))
                .mode(ChemblSynthonSpaceExpander.ExpansionMode.parse(cmd.getOptionValue("mode", "strict")))
                .descriptorShortName(cmd.getOptionValue("descriptor", "FragFp"))
                .similarityWeight(Double.parseDouble(cmd.getOptionValue("similarityWeight", "0.55")))
                .smallnessWeight(Double.parseDouble(cmd.getOptionValue("smallnessWeight", "0.35")))
                .occurrenceWeight(Double.parseDouble(cmd.getOptionValue("occurrenceWeight", "0.10")))
                .sizeBiasPower(Double.parseDouble(cmd.getOptionValue("sizeBiasPower", "1.0")))
                .progressInterval(Integer.parseInt(cmd.getOptionValue("progressInterval", "1000")))
                .build();

        ChemblSynthonSpaceExpander.Result result = new ChemblSynthonSpaceExpander().expand(seed, options);
        Path output = Path.of(cmd.getOptionValue("rawOut"));
        RawSynthonSpaceIO.write(result.rawSpace(), output);

        System.out.println("Expanded raw synthon space written to: " + output);
        System.out.println("Catalog input rows: " + result.catalogStats().inputRows());
        System.out.println("Catalog processed molecules: " + result.catalogStats().processedMolecules());
        System.out.println("Catalog accepted splits: " + result.catalogStats().acceptedSplits());
        System.out.println("Catalog rejected splits: " + result.catalogStats().rejectedSplits());
        System.out.println("Catalog fragments: " + result.catalogStats().catalogFragments());
        System.out.println("Catalog 1-connector fragments: " + result.catalogStats().oneConnectorFragments());
        System.out.println("Catalog 2-connector fragments: " + result.catalogStats().twoConnectorFragments());
        for (ChemblSynthonSpaceExpander.ReactionSummary summary : result.reactions()) {
            System.out.println(summary.reactionId()
                    + " " + summary.beforeCounts()
                    + " -> " + summary.afterCounts()
                    + " added=" + summary.addedSynthons()
                    + " products~" + summary.estimatedProductCount());
        }
    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption(Option.builder().longOpt("rawIn").hasArg().required(true)
                .desc("Input seed RawSynthonSpace JSON file; .gz is supported").build());
        options.addOption(Option.builder().longOpt("chemblInput").hasArg().required(true)
                .desc("ChEMBL-style TSV file used to mine expansion synthons; .gz is supported").build());
        options.addOption(Option.builder().longOpt("rawOut").hasArg().required(true)
                .desc("Output expanded RawSynthonSpace JSON file; .gz enables gzip").build());
        options.addOption(Option.builder().longOpt("idColumn").hasArg()
                .desc("Input column containing molecule IDs (default chembl_id)").build());
        options.addOption(Option.builder().longOpt("structureColumn").hasArg()
                .desc("Input column containing OCL idcodes (default idcode)").build());
        options.addOption(Option.builder().longOpt("catalogMaxMolecules").hasArg()
                .desc("Maximum catalog molecules to process; 0 means all (default 10000)").build());
        options.addOption(Option.builder().longOpt("catalogSplitFilter").hasArg()
                .desc("Catalog split mining: 1, 2, or both (default both)").build());
        options.addOption(Option.builder().longOpt("targetConnectorCount").hasArg()
                .desc("Synthon sets to expand by connector count: 1, 2, or both (default both)").build());
        options.addOption(Option.builder().longOpt("maxCutsetsPerMolecule").hasArg()
                .desc("Cap accepted catalog cutsets per source molecule (default 500)").build());
        options.addOption(Option.builder().longOpt("maxAddedPerSet").hasArg()
                .desc("Maximum expansion synthons added per synthon set (default 1000)").build());
        options.addOption(Option.builder().longOpt("maxCandidateHeavyAtoms").hasArg()
                .desc("Maximum non-connector heavy atoms per added synthon (default 24)").build());
        options.addOption(Option.builder().longOpt("minFragFpSimilarity").hasArg()
                .desc("Minimum FragFp Tanimoto to a seed-set anchor (default 0.25)").build());
        options.addOption(Option.builder().longOpt("connectorRegionSize").hasArg()
                .desc("Connector-proximal region size metadata/feature radius (default 3)").build());
        options.addOption(Option.builder().longOpt("seed").hasArg()
                .desc("Sampling/ranking seed (default 7)").build());
        options.addOption(Option.builder().longOpt("allowAdjacentTwoCuts")
                .desc("Allow two catalog cuts sharing an atom; disabled by default").build());
        options.addOption(Option.builder().longOpt("mode").hasArg()
                .desc("Expansion mode: strict (default strict)").build());
        options.addOption(Option.builder().longOpt("descriptor").hasArg()
                .desc("Descriptor short name for similarity ranking (default FragFp)").build());
        options.addOption(Option.builder().longOpt("similarityWeight").hasArg()
                .desc("Similarity score weight; normalized with other weights (default 0.55)").build());
        options.addOption(Option.builder().longOpt("smallnessWeight").hasArg()
                .desc("Small-fragment score weight; normalized with other weights (default 0.35)").build());
        options.addOption(Option.builder().longOpt("occurrenceWeight").hasArg()
                .desc("Occurrence score weight; normalized with other weights (default 0.10)").build());
        options.addOption(Option.builder().longOpt("sizeBiasPower").hasArg()
                .desc("Power applied to the smallness score; higher values favor smaller synthons (default 1.0)").build());
        options.addOption(Option.builder().longOpt("progressInterval").hasArg()
                .desc("Catalog progress interval in processed molecules; 0 disables progress output (default 1000)").build());
        return options;
    }
}
