package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.HydrogenAssembler;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.DescriptorHandlerShapeOneConf;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.actelion.research.chem.io.SDFileParser;
import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.localopt.LocalOptimizationRequest;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.screening.CandidateSampler;
import com.idorsia.research.chem.hyperspace.screening.ContinuousScreeningOrchestrator;
import com.idorsia.research.chem.hyperspace.screening.ScreeningMetrics;

import java.nio.file.Path;

public final class ContinuousScreeningCLI {

    public static void main(String[] args) throws Exception {
        Path configPath = parseConfigPath(args);
        ContinuousScreeningConfig cliConfig = ContinuousScreeningConfig.fromJson(configPath);
        RawSynthonSpace fullRaw = HyperspaceIOUtils.loadRawSynthonSpace(cliConfig.resolveRawFull(configPath).toString());
        RawSynthonSpace downsampledRaw =
                HyperspaceIOUtils.loadRawSynthonSpace(cliConfig.resolveRawDownsampled(configPath).toString());
        PheSAMolecule queryDescriptor = buildQueryDescriptor(cliConfig.resolveQuery(configPath));
        CandidateSampler.Config samplerConfig = cliConfig.toSamplerConfig();
        LocalOptimizationRequest fullRequest = cliConfig.toFullOptimizationRequest();
        boolean microEnabled = cliConfig.getMicroOptimization().isEnabled();
        LocalOptimizationRequest microRequest = microEnabled ? cliConfig.toMicroOptimizationRequest() : null;

        ContinuousScreeningOrchestrator.Config config = new ContinuousScreeningOrchestrator.Config()
                .withFullRaw(fullRaw)
                .withDownsampledRaw(downsampledRaw)
                .withQueryDescriptor(queryDescriptor)
                .withSamplerConfig(samplerConfig)
                .withFullOptimizationRequest(fullRequest)
                .withFullOptimizerThreads(cliConfig.getOrchestration().getWorkerThreads())
                .withQueueCapacity(cliConfig.getOrchestration().getQueueCapacity())
                .withProgressIntervalSeconds(cliConfig.getOrchestration().getProgressIntervalSeconds())
                .withReactionWeightExponent(cliConfig.getOrchestration().getReactionWeightExponent())
                .withReactionMinWeight(cliConfig.getOrchestration().getReactionMinWeight())
                .withHitOutput(cliConfig.resolveOutputHits(configPath))
                .withMinReportedSimilarity(cliConfig.getMinReportedSimilarity())
                .withDuplicateCacheSize(cliConfig.getOrchestration().getDuplicateCacheSize())
                .withRandomSeed(cliConfig.getRun().getRandomSeed());
        if (microEnabled) {
            config.withMicroEnabled(true).withMicroOptimizationRequest(microRequest);
        }
        config.validate();

        long iterations = cliConfig.getRun().getIterations();
        try (ContinuousScreeningOrchestrator orchestrator = new ContinuousScreeningOrchestrator(config)) {
            orchestrator.run(iterations);
            ScreeningMetrics metrics = orchestrator.getMetrics();
            System.out.println("Sampled evaluations: " + metrics.getSampled());
            System.out.println("Duplicate seeds skipped: " + metrics.getDuplicateSeeds());
            System.out.println("Submitted to optimizer: " + metrics.getSubmitted());
            System.out.println("Total hits: " + metrics.getHits());
        }
    }

    private static Path parseConfigPath(String[] args) {
        if (args.length == 1 && !args[0].startsWith("--")) {
            return Path.of(args[0]);
        }
        if (args.length == 2 && "--config".equals(args[0])) {
            return Path.of(args[1]);
        }
        if (args.length == 1 && ("--help".equals(args[0]) || "-h".equals(args[0]))) {
            printUsage();
            System.exit(0);
        }
        printUsage();
        throw new IllegalArgumentException("Provide a JSON config file path");
    }

    private static void printUsage() {
        System.err.println("Usage: ContinuousScreeningCLI <config.json>");
        System.err.println("   or: ContinuousScreeningCLI --config <config.json>");
    }

    private static PheSAMolecule buildQueryDescriptor(ContinuousScreeningConfig.QueryInput queryInput) throws Exception {
        if (hasText(queryInput.getSdfFile())) {
            return buildSdfQueryDescriptor(queryInput);
        }
        StereoMolecule molecule = parseQuery(queryInput);
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

    private static PheSAMolecule buildSdfQueryDescriptor(ContinuousScreeningConfig.QueryInput queryInput) {
        int recordIndex = queryInput.getSdfRecordIndex();
        Path sdfPath = Path.of(queryInput.getSdfFile());

        SDFileParser parser = new SDFileParser(sdfPath.toString());
        if (!parser.isOpen()) {
            throw new IllegalArgumentException("Unable to open query.sdfFile: " + sdfPath);
        }

        try {
            StereoMolecule molecule = null;
            int currentIndex = -1;
            while (parser.next()) {
                currentIndex++;
                if (currentIndex == recordIndex) {
                    molecule = parser.getMolecule();
                    if (molecule == null) {
                        throw new IllegalArgumentException("Unable to parse query molecule at sdf record "
                                + recordIndex + " from file: " + sdfPath);
                    }
                    break;
                }
            }

            if (molecule == null) {
                throw new IllegalArgumentException("query.sdfRecordIndex " + recordIndex
                        + " not found in sdfFile: " + sdfPath);
            }

            StereoMolecule prepared = new StereoMolecule(molecule);
            new HydrogenAssembler(prepared).addImplicitHydrogens();
            prepared.ensureHelperArrays(Molecule.cHelperCIP);

            DescriptorHandlerShapeOneConf handler = new DescriptorHandlerShapeOneConf();
            PheSAMolecule descriptor = handler.createDescriptor(prepared);
            if (handler.calculationFailed(descriptor)) {
                throw new IllegalStateException("Unable to create single-conformation PheSA descriptor from sdfFile: "
                        + sdfPath + " record " + recordIndex);
            }
            return descriptor;
        } finally {
            parser.close();
        }
    }

    private static StereoMolecule parseQuery(ContinuousScreeningConfig.QueryInput queryInput) throws Exception {
        String smiles = queryInput.getSmiles();
        String idcode = queryInput.getIdcode();
        boolean hasSmiles = smiles != null && !smiles.isBlank();
        boolean hasIdcode = idcode != null && !idcode.isBlank();
        if (hasSmiles == hasIdcode) {
            throw new IllegalArgumentException("Specify exactly one of query.smiles or query.idcode");
        }
        StereoMolecule molecule = new StereoMolecule();
        if (hasSmiles) {
            new SmilesParser().parse(molecule, smiles);
        } else {
            new IDCodeParser().parse(molecule, idcode);
        }
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        return molecule;
    }

    private static boolean hasText(String value) {
        return value != null && !value.isBlank();
    }
}
