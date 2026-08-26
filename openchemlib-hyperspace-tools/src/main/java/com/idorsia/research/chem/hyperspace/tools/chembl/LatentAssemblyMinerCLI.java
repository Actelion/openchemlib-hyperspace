package com.idorsia.research.chem.hyperspace.tools.chembl;

import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

/** Command-line entry point for the latent assembly ChEMBL miner. */
public final class LatentAssemblyMinerCLI {
    private LatentAssemblyMinerCLI() {}

    public static void main(String[] args) throws Exception {
        if (args.length == 0 || hasFlag(args, "--help")) {
            usage();
            return;
        }
        Map<String, String> values = parseArgs(args);
        Path input = Path.of(required(values, "--input"));
        Path output = Path.of(required(values, "--output"));

        LatentAssemblyMiner.Builder builder = LatentAssemblyMiner.Options.builder()
                .input(input)
                .idColumn(values.getOrDefault("--id-column", "chembl_id"))
                .structureColumn(values.getOrDefault("--structure-column", "idcode"))
                .maxMolecules(integer(values, "--max-molecules", 0))
                .minProductNonHydrogenAtoms(integer(values, "--min-product-atoms", 12))
                .maxProductNonHydrogenAtoms(integer(values, "--max-product-atoms", 32))
                .minScaffoldNonHydrogenAtoms(integer(values, "--min-scaffold-atoms", 6))
                .maxScaffoldNonHydrogenAtoms(integer(values, "--max-scaffold-atoms", 24))
                .minScaffoldRingAtoms(integer(values, "--min-scaffold-ring-atoms", 5))
                .minScaffoldFraction(decimal(values, "--min-scaffold-fraction", 0.30))
                .maxScaffoldFraction(decimal(values, "--max-scaffold-fraction", 0.75))
                .minArmNonHydrogenAtoms(integer(values, "--min-arm-atoms", 3))
                .maxArmNonHydrogenAtoms(integer(values, "--max-arm-atoms", 12))
                .minAttachmentDistance(integer(values, "--min-attachment-distance", 2))
                .maxDecompositionsPerMolecule(integer(values, "--max-per-molecule", 1))
                .onlyCommonDrugElements(!values.containsKey("--include-uncommon-elements"))
                .armMode(LatentAssemblyMiner.ArmMode.parse(values.getOrDefault("--arms", "2")));

        LatentAssemblyMiner.Stats stats;
        try (LatentAssemblyJsonlWriter writer = new LatentAssemblyJsonlWriter(output)) {
            stats = new LatentAssemblyMiner().mine(builder.build(), writer::write);
        }
        printStats(stats, output);
    }

    private static Map<String, String> parseArgs(String[] args) {
        Map<String, String> values = new LinkedHashMap<>();
        for (int i = 0; i < args.length; i++) {
            String option = args[i];
            if (!option.startsWith("--")) {
                throw new IllegalArgumentException("Expected an option, got: " + option);
            }
            if ("--include-uncommon-elements".equals(option)) {
                values.put(option, "true");
                continue;
            }
            if (i + 1 >= args.length || args[i + 1].startsWith("--")) {
                throw new IllegalArgumentException("Missing value for option: " + option);
            }
            values.put(option, args[++i]);
        }
        return values;
    }

    private static String required(Map<String, String> values, String option) {
        String value = values.get(option);
        if (value == null || value.isBlank()) {
            throw new IllegalArgumentException("Missing required option: " + option);
        }
        return value;
    }

    private static int integer(Map<String, String> values, String option, int defaultValue) {
        return values.containsKey(option) ? Integer.parseInt(values.get(option)) : defaultValue;
    }

    private static double decimal(Map<String, String> values, String option, double defaultValue) {
        return values.containsKey(option) ? Double.parseDouble(values.get(option)) : defaultValue;
    }

    private static boolean hasFlag(String[] args, String flag) {
        for (String arg : args) {
            if (flag.equals(arg)) return true;
        }
        return false;
    }

    private static void printStats(LatentAssemblyMiner.Stats stats, Path output) {
        System.out.println("Latent assembly mining complete");
        System.out.println("output=" + output.toAbsolutePath());
        System.out.println("input_rows=" + stats.inputRows());
        System.out.println("parsed_molecules=" + stats.parsedMolecules());
        System.out.println("eligible_products=" + stats.eligibleProducts());
        System.out.println("rejected_products_too_large=" + stats.rejectedProductsTooLarge());
        System.out.println("products_with_decompositions=" + stats.productsWithDecompositions());
        System.out.println("accepted_two_arm_candidates=" + stats.acceptedTwoArmCandidates());
        System.out.println("accepted_three_arm_candidates=" + stats.acceptedThreeArmCandidates());
        System.out.println("emitted_decompositions=" + stats.emittedDecompositions());
    }

    private static void usage() {
        System.out.println("Usage: LatentAssemblyMinerCLI --input <chembl.tsv[.gz]> --output <records.jsonl[.gz]> [options]");
        System.out.println("  --arms <2|3|both>                 default: 2");
        System.out.println("  --max-molecules <n>               default: 0 (all)");
        System.out.println("  --max-per-molecule <n>            default: 1");
        System.out.println("  --min-product-atoms <n>            default: 12");
        System.out.println("  --max-product-atoms <n>            default and hard maximum: 32");
        System.out.println("  --min-scaffold-atoms <n>           default: 6");
        System.out.println("  --max-scaffold-atoms <n>           default: 24");
        System.out.println("  --min-scaffold-ring-atoms <n>      default: 5");
        System.out.println("  --min-scaffold-fraction <f>        default: 0.30");
        System.out.println("  --max-scaffold-fraction <f>        default: 0.75");
        System.out.println("  --min-arm-atoms <n>                default: 3");
        System.out.println("  --max-arm-atoms <n>                default: 12");
        System.out.println("  --min-attachment-distance <n>      default: 2");
        System.out.println("  --id-column <name>                 default: chembl_id");
        System.out.println("  --structure-column <name>          default: idcode");
        System.out.println("  --include-uncommon-elements        disable the common drug-element filter");
    }
}
