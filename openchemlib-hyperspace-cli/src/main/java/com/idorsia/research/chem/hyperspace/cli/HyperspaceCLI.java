package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.StereoMolecule;
import org.apache.commons.cli.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

//

public class HyperspaceCLI {

    public static Map<String,CLIMethod> methods = new HashMap<>();

    private static String inputFile = "";
    private static String outputFile = "";



    // Required options:
    // "input"
    //
    public static void main(String args[]) {

        System.out.println("Hyperspace CLI 0.1");
        System.out.println("Idorsia SciComp 2022");

        // load CLIMethods:
        initMethod(new HyperspaceSubstructureSearchCLIMethod());

        // init command line parsing:
        Option input = new Option("input", true, "Input file, ending indicates input type. Supported are .idc, .dwar, .smi");
        Option output = new Option("output", true, "Output file");
        Option method = new Option("method", true, "Method to use, supported are: " + String.join(" ", methods.keySet().stream().collect(Collectors.toList())));

        Options options_a = new Options();
        options_a.addOption(input);
        options_a.addOption(output);
        options_a.addOption(method);

        // init method options:
        System.out.println("");
        for (CLIMethod mi : methods.values()) {
            for (Option oi : mi.getCommandLineOptions()) {
                Option oi_2 = new Option(oi.getArgName(), oi.hasArg(), mi.getMethodName() + " -> " + oi.getDescription());
                options_a.addOption(oi_2);
            }
        }

        // parse arguments, and try to parse input
        CommandLineParser cpa = new DefaultParser(true);
        CommandLine pcli      = null;
        String val_method = "";
        try {
            pcli       =  cpa.parse(options_a, args, true);
            inputFile  = pcli.getOptionValue(input);
            outputFile = pcli.getOptionValue(output);
            val_method = pcli.getOptionValue(method);
        } catch (ParseException e) {
            e.printStackTrace();
        }

        // Load input molecules:
        System.out.println("Try to load input structures:");
        List<StereoMolecule> in_structures = new ArrayList<>();
        if (inputFile.toLowerCase().endsWith(".txt") ||
                inputFile.toLowerCase().endsWith(".idc")) {
            in_structures.addAll(InputStructureParser.parseIDC(inputFile));
        }
        System.out.println("Try to load input structures: done -> parsed " + in_structures.size() + " structures");

        // run method:
        if (val_method.equals("sss")) {
            // 1. config
            System.out.println("Configure sss:");
            methods.get("sss").parseMethodConfigFromCommandline(pcli);

            // 2. run:
            methods.get("sss").runMethod(in_structures,outputFile);
        }





    }


    public static void initMethod(CLIMethod method) {
        methods.put(method.getMethodName(),method);
    }

}
