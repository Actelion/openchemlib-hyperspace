package com.idorsia.research.chem.hyperspace.cli;

import com.actelion.research.chem.StereoMolecule;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;

import java.util.List;

public interface CLIMethod {
    public static interface CLIMethodConfig{};

    public String getMethodName();
    public List<Option> getCommandLineOptions();
    public void parseMethodConfigFromCommandline(CommandLine cli);
    public void runMethod(List<StereoMolecule> structures, String output);

}
