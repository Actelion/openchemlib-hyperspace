package com.idorsia.research.chem.hyperspace.gui.action.stats;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import com.idorsia.research.chem.hyperspace.gui.search.HyperspaceSubstructureSearch;

import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.awt.event.ActionEvent;
import java.io.File;
import java.nio.file.Path;

public final class ExportRawSynthonSpaceStatsAction extends AbstractAction {
    private final HyperspaceSearchGUI gui;

    public ExportRawSynthonSpaceStatsAction(HyperspaceSearchGUI gui) {
        super("Export Rawspace Statistics...");
        this.gui = gui;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        File defaultInput = findDefaultRawspaceFile();
        JFileChooser inputChooser = new JFileChooser();
        inputChooser.setDialogTitle("Select RawSynthonSpace file");
        if (defaultInput != null) {
            inputChooser.setSelectedFile(defaultInput);
        }
        int inputSelection = inputChooser.showOpenDialog(gui.getFrame());
        if (inputSelection != JFileChooser.APPROVE_OPTION) {
            return;
        }

        JFileChooser outputChooser = new JFileChooser();
        outputChooser.setDialogTitle("Select output directory");
        outputChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        File selectedInput = inputChooser.getSelectedFile();
        File parent = selectedInput == null ? null : selectedInput.getParentFile();
        if (parent != null) {
            outputChooser.setSelectedFile(new File(parent, stripRawspaceExtension(selectedInput.getName()) + "_stats"));
        }
        int outputSelection = outputChooser.showSaveDialog(gui.getFrame());
        if (outputSelection != JFileChooser.APPROVE_OPTION) {
            return;
        }

        Path rawInput = selectedInput.toPath().toAbsolutePath().normalize();
        Path outputDirectory = outputChooser.getSelectedFile().toPath().toAbsolutePath().normalize();
        RawSynthonSpaceStatsExportProcess process = new RawSynthonSpaceStatsExportProcess(
                rawInput,
                outputDirectory,
                10,
                100,
                13L,
                Runtime.getRuntime().availableProcessors());
        process.addSearchProviderListener(() -> {
            if (process.getProcessStatus() == RawSynthonSpaceStatsExportProcess.ProcessStatus.DONE) {
                JOptionPane.showMessageDialog(gui.getFrame(),
                        "Rawspace statistics written to:\n" + outputDirectory,
                        "Export complete",
                        JOptionPane.INFORMATION_MESSAGE);
            } else if (process.getProcessStatus() == RawSynthonSpaceStatsExportProcess.ProcessStatus.FAILED) {
                JOptionPane.showMessageDialog(gui.getFrame(),
                        "Rawspace statistics export failed:\n" + process.getProcessStatusMessage(),
                        "Export failed",
                        JOptionPane.ERROR_MESSAGE);
            }
        });
        gui.getProcessListModel().addProcess(process);
        process.startAsync();
    }

    private File findDefaultRawspaceFile() {
        if (gui.getHyperspaceMainPanel() == null
                || gui.getHyperspaceMainPanel().getHyperspaceSearchPanel() == null
                || gui.getHyperspaceMainPanel().getHyperspaceSearchPanel().getSearchProviderListPanel() == null) {
            return null;
        }
        for (AbstractSearchProvider provider : gui.getHyperspaceMainPanel()
                .getHyperspaceSearchPanel()
                .getSearchProviderListPanel()
                .getSearchProviders()) {
            AbstractSearchProvider.SearchProviderConfiguration configuration = provider.getSearchProviderConfiguration();
            if (configuration instanceof HyperspaceSubstructureSearch.InitializationConfig init) {
                String file = init.getFile();
                if (looksLikeRawspace(file)) {
                    return new File(file);
                }
            }
        }
        return null;
    }

    private static boolean looksLikeRawspace(String file) {
        if (file == null) {
            return false;
        }
        String lower = file.toLowerCase();
        return lower.endsWith(".rawspace") || lower.endsWith(".rawspace.gz");
    }

    private static String stripRawspaceExtension(String name) {
        if (name.endsWith(".rawspace.gz")) {
            return name.substring(0, name.length() - ".rawspace.gz".length());
        }
        if (name.endsWith(".rawspace")) {
            return name.substring(0, name.length() - ".rawspace".length());
        }
        return name;
    }
}
