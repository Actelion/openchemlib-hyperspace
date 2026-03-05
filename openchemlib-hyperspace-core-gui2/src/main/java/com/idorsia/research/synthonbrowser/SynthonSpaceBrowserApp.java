package com.idorsia.research.synthonbrowser;

import javax.swing.*;
import java.nio.file.Path;

public class SynthonSpaceBrowserApp {

    public static void main(String[] args) {
        Path initialPath = null;
        if (args != null && args.length > 0 && args[0] != null && !args[0].isBlank()) {
            initialPath = Path.of(args[0]);
        }
        Path finalInitialPath = initialPath;

        SwingUtilities.invokeLater(() -> {
            SynthonSpaceBrowserFrame frame = new SynthonSpaceBrowserFrame();
            frame.setVisible(true);
            if (finalInitialPath != null) {
                frame.openPath(finalInitialPath);
            }
        });
    }
}
