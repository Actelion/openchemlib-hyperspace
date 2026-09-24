package com.idorsia.research.chem.hyperspace.gui2;

import com.formdev.flatlaf.FlatLightLaf;
import com.idorsia.research.chem.hyperspace.gui2.action.LoadSynthonSpaceAction;
import com.idorsia.research.chem.hyperspace.gui2.action.RunSubstructureSearchAction;
import com.idorsia.research.chem.hyperspace.gui2.model.LeetHyperspaceModel;
import com.idorsia.research.chem.hyperspace.gui2.task.LoadSynthonSpaceTask;
import com.idorsia.research.chem.hyperspace.gui2.view.LeetHyperspaceView;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class LeetHyperspaceMain {

    private static JFrame fi;

    private static LeetHyperspaceModel model;
    private static LeetHyperspaceView  view;


    private static void initLAF() {
        try {
            UIManager.setLookAndFeel( new FlatLightLaf() );
        } catch( Exception ex ) {
            System.err.println( "Failed to initialize LaF" );
        }
    }

    public static void main(String[] args) {
        java.nio.file.Path config = args.length > 0 ? java.nio.file.Path.of(args[0])
                : java.nio.file.Path.of("spaces.conf");
        final Gui2Configuration.Configuration configuration;
        try {
            configuration = args.length == 0 && !java.nio.file.Files.exists(config)
                    ? null : Gui2Configuration.read(config);
        } catch (Exception ex) {
            System.err.println("GUI 2 configuration error: " + ex.getMessage());
            if (!GraphicsEnvironment.isHeadless()) SwingUtilities.invokeLater(() ->
                    JOptionPane.showMessageDialog(null, ex.getMessage(), "Cannot open GUI 2",
                            JOptionPane.ERROR_MESSAGE));
            throw new IllegalArgumentException("Cannot load GUI 2 configuration", ex);
        }
        SwingUtilities.invokeLater(() -> open(configuration));
    }

    public static JFrame open(Gui2Configuration.Configuration configuration) {
        if (!SwingUtilities.isEventDispatchThread()) throw new IllegalStateException("Use the Swing EDT");
        initLAF();
        model = new LeetHyperspaceModel();
        view = new LeetHyperspaceView(model);
        fi = new JFrame("Hyperspace - GUI 2 - Substructure Search");
        fi.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        fi.getContentPane().add(view, BorderLayout.CENTER);
        initMenu(fi);
        Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
        fi.setSize(Math.min(1200, screen.width), Math.min(850, screen.height));
        fi.setLocationByPlatform(true);
        fi.setVisible(true);
        if (configuration != null) {
            for (Gui2Configuration.Space space : configuration.spaces()) {
                LoadSynthonSpaceTask task = new LoadSynthonSpaceTask(
                        model, space.file().toString(), space.name(), space.threads());
                model.addTask(task);
                task.execute();
            }
            if (!configuration.warnings().isEmpty()) JOptionPane.showMessageDialog(fi,
                    String.join("\n", configuration.warnings()), "GUI 2: substructure only",
                    JOptionPane.INFORMATION_MESSAGE);
        }
        return fi;
    }

    public static void initMenu(JFrame fi) {
        JMenuBar jb = new JMenuBar();

        JMenu jmFile = new JMenu("File");
        jmFile.add(new LoadSynthonSpaceAction(model,view));

        jb.add(jmFile);

        JMenu jmSearch = new JMenu("Search");
        List<JMenuItem> jmSearchItems = new ArrayList<>();
        jb.add(jmSearch);

        model.addListener(new LeetHyperspaceModel.LeetHyperspaceModelListener() {
            @Override
            public void synthonSpacesChanged() {
                jmSearchItems.stream().forEach( xi -> jmSearch.remove(xi));
                jmSearchItems.clear();
                //model.getSynthonSpaces().stream().forEach( xi -> jmSearchItems.add(new JMenuItem( new RunSubstructureSearchAction(model,view,xi,() -> model.getQuery()))));
                model.getSynthonSpaces().stream().forEach( xi -> jmSearchItems.add(new JMenuItem( new RunSubstructureSearchAction(model,view,xi.getSpace(),xi.getName(),() -> model.getQuery(),xi.getThreads()))));
                jmSearchItems.stream().forEach( xi -> jmSearch.add(xi));
            }
        });

        fi.setJMenuBar(jb);
    }

}
