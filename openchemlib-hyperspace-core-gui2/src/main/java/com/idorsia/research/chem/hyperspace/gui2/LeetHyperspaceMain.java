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

    public static void main(String args[]) {

        initLAF();

        model = new LeetHyperspaceModel();
        view  = new LeetHyperspaceView(model);

        fi = new JFrame("Hyperspace 3.0");
        fi.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

        fi.getContentPane().setLayout(new BorderLayout());
        fi.getContentPane().add(view,BorderLayout.CENTER);
        initMenu(fi);

        fi.setSize(600,400);
        fi.setVisible(true);


        parseConfig("spaces.conf", model);

        //LoadSynthonSpaceTask load_a = new LoadSynthonSpaceTask(model,view,"C:\\dev_github\\openchemlib-hyperspace\\test_space.data");
        //LoadSynthonSpaceTask load_a = new LoadSynthonSpaceTask(model,view,"C:\\dev_github\\openchemlib-hyperspace\\divchem_1s_FragFp.data");
        //LoadSynthonSpaceTask load_a = new LoadSynthonSpaceTask(model,"C:\\dev_github\\openchemlib-hyperspace\\realspace_FragFp.data");

        //model.addTask(load_a);
        //load_a.execute();
    }

    /**
     * TODO: extend at some point to some reasonable config..
     * @param file
     * @param model
     */
    private static void parseConfig(String file, LeetHyperspaceModel model) {
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = null;
            while( (line=in.readLine()) != null ) {
                try{
                    File fi = new File(line);
                    LoadSynthonSpaceTask load_a = new LoadSynthonSpaceTask(model,fi.getAbsolutePath());
                    model.addTask(load_a);
                    load_a.execute();
                }
                catch(Exception ex2) {
                    ex2.printStackTrace();
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
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
                model.getSynthonSpaces().stream().forEach( xi -> jmSearchItems.add(new JMenuItem( new RunSubstructureSearchAction(model,view,xi.getSpace(),xi.getName(),() -> model.getQuery()))));
                jmSearchItems.stream().forEach( xi -> jmSearch.add(xi));
            }
        });

        fi.setJMenuBar(jb);
    }

}
