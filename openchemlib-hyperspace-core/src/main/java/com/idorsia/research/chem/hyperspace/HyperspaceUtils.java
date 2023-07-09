package com.idorsia.research.chem.hyperspace;

import com.actelion.research.chem.*;
import com.actelion.research.gui.DefaultCompoundCollectionModel;
import com.actelion.research.gui.JStructureView;
import com.actelion.research.gui.generic.GenericRectangle;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;

public class HyperspaceUtils {



    public static int setHighlightedSubstructure(StereoMolecule M, StereoMolecule query) {
        return setHighlightedSubstructure(M,query,true);
    }

    /**
     * Sets any set of bonds that match the query substructure
     * to highlighted.
     *
     * @param M
     * @param query
     */
    public static int setHighlightedSubstructure(StereoMolecule M, StereoMolecule query, boolean highlighAllMatches) {
        StereoMolecule fi = query;

        //fi.setFragment(true);
        SSSearcher sss = new SSSearcher();
        sss.setFragment(fi);
        sss.setMolecule(M);
        sss.findFragmentInMolecule();
        for (int zi = 0; zi < M.getBonds(); zi++) {
            M.setBondBackgroundHiliting(zi, false);
        }
        ArrayList<int[]> matches = sss.getMatchList();
        for (int[] match_i : matches) {
            // 1. determine bonds:
            BitSet atoms = new BitSet(M.getAtoms());
            Arrays.stream(match_i).forEach(ai -> atoms.set(ai));
            BitSet bonds = new BitSet(M.getBonds());
            for (int zi = 0; zi < M.getBonds(); zi++) {
                int baa = M.getBondAtom(0, zi);
                int bab = M.getBondAtom(1, zi);
                if (atoms.get(baa) && atoms.get(bab)) {
                    bonds.set(zi, true);
                }
            }
            // 2. highlight bonds:
            for (int zi = 0; zi < M.getBonds(); zi++) {
                if (bonds.get(zi)) {
                    M.setBondBackgroundHiliting(zi, true);
                }
            }
            if(!highlighAllMatches) {
                break;
            }
        }
        return matches.size();
    }

    public static int setSelectedSubstructure(StereoMolecule M, StereoMolecule query) {
        StereoMolecule fi = query;

        //fi.setFragment(true);
        SSSearcher sss = new SSSearcher();
        sss.setFragment(fi);
        sss.setMolecule(M);
        sss.findFragmentInMolecule();
        for (int zi = 0; zi < M.getBonds(); zi++) {
            M.setBondBackgroundHiliting(zi, false);
        }
        ArrayList<int[]> matches = sss.getMatchList();
        for (int[] match_i : matches) {
            // 1. select atoms:
            //BitSet atoms = new BitSet(M.getAtoms());
            Arrays.stream(match_i).forEach(ai -> M.setAtomSelection(ai,true));
        }
        return matches.size();
    }


    public static String readFileIntoString(String file) {
        String content = null;
        try {
            content = new String(Files.readAllBytes(Paths.get(file)));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return content;
    }

    public static List<StereoMolecule> parseInputFile(String in) throws IOException {
        List<StereoMolecule> parsed = new ArrayList<>();

        File fi = new File(in);
        BufferedReader bri = new BufferedReader( new FileReader( fi ) );

        if( in.endsWith(".smi") || in.endsWith(".smiles") || in.endsWith(".SMI") ) {
            SmilesParser spi = new SmilesParser();
            System.out.println("Parse smiles input file "+fi.getName());
            String line = null;
            while( (line=bri.readLine()) != null ) {
                try{
                    StereoMolecule mi = new StereoMolecule();
                    spi.parse(mi,line);
                    parsed.add(mi);
                }
                catch(Exception ex){
                }
            }
        }
        else {
            System.out.println("Parse idcodes input file "+fi.getName());
            IDCodeParser idp = new IDCodeParser();
            String line = null;
            while( (line=bri.readLine()) != null ) {
                try{
                    StereoMolecule mi = new StereoMolecule();
                    idp.parse(mi,line);
                    parsed.add(mi);
                }
                catch(Exception ex){
                }
            }
        }

        return parsed;
    }


    /**
     * Parses for each line two parts: first must be smiles or idcode,
     * then second (separated by any whitespace-like characters) is
     * an ID.
     *
     * If there is no ID supplied, then a random number will be generated.
     * The ID will be sanitized, such that it does not contain any characters
     * that cannot be used as a folder name in file systems.
     *
     * @param in
     * @return
     * @throws IOException
     */
    public static Map<String,StereoMolecule> parseInputFile_2(String in) throws IOException {
        Map<String,StereoMolecule> parsed = new HashMap<>();

        Random r = new Random();

        File fi = new File(in);
        BufferedReader bri = new BufferedReader( new FileReader( fi ) );

        if( in.endsWith(".smi") || in.endsWith(".smiles") || in.endsWith(".SMI") ) {
            SmilesParser spi = new SmilesParser();
            System.out.println("Parse smiles input file "+fi.getName());
            String line = null;
            while( (line=bri.readLine()) != null ) {
                try{
                    String parsed_line[] = parseInputFile_2_processLine(r,line,parsed);
                    String mol = parsed_line[0];
                    String sisi = parsed_line[1];
                    StereoMolecule mi = new StereoMolecule();
                    spi.parse(mi,mol);
                    parsed.put(sisi,mi);
                }
                catch(Exception ex){
                }
            }
        }
        else {
            System.out.println("Parse idcodes input file "+fi.getName());
            IDCodeParser idp = new IDCodeParser();
            String line = null;
            while( (line=bri.readLine()) != null ) {
                try{
                    String parsed_line[] = parseInputFile_2_processLine(r,line,parsed);
                    String mol = parsed_line[0];
                    String sisi = parsed_line[1];

                    StereoMolecule mi = new StereoMolecule();
                    idp.parse(mi,mol);
                    parsed.put(sisi,mi);
                }
                catch(Exception ex){
                }
            }
        }

        return parsed;
    }

    /**
     * returns first part in result[0], and a sanitized and unique id in result[2]
     *
     * @param line
     * @return
     */
    private static String[] parseInputFile_2_processLine(Random r, String line, Map<String,StereoMolecule> parsed) {
        String splits[] = line.split("\\s+");
        String mol = splits[0];
        String idi = "";
        //now build id:
        if(splits.length==1) {}
        else {
            for(int zi=1;zi<splits.length;zi++){ idi += splits[zi]; }
        }
        // now sanitize:
        char[] cidi = idi.toCharArray();
        for(int zi=0;zi<cidi.length;zi++) {
            if(! ( Character.isLetterOrDigit(cidi[zi]) || cidi[zi]=='-' || cidi[zi]=='.' ) ){ cidi[zi]='_';}
        }
        String sisi = new String(cidi);
        // now check that we have unique string..
        while(parsed.containsKey(sisi)) {
            sisi += "_"+String.format( "%06d",r.nextInt(1000000) );
        }

        return new String[]{ mol,sisi };
    }


    public static String idcodeToSmiles(String idcode) {
        StereoMolecule m = new StereoMolecule();
        IDCodeParser p = new IDCodeParser();
        p.parse(m,idcode);

        String smiles = "exception_in_smiles_creator";
        try{
            IsomericSmilesCreator smiles_creator = new IsomericSmilesCreator(m);
            smiles=smiles_creator.getSmiles();
        }
        catch(Exception ex) {
            System.out.println("Exception in idcodeToSmiles..");
        }
        return smiles;
    }

    public static String stereoMoleculeToSmiles(StereoMolecule sm) {
        IsomericSmilesCreator smiles_creator = new IsomericSmilesCreator(sm);
        return smiles_creator.getSmiles();
    }


    public static StereoMolecule parseIDCode(String idc) {
        IDCodeParser icp = new IDCodeParser();
        //StereoMolecule mi = new StereoMolecule();
        //icp.parse(mi,idc);
        StereoMolecule mi = icp.getCompactMolecule(idc);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        return mi;
    }

    public static StereoMolecule removeConnectorsFromStructure(StereoMolecule mol, boolean setRequireSubstitution) {
        StereoMolecule mol2 = new StereoMolecule(mol);
        mol2.ensureHelperArrays(Molecule.cHelperCIP);

        if(setRequireSubstitution) {
            // remove all potential connectors and exit vector atoms, and set query features accordingly
            mol2.setFragment(true);
            mol2.ensureHelperArrays(Molecule.cHelperCIP);

            for (int za = 0; za < mol2.getAtoms(); za++) {
                if (mol2.getAtomicNo(za) >= 88) {
                    // and set require substitution query feature for neighbor atom
                    if (mol2.getConnAtoms(za) != 1) {
                        // this is an error..
                        System.out.println("[ERROR] encountered unexpected scaffold, neighbor count of " + za + " = " + mol2.getNonHydrogenNeighbourCount(za) + " but must be 1");
                        //List<RepresentativeCombinatorialHits> hits_empty = new ArrayList<>();
                        //return hits_empty;
                        //return new RepresentativeCombinatorialHits("","",null,null);
                    } else {
                        // set requires substitution query feature
                        int atom_a = mol2.getConnAtom(za, 0);
                        mol2.setAtomQueryFeature(atom_a, StereoMolecule.cAtomQFMoreNeighbours, true);
                    }
                }
            }
        }

        for(int zi=0;zi<mol2.getAtoms();zi++) {
            if(mol2.getAtomicNo(zi)>= 89) {
                mol2.markAtomForDeletion(zi);
            }
        }

        mol2.deleteMarkedAtomsAndBonds();
        mol2.ensureHelperArrays(Molecule.cHelperCIP);

        return mol2;
    }


    /**
     * prints:
     * Structure[idcode],SynthonRxn,FragId,Frag
     *
     * @param separator
     * @param fid
     * @return
     */
    public static String printSynthonToCSV(String separator, SynthonSpace.FragId fid) {
        List<String> parts = new ArrayList<>();
        parts.add(fid.idcode);
        parts.add(fid.rxn_id);
        parts.add(fid.fragment_id);
        parts.add(""+fid.frag);
        return String.join(separator,parts);
    }

    public static void createPNGFromStructure(StereoMolecule mi, String path, int width, int height) throws IOException, InvocationTargetException, InterruptedException {
        BufferedImage bi = new BufferedImage(width,height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g_1 = bi.createGraphics();

        g_1.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        //g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
        g_1.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
        g_1.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);

        Depictor2D pic = new Depictor2D(mi);
        //Rectangle.Double rect = new Rectangle2D.Double(0,0,width,height);
        GenericRectangle rect = new GenericRectangle(0,0,width,height);
        //pic.setTransformation(new DepictorTransformation(4,0,0));
        pic.validateView(g_1,rect,Depictor2D.cModeInflateToMaxAVBL | Depictor2D.cDModeSuppressChiralText);

        //pic.setDisplayMode(Depictor2D.cDModeNoTabus);
        pic.paint(g_1);



        //Rectangle.Double bounds = pic.getBoundingRect();
        GenericRectangle bounds = pic.getBoundingRect();
        g_1.setColor(Color.blue);
        g_1.drawRect((int)bounds.x,(int)bounds.y,(int)bounds.width,(int)bounds.height);

        Depictor2D pic2 = new Depictor2D(mi);
        //pic.setTransformation(new DepictorTransformation(4,0,0));

        //Rectangle.Double bounds_extended = new Rectangle2D.Double(0,0,bounds.width+1,bounds.height+1);
        GenericRectangle bounds_extended = new GenericRectangle(0,0,bounds.width+1,bounds.height+1);

        BufferedImage bi_2 = new BufferedImage((int)Math.ceil(bounds.width)+1,(int)Math.ceil(bounds.height)+1, BufferedImage.TYPE_INT_ARGB);
        Graphics2D    g_2  = bi_2.createGraphics();
        g_2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g_2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
        g_2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);

        pic2.validateView(g_2,bounds_extended,Depictor2D.cModeInflateToMaxAVBL | Depictor2D.cDModeSuppressChiralText);
        //pic.setDisplayMode(Depictor2D.cDModeNoTabus);
        pic2.paint(g_2);

        // BufferedImage bi_2 = new BufferedImage((int)Math.ceil(bounds.width),(int)Math.ceil(bounds.height), BufferedImage.TYPE_INT_ARGB);
        // Graphics2D    g_2  = bi_2.createGraphics();
        // pic.validateView(g_2,bounds,Depictor2D.cModeInflateToMaxAVBL);
        // pic.paint(g_2);

        g_1.dispose();
        g_2.dispose();
        try {
            if(true) {
                ImageIO.write(bi, "png", new File(path.replace(".png","_a.png")));
                ImageIO.write(bi_2, "png", new File(path));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static StereoMolecule parseSmiles(String smiles) throws Exception {
        SmilesParser spi   = new SmilesParser();
        StereoMolecule mi  = new StereoMolecule();
        spi.parse(mi,smiles);
        mi.ensureHelperArrays(Molecule.cHelperCIP);
        return mi;
    }


    public static class DebugOutput {

        static JFrame frame = null;
        static JTabbedPane tp = null;

        public static void plotMolecules(String title, List<String> idcodes, int sx, int sy) {
            plotMolecules(title,idcodes.stream().toArray(String[]::new),sx,sy);
        }

        public static void plotMolecules(String title, String idcodes[], int sx, int sy) {

            if (frame == null) {
                JFrame f = new JFrame();
                frame = f;
                tp = new JTabbedPane();
                f.getContentPane().setLayout(new BorderLayout());
                f.getContentPane().add(tp);
                f.setVisible(true);

                int w = Math.min(800, sx * 300);
                int h = Math.min(800, sx * 300);
                f.setSize(w, h);
            }

            JScrollPane sp = new JScrollPane();
            JPanel pi = new JPanel();
            pi.setLayout(new GridLayout(sx, sy));
            sp.setViewportView(pi);

            for (int zi = 0; zi < idcodes.length; zi++) {
                String si = idcodes[zi];
                JStructureView sv = new JStructureView();
                sv.setIDCode(si);
                pi.add(sv);
                sv.setBorder(new LineBorder(Color.black));
            }

            tp.add(title, pi);
        }


        public static void main(String args[]) {
            String sa = "didHPD@zxHR[Y^FZZX@`";
            String sb = "dmuHPHF`neNdefuQfjjj`B";
            String idcodes[] = new String[]{sa, sb};
            plotMolecules("Out_A", idcodes, 2, 2);
            plotMolecules("Out_b", idcodes, 2, 2);
            plotMolecules("Out_c", idcodes, 2, 2);
        }
    }


}
