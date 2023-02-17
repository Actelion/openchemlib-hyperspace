package com.idorsia.research.chem.hyperspace.util;

import com.actelion.research.chem.DiversitySelector;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class RelevantSubgraphSampler {

    private StereoMolecule M;

    public RelevantSubgraphSampler(StereoMolecule m) {
        this.M = new StereoMolecule(m);
        this.M.ensureHelperArrays(Molecule.cHelperCIP);
    }


    public List<RelevantSubgraph> sampleSubgraphs_A(int min_size, int max_size, int min_bits_ffp, int max_bits_ffp, int numSubgraphs) {

        DescriptorHandlerLongFFP512 dh = (DescriptorHandlerLongFFP512) DescriptorHandlerLongFFP512.getDefaultInstance();

        //List<Pair<BitSet, StereoMolecule>> result = new ArrayList<>();
        List<RelevantSubgraph> result = new ArrayList<>();

        FragmentGenerator fg = new FragmentGenerator(M, min_size, max_size);
        List<BitSet> frags = fg.getFragments();
        Map<BitSet, StereoMolecule> fragMolecules = new HashMap<>();
        Map<BitSet, long[]> fragFFPs = new HashMap<>();

        List<BitSet> frags_filtered = new ArrayList<>();

        // compute ffps and filter by number of ffp bits
        for (int zi = 0; zi < frags.size(); zi++) {
            StereoMolecule mi = new StereoMolecule();
            M.copyMoleculeByAtoms(mi, ConnectedPairSampler.ChemUtils.toBooleanArray(frags.get(zi),M.getAtoms()), true, null);
            mi.ensureHelperArrays(Molecule.cHelperCIP);
            fragMolecules.put(frags.get(zi), mi);

            long[] ffpi = dh.createDescriptor(mi);
            BitSet bsffpi = BitSet.valueOf(ffpi);
            int card_ffpi = bsffpi.cardinality();
            if(card_ffpi>=min_bits_ffp && card_ffpi <= max_bits_ffp) {
                frags_filtered.add(frags.get(zi));
                fragFFPs.put(frags.get(zi),ffpi);
            }
        }

        // prepare descriptor array for :
        long[][] descriptors = new long[frags_filtered.size()][];

        for(int zi=0;zi<frags_filtered.size();zi++) {
            descriptors[zi] = fragFFPs.get(frags_filtered.get(zi));
        }

        DiversitySelector selector = new DiversitySelector();
        selector.initializeExistingSet(512);

        int[] selected = selector.select(descriptors, numSubgraphs);

        for (int si : selected) {
            result.add( new RelevantSubgraph(M, fragMolecules.get(frags_filtered.get(si)), frags_filtered.get(si), BitSet.valueOf( fragFFPs.get(frags_filtered.get(si))) ));
            if (true) {
                System.out.println(fragMolecules.get(frags.get(si)).getIDCode());
            }
        }

        return result;
    }


    public static class RelevantSubgraph {
        private StereoMolecule m;
        private StereoMolecule frag;
        private BitSet fragBS;

        private BitSet fragFP;

        public RelevantSubgraph(StereoMolecule m, StereoMolecule frag, BitSet fragBS, BitSet fragfp) {
            this.m = m;
            this.frag = frag;
            this.fragBS = fragBS;
            this.fragFP = fragfp;
        }

        public int getBitsFragFP() { return this.getFragFP().cardinality(); }


        public StereoMolecule getM() {
            return m;
        }

        public StereoMolecule getFrag() {
            return frag;
        }

        public BitSet getFragBS() {
            return fragBS;
        }

        public BitSet getFragFP() {
            return fragFP;
        }
    }

    public static RelevantSubgraph createRelevantSubgraph(StereoMolecule m, BitSet atoms) {
        StereoMolecule fa = new StereoMolecule();
        int map_fa[] = new int[m.getAtoms()];
        m.copyMoleculeByAtoms(fa, ConnectedPairSampler.ChemUtils.toBooleanArray(atoms, m.getAtoms()), true, map_fa);

        fa.ensureHelperArrays(Molecule.cHelperCIP);
        BitSet fragfp = BitSet.valueOf(DescriptorHandlerLongFFP512.getDefaultInstance().createDescriptor(fa));

        return new RelevantSubgraph(m,fa,atoms,fragfp);
    }

    public static void main(String args[]) {
        String mol_a = "O=C(Nc1ncnc2nc[nH]c12)c1cccc(-c2cccc(-c3cn(CP(=O)(O)O)cn3)c2)c1"; // from chembl

        StereoMolecule ma = null;
        try {
            ma = HyperspaceUtils.parseSmiles(mol_a);
        } catch (Exception e) {
            e.printStackTrace();
        }

        RelevantSubgraphSampler rss = new RelevantSubgraphSampler(ma);

        //List<ConnectedPairOfSubgraphs> pairs = cps.computePairs(10, 10, 3, 7);
        List<RelevantSubgraph> subgraphs = rss.sampleSubgraphs_A(16, 20, 60, 120, 16);

        System.out.println("mkay..");

        //Collections.shuffle(pairs);
        //for (int zi = 0; zi < Math.min(1000, subgraphs.size()); zi++) {
        //    System.out.println(subgraphs.get(zi).getFrag().getIDCode() + " : ffpbits="+subgraphs.get(zi).getBitsFragFP());
        //}

        try(BufferedWriter out = new BufferedWriter(new FileWriter("out_relevant_subgraph_sampler_a.csv"))) {
            out.write("Structure[idcode],FFPBits\n");
            for (int zi = 0; zi < Math.min(1000, subgraphs.size()); zi++) {
                out.write(subgraphs.get(zi).frag.getIDCode()+","+subgraphs.get(zi).getBitsFragFP()+"\n");
            }
            out.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

}
