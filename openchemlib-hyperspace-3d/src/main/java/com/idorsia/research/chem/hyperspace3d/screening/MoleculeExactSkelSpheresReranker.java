package com.idorsia.research.chem.hyperspace3d.screening;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Exact OCL BinarySkelSpheres reranking for resolved flat-molecule hits. */
public final class MoleculeExactSkelSpheresReranker {
    public List<MoleculeSkelSpheresHit> rerank(StereoMolecule query,
            List<MoleculeSkelSpheresHit> hits) throws Exception {
        var handler = new DescriptorHandlerBinarySkelSpheres();
        query.ensureHelperArrays(Molecule.cHelperCIP);
        int[] queryDescriptor = handler.createDescriptor(query);
        var parser = new SmilesParser();
        List<MoleculeSkelSpheresHit> result = new ArrayList<>(hits.size());
        for (var hit : hits) {
            StereoMolecule candidate = new StereoMolecule();
            parser.parse(candidate, hit.molecule().smiles());
            candidate.ensureHelperArrays(Molecule.cHelperCIP);
            double exact = handler.getSimilarity(queryDescriptor,
                    handler.createDescriptor(candidate));
            result.add(hit.withExact(exact));
        }
        result.sort(Comparator.comparing(MoleculeSkelSpheresHit::exactSimilarity,
                Comparator.nullsLast(Comparator.reverseOrder()))
                .thenComparing(Comparator.comparingDouble(
                        MoleculeSkelSpheresHit::predictedSimilarity).reversed())
                .thenComparingInt(hit -> hit.reference().shardIndex())
                .thenComparingLong(hit -> hit.reference().localRow()));
        return result;
    }
}
