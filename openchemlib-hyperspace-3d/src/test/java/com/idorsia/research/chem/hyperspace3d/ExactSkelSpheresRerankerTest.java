package com.idorsia.research.chem.hyperspace3d;

import static org.junit.jupiter.api.Assertions.assertEquals;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import com.idorsia.research.chem.hyperspace3d.index.ProductVectorReference;
import com.idorsia.research.chem.hyperspace3d.screening.CompactSkelSpheresHit;
import com.idorsia.research.chem.hyperspace3d.screening.ExactSkelSpheresReranker;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class ExactSkelSpheresRerankerTest {
    @Test void authoritativeSynthonsAreAssembledAndScored() throws Exception {
        RawSynthon a = synthon(0, "a", "[U]CCC");
        RawSynthon b = synthon(1, "b", "[U]CCO");
        RawSynthonSpace space = RawSynthonSpace.builder("space")
                .addRawFragments("r", 0, List.of(a))
                .addRawFragments("r", 1, List.of(b)).build();
        List<StereoMolecule> parts = new ArrayList<>();
        IDCodeParser parser = new IDCodeParser();
        for (RawSynthon raw : List.of(a, b)) {
            StereoMolecule part = new StereoMolecule(); parser.parse(part, raw.getIdcode()); parts.add(part);
        }
        StereoMolecule query = SynthonAssembler.assembleSynthons_faster(parts);
        var tuple = new ProductTuple("r", List.of("a", "b"), List.of(0, 0));
        var hit = new CompactSkelSpheresHit(new ProductVectorReference(0, 0, 0, 24),
                tuple, 0.5, 0.6, null);
        var result = new ExactSkelSpheresReranker().rerank(query, List.of(hit), space);
        assertEquals(1.0, result.get(0).exactSimilarity(), 1e-12);
    }
    private static RawSynthon synthon(int position, String id, String smiles) throws Exception {
        StereoMolecule molecule = new StereoMolecule(); new SmilesParser().parse(molecule, smiles);
        return RawSynthon.fromMolecule("r", position, id, molecule);
    }
}
