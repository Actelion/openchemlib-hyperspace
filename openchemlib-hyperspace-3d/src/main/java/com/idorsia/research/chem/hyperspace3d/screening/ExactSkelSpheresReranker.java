package com.idorsia.research.chem.hyperspace3d.screening;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerBinarySkelSpheres;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Optional CPU reranking boundary; it is deliberately not part of the default learned scan. */
public final class ExactSkelSpheresReranker {
    public List<CompactSkelSpheresHit> rerank(StereoMolecule query,
            List<CompactSkelSpheresHit> hits, RawSynthonSpace fullSpace) {
        var handler = new DescriptorHandlerBinarySkelSpheres();
        int[] queryDescriptor = handler.createDescriptor(query);
        List<CompactSkelSpheresHit> result = new ArrayList<>(hits.size());
        IDCodeParser parser = new IDCodeParser();
        for (CompactSkelSpheresHit hit : hits) {
            var reaction = fullSpace.getReactions().get(hit.tuple().reactionId());
            if (reaction == null) throw new IllegalArgumentException("unknown reaction " + hit.tuple().reactionId());
            List<StereoMolecule> parts = new ArrayList<>();
            for (String id : hit.tuple().synthonIds()) {
                RawSynthon raw = reaction.findRawFragment(id);
                if (raw == null) throw new IllegalArgumentException("unknown synthon " + id);
                StereoMolecule part = new StereoMolecule(); parser.parse(part, raw.getIdcode());
                part.ensureHelperArrays(Molecule.cHelperCIP); parts.add(part);
            }
            StereoMolecule product = SynthonAssembler.assembleSynthons_faster(parts);
            double exact = handler.getSimilarity(queryDescriptor, handler.createDescriptor(product));
            result.add(hit.withExact(exact));
        }
        result.sort(Comparator.comparing(CompactSkelSpheresHit::exactSimilarity,
                Comparator.nullsLast(Comparator.reverseOrder()))
                .thenComparing(Comparator.comparingDouble(
                        CompactSkelSpheresHit::predictedSimilarity).reversed()));
        return result;
    }
}
