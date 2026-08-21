package com.idorsia.research.chem.hyperspace3d.batch;

import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace3d.index.ProductTuple;
import java.util.List;

@FunctionalInterface
public interface SynthonTupleResolver {
    List<StereoMolecule> resolve(ProductTuple tuple);
}
