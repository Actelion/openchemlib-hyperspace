package com.idorsia.research.chem.hyperspace3d.mining;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;

/** Thread-confined exact PheSA descriptor and pair scorer. */
public final class ExactPheSALabeler {
    private static final double COMPONENT_TOLERANCE = 2e-6;
    private final int maxConformers;
    private final double pharmacophoreWeight;
    private final DescriptorHandlerShape handler;

    public ExactPheSALabeler(int maxConformers, double pharmacophoreWeight) {
        if (maxConformers < 1) throw new IllegalArgumentException("maxConformers must be positive");
        if (!Double.isFinite(pharmacophoreWeight) || pharmacophoreWeight < 0.0
                || pharmacophoreWeight > 1.0) {
            throw new IllegalArgumentException("pharmacophoreWeight must be within [0,1]");
        }
        this.maxConformers = maxConformers;
        this.pharmacophoreWeight = pharmacophoreWeight;
        this.handler = new DescriptorHandlerShape(false, maxConformers, pharmacophoreWeight);
        this.handler.getPhesaSetting().setPpWeight(pharmacophoreWeight);
    }

    public PheSAMolecule describe(StereoMolecule molecule) {
        if (molecule == null) return null;
        PheSAMolecule result = handler.createDescriptor(new StereoMolecule(molecule));
        return result == null || handler.calculationFailed(result) ? null : result;
    }

    public String encode(PheSAMolecule descriptor) {
        if (descriptor == null || handler.calculationFailed(descriptor)) return null;
        return handler.encode(descriptor);
    }

    public PheSAMolecule decode(String encoded) {
        if (encoded == null || encoded.isBlank()) return null;
        PheSAMolecule value = handler.decode(encoded);
        return value == null || handler.calculationFailed(value) ? null : value;
    }

    public ExactPheSALabel score(PheSAMolecule query, PheSAMolecule candidate) {
        long started = System.nanoTime();
        if (query == null || candidate == null) {
            return ExactPheSALabel.failure(PheSAQueryPairRecord.EXACT_DESCRIPTOR_FAILED,
                    System.nanoTime() - started);
        }
        try {
            float total = handler.getSimilarity(query, candidate);
            double[] components = handler.getPreviousPheSAResult();
            if (components == null || components.length < 3) {
                return ExactPheSALabel.failure(PheSAQueryPairRecord.EXACT_ALIGNMENT_FAILED,
                        System.nanoTime() - started);
            }
            float pharmacophore = (float) components[1];
            float shape = (float) components[2];
            if (!Float.isFinite(total) || !Float.isFinite(shape)
                    || !Float.isFinite(pharmacophore)) {
                return ExactPheSALabel.failure(PheSAQueryPairRecord.EXACT_NONFINITE,
                        System.nanoTime() - started);
            }
            double reconstructed = (1.0 - pharmacophoreWeight) * shape
                    + pharmacophoreWeight * pharmacophore;
            if (Math.abs(total - reconstructed) > COMPONENT_TOLERANCE) {
                return ExactPheSALabel.failure(PheSAQueryPairRecord.EXACT_COMPONENT_MISMATCH,
                        System.nanoTime() - started);
            }
            return new ExactPheSALabel(total, shape, pharmacophore,
                    PheSAQueryPairRecord.EXACT_OK, System.nanoTime() - started);
        } catch (RuntimeException error) {
            return ExactPheSALabel.failure(PheSAQueryPairRecord.EXACT_ALIGNMENT_FAILED,
                    System.nanoTime() - started);
        }
    }

    public int maxConformers() { return maxConformers; }
    public double pharmacophoreWeight() { return pharmacophoreWeight; }
}
