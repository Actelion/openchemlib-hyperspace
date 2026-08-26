package com.idorsia.research.chem.hyperspace.tools.chembl;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;
import java.util.zip.GZIPOutputStream;

/** Writes latent assembly decompositions as one self-contained JSON object per line. */
public final class LatentAssemblyJsonlWriter implements AutoCloseable {
    private final BufferedWriter writer;

    public LatentAssemblyJsonlWriter(Path output) throws IOException {
        Path parent = output.toAbsolutePath().getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        OutputStream stream = Files.newOutputStream(output);
        if (output.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".gz")) {
            stream = new GZIPOutputStream(stream);
        }
        writer = new BufferedWriter(new OutputStreamWriter(stream, StandardCharsets.UTF_8));
    }

    public synchronized void write(LatentAssemblyDecomposition decomposition) throws IOException {
        writer.write(toJson(decomposition).toString());
        writer.newLine();
    }

    static JSONObject toJson(LatentAssemblyDecomposition decomposition) {
        JSONObject json = new JSONObject();
        json.put("schema_version", decomposition.schemaVersion());
        json.put("source_molecule_id", decomposition.sourceMoleculeId());
        json.put("product_idcode", decomposition.productIdcode());
        json.put("product_atom_count", decomposition.productAtomCount());
        json.put("product_non_hydrogen_atom_count", decomposition.productNonHydrogenAtomCount());
        json.put("scaffold_product_atoms", array(decomposition.scaffoldProductAtoms()));

        JSONArray arms = new JSONArray();
        for (LatentAssemblyDecomposition.Arm arm : decomposition.arms()) {
            arms.put(new JSONObject()
                    .put("index", arm.index())
                    .put("label", arm.label())
                    .put("product_atoms", array(arm.productAtoms())));
        }
        json.put("arms", arms);
        json.put("scaffold_view", viewJson(decomposition.scaffoldView()));

        JSONArray views = new JSONArray();
        decomposition.assemblyViews().forEach(view -> views.put(viewJson(view)));
        json.put("assembly_views", views);

        JSONArray cutBonds = new JSONArray();
        for (LatentAssemblyDecomposition.CutBond cut : decomposition.cutBonds()) {
            cutBonds.put(new JSONObject()
                    .put("product_bond_index", cut.productBondIndex())
                    .put("scaffold_product_atom", cut.scaffoldProductAtom())
                    .put("arm_product_atom", cut.armProductAtom())
                    .put("arm_index", cut.armIndex())
                    .put("bond_type", cut.bondType())
                    .put("bond_order", cut.bondOrder()));
        }
        json.put("cut_bonds", cutBonds);

        LatentAssemblyDecomposition.SelectionMetrics metrics = decomposition.selectionMetrics();
        json.put("selection_metrics", new JSONObject()
                .put("scaffold_non_hydrogen_atoms", metrics.scaffoldNonHydrogenAtoms())
                .put("scaffold_ring_atoms", metrics.scaffoldRingAtoms())
                .put("arm_non_hydrogen_atoms", array(metrics.armNonHydrogenAtoms()))
                .put("minimum_attachment_distance", metrics.minimumAttachmentDistance())
                .put("scaffold_fraction", metrics.scaffoldFraction()));
        return json;
    }

    private static JSONObject viewJson(LatentAssemblyDecomposition.View view) {
        return new JSONObject()
                .put("label", view.label())
                .put("idcode", view.idcode())
                .put("product_atom_by_view_atom", array(view.productAtomByViewAtom()))
                .put("open_scaffold_attachment_view_atoms", array(view.openScaffoldAttachmentViewAtoms()));
    }

    private static JSONArray array(List<Integer> values) {
        JSONArray array = new JSONArray();
        values.forEach(array::put);
        return array;
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
