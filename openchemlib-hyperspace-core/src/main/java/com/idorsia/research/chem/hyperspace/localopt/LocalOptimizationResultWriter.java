package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class LocalOptimizationResultWriter implements AutoCloseable {

    private final BufferedWriter writer;
    private boolean headerWritten = false;

    public LocalOptimizationResultWriter(Path outputPath) throws IOException {
        Path parent = outputPath.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        this.writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8);
    }

    public void write(LocalOptimizationResult result) throws IOException {
        if (!headerWritten) {
            writer.write("rxnId	fragIds	assembledIdcode	atoms	rotatableBonds	phesaSimilarity	attemptIndex	seedFragIds");
            writer.newLine();
            headerWritten = true;
        }
        for (LocalOptimizationResult.BeamEntry entry : result.getBeamEntries()) {
            writer.write(result.getReactionId());
            writer.write('	');
            writer.write(joinFragments(entry.getFragments()));
            writer.write('	');
            writer.write(entry.getIdcode());
            writer.write('	');
            writer.write(Integer.toString(entry.getAtomCount()));
            writer.write('	');
            writer.write(Integer.toString(entry.getRotatableBonds()));
            writer.write('	');
            writer.write(Double.toString(entry.getPhesaSimilarity()));
            writer.write('	');
            writer.write(Integer.toString(entry.getOriginatingRound()));
            writer.write('	');
            writer.write(String.join(";", result.getSeedFragments()));
            writer.newLine();
        }
        writer.flush();
    }

    private String joinFragments(List<SynthonSpace.FragId> fragments) {
        return fragments.stream().map(f -> f.fragment_id).collect(Collectors.joining(";"));
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
