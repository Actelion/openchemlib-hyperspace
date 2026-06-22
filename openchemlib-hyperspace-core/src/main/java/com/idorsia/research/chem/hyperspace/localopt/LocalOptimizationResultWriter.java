package com.idorsia.research.chem.hyperspace.localopt;

import com.idorsia.research.chem.hyperspace.SynthonSpace;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public class LocalOptimizationResultWriter implements AutoCloseable {

    private final BufferedWriter writer;
    private final Function<String, String> reactionSourceProvider;
    private final boolean includeSourceSpaceColumn;
    private boolean headerWritten = false;

    public LocalOptimizationResultWriter(Path outputPath) throws IOException {
        this(outputPath, null);
    }

    public LocalOptimizationResultWriter(Path outputPath, Function<String, String> reactionSourceProvider) throws IOException {
        Path parent = outputPath.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        this.writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8);
        this.reactionSourceProvider = reactionSourceProvider;
        this.includeSourceSpaceColumn = reactionSourceProvider != null;
    }

    public void write(LocalOptimizationResult result) throws IOException {
        if (!headerWritten) {
            writer.write("rxnId");
            if (includeSourceSpaceColumn) {
                writer.write("	sourceSpace");
            }
            writer.write("	fragIds	Structure [idcode]	atoms	rotatableBonds	phesaSimilarity	attemptIndex	seedFragIds");
            writer.newLine();
            headerWritten = true;
        }
        for (LocalOptimizationResult.BeamEntry entry : result.getBeamEntries()) {
            writer.write(result.getReactionId());
            if (includeSourceSpaceColumn) {
                writer.write('	');
                writer.write(sanitizeCell(reactionSourceProvider.apply(result.getReactionId())));
            }
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

    private String sanitizeCell(String value) {
        if (value == null || value.isBlank()) {
            return "";
        }
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }
}
