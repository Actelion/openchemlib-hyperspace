package com.idorsia.research.chem.hyperspace3d.archive;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public final class CandidateArchiveWriter implements Closeable {
    private final ObjectMapper mapper = new ObjectMapper();
    private final BufferedWriter output;

    public CandidateArchiveWriter(Path path) throws IOException {
        output = Files.newBufferedWriter(path);
    }

    public void write(CandidateArchiveRecord record) throws IOException {
        output.write(mapper.writeValueAsString(record));
        output.newLine();
    }

    @Override public void close() throws IOException { output.close(); }
}
