package com.idorsia.research.chem.hyperspace.seed;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Thread-safe TSV writer with optional file rotation after a configurable number of hits.
 */
class SeedHitWriter implements AutoCloseable {

    private final Path basePath;
    private final int maxHitsPerFile;
    private BufferedWriter currentWriter;
    private final AtomicInteger fileIndex = new AtomicInteger();
    private int hitsInCurrentFile = 0;

    SeedHitWriter(Path basePath, int maxHitsPerFile) {
        this.basePath = basePath;
        this.maxHitsPerFile = maxHitsPerFile;
    }

    synchronized void writeHit(HitRecord hit) {
        try {
            rotateWriterIfNeeded();
            currentWriter.write(format(hit));
            currentWriter.newLine();
            hitsInCurrentFile++;
            currentWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException("Unable to write seed hit", e);
        }
    }

    private String format(HitRecord hit) {
        String fragmentIds = String.join(";", hit.fragmentIds);
        return String.join("\t",
                hit.reactionId,
                fragmentIds,
                hit.assembledIdCode,
                Integer.toString(hit.atomCount),
                Integer.toString(hit.rotatableBonds),
                Double.toString(hit.phesaSimilarity),
                Long.toString(hit.attemptIndex));
    }

    private void rotateWriterIfNeeded() throws IOException {
        if (currentWriter == null) {
            currentWriter = newWriter(basePath);
            writeHeader();
            return;
        }
        if (maxHitsPerFile > 0 && hitsInCurrentFile >= maxHitsPerFile) {
            currentWriter.close();
            Path rotated = rotatePath(basePath, fileIndex.incrementAndGet());
            currentWriter = newWriter(rotated);
            hitsInCurrentFile = 0;
            writeHeader();
        }
    }

    private void writeHeader() throws IOException {
        currentWriter.write("rxnId\tfragIds\tassembledIdcode\tatoms\trotatableBonds\tphesaSimilarity\tattemptIndex");
        currentWriter.newLine();
    }

    private static BufferedWriter newWriter(Path path) throws IOException {
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        return Files.newBufferedWriter(path, StandardCharsets.UTF_8);
    }

    private static Path rotatePath(Path base, int index) {
        if (index == 0) {
            return base;
        }
        String fileName = base.getFileName().toString();
        int dot = fileName.lastIndexOf('.');
        String prefix = dot >= 0 ? fileName.substring(0, dot) : fileName;
        String suffix = dot >= 0 ? fileName.substring(dot) : "";
        String rotatedName = prefix + "_part" + index + suffix;
        return base.getParent() == null ? Path.of(rotatedName) : base.getParent().resolve(rotatedName);
    }

    @Override
    public synchronized void close() throws IOException {
        if (currentWriter != null) {
            currentWriter.close();
            currentWriter = null;
        }
    }

    static final class HitRecord {
        final String reactionId;
        final List<String> fragmentIds;
        final String assembledIdCode;
        final int atomCount;
        final int rotatableBonds;
        final double phesaSimilarity;
        final long attemptIndex;

        HitRecord(String reactionId,
                  List<String> fragmentIds,
                  String assembledIdCode,
                  int atomCount,
                  int rotatableBonds,
                  double phesaSimilarity,
                  long attemptIndex) {
            this.reactionId = reactionId;
            this.fragmentIds = fragmentIds;
            this.assembledIdCode = assembledIdCode;
            this.atomCount = atomCount;
            this.rotatableBonds = rotatableBonds;
            this.phesaSimilarity = phesaSimilarity;
            this.attemptIndex = attemptIndex;
        }
    }
}
