package com.idorsia.research.chem.hyperspace.io;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongPFP512;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerLongFFP1024_plus;
import com.idorsia.research.chem.hyperspace.descriptor.DescriptorHandlerPPCore;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthon;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.SynthonReactionValidator;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

/**
 * Parser capable of ingesting classic Enamine TSV files as well as directories containing
 * per-reaction CSV listings with optional metadata columns.
 */
public final class SynthonSpaceParser3 {

    private SynthonSpaceParser3() {
    }

    public static ParsedSpace parseEnamine(EnamineOptions options) throws Exception {
        Objects.requireNonNull(options, "options");
        DescriptorContext descriptor = DescriptorContext.forMode(options.mode());

        Map<String, ReactionRecord> records = parseEnamineFile(options, descriptor);
        List<ReactionRecord> filtered = filterReactions(records, options.maxSynthonSets());

        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(options.spaceName())
                .putMetadata(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT, "enamine-tsv")
                .putMetadata("source.file", options.input().toAbsolutePath().toString())
                .putMetadata("parser.mode", descriptor.modeLabel())
                .putMetadata("parser.maxSynthonSets", Integer.toString(options.maxSynthonSets()))
                .putMetadata("parser.threads", Integer.toString(options.threads()))
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_SHORT_NAME, descriptor.shortName())
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_BITS, Integer.toString(descriptor.bits()));
        addDescriptorTags(builder, descriptor, options.descriptorTags());

        validateAndApply(filtered, builder);
        RawSynthonSpace raw = builder.build();
        SynthonSpace synthonSpace = maybeBuildSynthonSpace(raw, descriptor, options.buildSynthonSpace());
        return new ParsedSpace(raw, synthonSpace);
    }

    public static ParsedSpace parseCsvDirectory(CsvDirectoryOptions options) throws Exception {
        Objects.requireNonNull(options, "options");
        if (!Files.isDirectory(options.directory())) {
            throw new IllegalArgumentException("Input is not a directory: " + options.directory());
        }
        DescriptorContext descriptor = DescriptorContext.forMode(options.mode());

        List<Path> reactionFiles = listReactionFiles(options);
        if (reactionFiles.isEmpty()) {
            throw new IllegalArgumentException("No reaction files found in " + options.directory());
        }

        List<ReactionRecord> records = parseCsvFiles(reactionFiles, options, descriptor);

        RawSynthonSpace.Builder builder = RawSynthonSpace.builder(options.spaceName())
                .putMetadata(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT, "csv-per-reaction")
                .putMetadata("source.directory", options.directory().toAbsolutePath().toString())
                .putMetadata("parser.mode", descriptor.modeLabel())
                .putMetadata("parser.smilesColumn", options.smilesColumn())
                .putMetadata("parser.idColumn", options.idColumn())
                .putMetadata("parser.priceColumn", Optional.ofNullable(options.priceColumn()).orElse(""))
                .putMetadata("parser.priceAttributeKey", Optional.ofNullable(options.priceAttributeKey()).orElse(""))
                .putMetadata("parser.synthonSetColumn", Optional.ofNullable(options.synthonSetColumn()).orElse(""))
                .putMetadata("parser.defaultSynthonSet", Integer.toString(options.defaultSynthonSet()))
                .putMetadata("parser.includeAllFiles", Boolean.toString(options.includeAllFiles()))
                .putMetadata("parser.threads", Integer.toString(options.threads()))
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_SHORT_NAME, descriptor.shortName())
                .putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_BITS, Integer.toString(descriptor.bits()));
        addDescriptorTags(builder, descriptor, options.descriptorTags());

        validateAndApply(records, builder);
        RawSynthonSpace raw = builder.build();
        SynthonSpace synthonSpace = maybeBuildSynthonSpace(raw, descriptor, options.buildSynthonSpace());
        return new ParsedSpace(raw, synthonSpace);
    }

    private static final int ENAMINE_BATCH_SIZE = 2_000;

    private static Map<String, ReactionRecord> parseEnamineFile(EnamineOptions options,
                                                                DescriptorContext descriptor) throws Exception {
        Map<String, ReactionRecord> records = new LinkedHashMap<>();
        Object recordLock = new Object();
        ProgressTracker progress = ProgressTracker.create(options.input());
        int threads = Math.max(1, options.threads());
        ExecutorService executor = threads > 1 ? Executors.newFixedThreadPool(threads) : null;
        List<Future<ChunkStats>> futures = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(options.input(), StandardCharsets.UTF_8)) {
            String header = reader.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Empty input file: " + options.input());
            }
            List<String> batch = new ArrayList<>(ENAMINE_BATCH_SIZE);
            long batchBytes = 0L;
            String line;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }
                batch.add(trimmed);
                batchBytes += line.length() + 1L;
                if (batch.size() >= ENAMINE_BATCH_SIZE) {
                    submitEnamineBatch(batch,
                            batchBytes,
                            options,
                            descriptor,
                            records,
                            recordLock,
                            executor,
                            futures,
                            progress);
                    batch = new ArrayList<>(ENAMINE_BATCH_SIZE);
                    batchBytes = 0L;
                }
            }
            if (!batch.isEmpty()) {
                submitEnamineBatch(batch,
                        batchBytes,
                        options,
                        descriptor,
                        records,
                        recordLock,
                        executor,
                        futures,
                        progress);
            }
            if (executor != null) {
                executor.shutdown();
                for (Future<ChunkStats> future : futures) {
                    ChunkStats stats = future.get();
                    progress.increment(stats);
                }
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
            }
        } catch (IOException e) {
            if (executor != null) {
                executor.shutdownNow();
            }
            throw new IllegalStateException("Unable to parse file " + options.input(), e);
        } catch (Exception e) {
            if (executor != null) {
                executor.shutdownNow();
            }
            throw e;
        }
        progress.finish();
        return records;
    }

    private static void submitEnamineBatch(List<String> batch,
                                           long batchBytes,
                                           EnamineOptions options,
                                           DescriptorContext descriptor,
                                           Map<String, ReactionRecord> records,
                                           Object recordLock,
                                           ExecutorService executor,
                                           List<Future<ChunkStats>> futures,
                                           ProgressTracker progress) throws Exception {
        EnamineChunkTask task = new EnamineChunkTask(new ArrayList<>(batch),
                batchBytes,
                options.input(),
                descriptor,
                records,
                recordLock);
        if (executor == null) {
            ChunkStats stats = task.call();
            progress.increment(stats);
        } else {
            futures.add(executor.submit(task));
        }
    }

    private static List<Path> listReactionFiles(CsvDirectoryOptions options) throws IOException {
        try (var stream = Files.list(options.directory())) {
            return stream.filter(Files::isRegularFile)
                    .filter(path -> options.includeAllFiles()
                            || path.getFileName().toString().toLowerCase(Locale.ROOT).endsWith(".csv"))
                    .sorted(Comparator.comparing(path -> path.getFileName().toString().toLowerCase(Locale.ROOT)))
                    .collect(Collectors.toList());
        }
    }

    private static List<ReactionRecord> parseCsvFiles(List<Path> files,
                                                      CsvDirectoryOptions options,
                                                      DescriptorContext descriptor) throws Exception {
        List<ReactionRecord> records = new ArrayList<>();
        if (options.threads() <= 1) {
            SmilesParser parser = new SmilesParser();
            for (Path file : files) {
                ReactionRecord record = parseCsvFile(file, options, descriptor, parser);
                if (record != null) {
                    records.add(record);
                }
            }
        } else {
            ExecutorService executor = Executors.newFixedThreadPool(options.threads());
            List<Future<ReactionRecord>> futures = new ArrayList<>();
            for (Path file : files) {
                futures.add(executor.submit(() -> parseCsvFile(file, options, descriptor, new SmilesParser())));
            }
            executor.shutdown();
            for (Future<ReactionRecord> future : futures) {
                ReactionRecord record;
                try {
                    record = future.get();
                } catch (ExecutionException e) {
                    throw new IllegalStateException("Failed to parse reaction CSV file", e.getCause());
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new IllegalStateException("Parsing interrupted", e);
                }
                if (record != null) {
                    records.add(record);
                }
            }
        }
        return records;
    }

    private static ReactionRecord parseCsvFile(Path file,
                                               CsvDirectoryOptions options,
                                               DescriptorContext descriptor,
                                               SmilesParser parser) throws Exception {
        String reactionId = reactionIdFromFile(file);
        ReactionRecord record = new ReactionRecord(reactionId);
        List<String> headerTokens;
        try (BufferedReader reader = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IllegalArgumentException("Empty reaction file: " + file);
            }
            headerTokens = splitCsvLine(headerLine);
            Map<String, Integer> headerIndex = buildHeaderIndex(headerTokens);
            ColumnMapping mapping = ColumnMapping.fromHeader(headerIndex, options);

            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                List<String> parts = splitCsvLine(line);
                if (parts.isEmpty()) {
                    continue;
                }
                while (parts.size() < headerTokens.size()) {
                    parts.add("");
                }
                String smiles = parts.get(mapping.smilesIdx()).trim();
                if (smiles.isEmpty()) {
                    continue;
                }
                StereoMolecule molecule = parseSmiles(parser, smiles, file.toString());
                String fragmentId = parts.get(mapping.idIdx()).trim();
                int synthonIdx = mapping.synthonIdx(parts);

                RawSynthon fragment = createFragment(reactionId, synthonIdx, fragmentId, molecule);
                record.addFragment(synthonIdx, fragment);
                mapping.priceValue(parts).ifPresent(price -> record.addFragmentAttribute(fragmentId,
                        options.priceAttributeKey(), price));
            }
        } catch (IOException e) {
            throw new IllegalStateException("Unable to parse reaction file " + file, e);
        }
        if (record.isEmpty()) {
            System.out.println("Skipping empty reaction file: " + file);
            return null;
        }
        return record;
    }

    private static StereoMolecule parseSmiles(SmilesParser parser, String smiles, String source) throws Exception {
        try {
            StereoMolecule molecule = new StereoMolecule();
            parser.parse(molecule, smiles);
            molecule.ensureHelperArrays(Molecule.cHelperCIP);
            return molecule;
        } catch (Exception e) {
            throw new IllegalArgumentException("Failed to parse SMILES '" + smiles + "' from " + source, e);
        }
    }

    private static List<ReactionRecord> filterReactions(Map<String, ReactionRecord> records, int maxSynthonSets) {
        if (maxSynthonSets <= 0) {
            return new ArrayList<>(records.values());
        }
        List<ReactionRecord> filtered = new ArrayList<>();
        records.forEach((rxnId, record) -> {
            if (record.fragmentSetCount() > maxSynthonSets) {
                System.out.println("Skipping reaction " + rxnId + " with " + record.fragmentSetCount()
                        + " synthon sets (max=" + maxSynthonSets + ")");
            } else {
                filtered.add(record);
            }
        });
        return filtered;
    }

    private static void validateAndApply(List<ReactionRecord> records,
                                         RawSynthonSpace.Builder builder) {
        for (ReactionRecord record : records) {
            try {
                record.validate();
            } catch (Exception ex) {
                System.out.println("Skipping reaction " + record.reactionId()
                        + " due to validation failure: " + ex.getMessage());
                continue;
            }
            record.apply(builder);
        }
    }

    private static SynthonSpace maybeBuildSynthonSpace(RawSynthonSpace raw,
                                                       DescriptorContext descriptor,
                                                       boolean build) throws Exception {
        if (!build) {
            return null;
        }
        RawSynthonSpaceAssembler.BuildOptions buildOptions = RawSynthonSpaceAssembler.BuildOptions.builder()
                .descriptorShortName(descriptor.shortName())
                .descriptorBits(descriptor.bits())
                .build();
        return RawSynthonSpaceAssembler.buildSynthonSpace(raw, buildOptions);
    }

    private static void addDescriptorTags(RawSynthonSpace.Builder builder,
                                          DescriptorContext descriptor,
                                          List<String> requestedTags) {
        LinkedHashSet<String> tags = new LinkedHashSet<>();
        if (descriptor != null && descriptor.shortName() != null) {
            tags.add(descriptor.shortName());
        }
        if (requestedTags != null) {
            requestedTags.stream()
                    .filter(Objects::nonNull)
                    .map(String::trim)
                    .filter(token -> !token.isEmpty())
                    .forEach(tags::add);
        }
        if (!tags.isEmpty()) {
            builder.putMetadata(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS, String.join(",", tags));
        }
    }

    private static RawSynthon createFragment(String reactionId,
                                             int synthonIdx,
                                             String fragmentId,
                                             StereoMolecule molecule) {
        StereoMolecule working = new StereoMolecule(molecule);
        working.ensureHelperArrays(Molecule.cHelperCIP);
        return RawSynthon.fromMolecule(reactionId, synthonIdx, fragmentId, working);
    }

    private static String reactionIdFromFile(Path file) {
        String name = file.getFileName().toString();
        int idx = name.lastIndexOf('.');
        if (idx > 0) {
            return name.substring(0, idx);
        }
        return name;
    }

    private static Map<String, Integer> buildHeaderIndex(List<String> headerTokens) {
        Map<String, Integer> index = new HashMap<>();
        for (int i = 0; i < headerTokens.size(); i++) {
            String token = headerTokens.get(i);
            if (token == null) {
                continue;
            }
            index.put(token.trim().toLowerCase(Locale.ROOT), i);
        }
        return index;
    }

    private static List<String> splitCsvLine(String line) {
        List<String> tokens = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    current.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (ch == ',' && !inQuotes) {
                tokens.add(current.toString());
                current.setLength(0);
            } else {
                current.append(ch);
            }
        }
        tokens.add(current.toString());
        return tokens;
    }

    public static final class ParsedSpace {
        private final RawSynthonSpace rawSpace;
        private final SynthonSpace synthonSpace;

        public ParsedSpace(RawSynthonSpace rawSpace, SynthonSpace synthonSpace) {
            this.rawSpace = rawSpace;
            this.synthonSpace = synthonSpace;
        }

        public RawSynthonSpace rawSpace() {
            return rawSpace;
        }

        public SynthonSpace synthonSpace() {
            return synthonSpace;
        }
    }

    public static final class EnamineOptions {
        private final Path input;
        private final String spaceName;
        private final String mode;
        private final int threads;
        private final int maxSynthonSets;
        private final List<String> descriptorTags;
        private final boolean buildSynthonSpace;

        private EnamineOptions(Builder builder) {
            this.input = Objects.requireNonNull(builder.input, "input");
            this.spaceName = Objects.requireNonNull(builder.spaceName, "spaceName");
            this.mode = builder.mode == null ? "FragFp" : builder.mode;
            this.threads = builder.threads;
            this.maxSynthonSets = builder.maxSynthonSets;
            this.descriptorTags = List.copyOf(builder.descriptorTags);
            this.buildSynthonSpace = builder.buildSynthonSpace;
        }

        public Path input() {
            return input;
        }

        public String spaceName() {
            return spaceName;
        }

        public String mode() {
            return mode;
        }

        public int threads() {
            return threads;
        }

        public int maxSynthonSets() {
            return maxSynthonSets;
        }

        public List<String> descriptorTags() {
            return descriptorTags;
        }

        public boolean buildSynthonSpace() {
            return buildSynthonSpace;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static final class Builder {
            private Path input;
            private String spaceName;
            private String mode = "FragFp";
            private int threads = 4;
            private int maxSynthonSets = 3;
            private final LinkedHashSet<String> descriptorTags = new LinkedHashSet<>();
            private boolean buildSynthonSpace = false;

            public Builder input(Path input) {
                this.input = input;
                return this;
            }

            public Builder spaceName(String spaceName) {
                this.spaceName = spaceName;
                return this;
            }

            public Builder mode(String mode) {
                this.mode = mode;
                return this;
            }

            public Builder threads(int threads) {
                this.threads = threads;
                return this;
            }

            public Builder maxSynthonSets(int maxSynthonSets) {
                this.maxSynthonSets = maxSynthonSets;
                return this;
            }

            public Builder addDescriptorTag(String tag) {
                if (tag != null) {
                    String trimmed = tag.trim();
                    if (!trimmed.isEmpty()) {
                        descriptorTags.add(trimmed);
                    }
                }
                return this;
            }

            public Builder buildSynthonSpace(boolean build) {
                this.buildSynthonSpace = build;
                return this;
            }

            public EnamineOptions build() {
                return new EnamineOptions(this);
            }
        }
    }

    public static final class CsvDirectoryOptions {
        private final Path directory;
        private final String spaceName;
        private final String mode;
        private final String smilesColumn;
        private final String idColumn;
        private final String priceColumn;
        private final String priceAttributeKey;
        private final String synthonSetColumn;
        private final int defaultSynthonSet;
        private final boolean includeAllFiles;
        private final int threads;
        private final List<String> descriptorTags;
        private final boolean buildSynthonSpace;

        private CsvDirectoryOptions(Builder builder) {
            this.directory = Objects.requireNonNull(builder.directory, "directory");
            this.spaceName = Objects.requireNonNull(builder.spaceName, "spaceName");
            this.mode = builder.mode == null ? "FragFp" : builder.mode;
            this.smilesColumn = builder.smilesColumn;
            this.idColumn = builder.idColumn;
            this.priceColumn = builder.priceColumn;
            this.priceAttributeKey = builder.priceAttributeKey;
            this.synthonSetColumn = builder.synthonSetColumn;
            this.defaultSynthonSet = builder.defaultSynthonSet;
            this.includeAllFiles = builder.includeAllFiles;
            this.threads = builder.threads;
            this.descriptorTags = List.copyOf(builder.descriptorTags);
            this.buildSynthonSpace = builder.buildSynthonSpace;
        }

        public Path directory() {
            return directory;
        }

        public String spaceName() {
            return spaceName;
        }

        public String mode() {
            return mode;
        }

        public String smilesColumn() {
            return smilesColumn;
        }

        public String idColumn() {
            return idColumn;
        }

        public String priceColumn() {
            return priceColumn;
        }

        public String priceAttributeKey() {
            return priceAttributeKey;
        }

        public String synthonSetColumn() {
            return synthonSetColumn;
        }

        public int defaultSynthonSet() {
            return defaultSynthonSet;
        }

        public boolean includeAllFiles() {
            return includeAllFiles;
        }

        public int threads() {
            return threads;
        }

        public List<String> descriptorTags() {
            return descriptorTags;
        }

        public boolean buildSynthonSpace() {
            return buildSynthonSpace;
        }

        public static Builder builder() {
            return new Builder();
        }

        public static final class Builder {
            private Path directory;
            private String spaceName;
            private String mode = "FragFp";
            private String smilesColumn = "SMILES";
            private String idColumn = "bb1_parent_id";
            private String priceColumn = "Price";
            private String priceAttributeKey = "price";
            private String synthonSetColumn;
            private int defaultSynthonSet = 0;
            private boolean includeAllFiles = false;
            private int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
            private final LinkedHashSet<String> descriptorTags = new LinkedHashSet<>();
            private boolean buildSynthonSpace = false;

            public Builder directory(Path directory) {
                this.directory = directory;
                return this;
            }

            public Builder spaceName(String spaceName) {
                this.spaceName = spaceName;
                return this;
            }

            public Builder mode(String mode) {
                this.mode = mode;
                return this;
            }

            public Builder smilesColumn(String smilesColumn) {
                this.smilesColumn = smilesColumn;
                return this;
            }

            public Builder idColumn(String idColumn) {
                this.idColumn = idColumn;
                return this;
            }

            public Builder priceColumn(String priceColumn) {
                this.priceColumn = priceColumn;
                return this;
            }

            public Builder priceAttributeKey(String priceAttributeKey) {
                this.priceAttributeKey = priceAttributeKey;
                return this;
            }

            public Builder synthonSetColumn(String synthonSetColumn) {
                this.synthonSetColumn = synthonSetColumn;
                return this;
            }

            public Builder defaultSynthonSet(int defaultSynthonSet) {
                this.defaultSynthonSet = defaultSynthonSet;
                return this;
            }

            public Builder includeAllFiles(boolean includeAllFiles) {
                this.includeAllFiles = includeAllFiles;
                return this;
            }

            public Builder threads(int threads) {
                if (threads > 0) {
                    this.threads = threads;
                }
                return this;
            }

            public Builder addDescriptorTag(String tag) {
                if (tag != null) {
                    String trimmed = tag.trim();
                    if (!trimmed.isEmpty()) {
                        descriptorTags.add(trimmed);
                    }
                }
                return this;
            }

            public Builder buildSynthonSpace(boolean build) {
                this.buildSynthonSpace = build;
                return this;
            }

            public CsvDirectoryOptions build() {
                return new CsvDirectoryOptions(this);
            }
        }
    }

    private static ReactionRecord getOrCreateRecord(Map<String, ReactionRecord> records,
                                                    Object lock,
                                                    String reactionId) {
        synchronized (lock) {
            return records.computeIfAbsent(reactionId, ReactionRecord::new);
        }
    }

    private static final class EnamineChunkTask implements Callable<ChunkStats> {
        private final List<String> lines;
        private final long bytes;
        private final Path source;
        private final DescriptorContext descriptor;
        private final Map<String, ReactionRecord> records;
        private final Object recordLock;

        private EnamineChunkTask(List<String> lines,
                                 long bytes,
                                 Path source,
                                 DescriptorContext descriptor,
                                 Map<String, ReactionRecord> records,
                                 Object recordLock) {
            this.lines = lines;
            this.bytes = bytes;
            this.source = source;
            this.descriptor = descriptor;
            this.records = records;
            this.recordLock = recordLock;
        }

        @Override
        public ChunkStats call() throws Exception {
            SmilesParser parser = new SmilesParser();
            long fragments = 0L;
            for (String line : lines) {
                if (line.isEmpty()) {
                    continue;
                }
                String[] tokens = line.split("\\t");
                if (tokens.length < 4) {
                    throw new IllegalArgumentException("Malformed TSV line: " + line);
                }
                String smiles = tokens[0];
                String fragmentId = tokens[1];
                int synthonIdx = Integer.parseInt(tokens[2]);
                String reactionId = tokens[3];

                StereoMolecule molecule = parseSmiles(parser, smiles, source.toString());
                RawSynthon fragment = createFragment(reactionId, synthonIdx, fragmentId, molecule);
                ReactionRecord record = getOrCreateRecord(records, recordLock, reactionId);
                record.addFragment(synthonIdx, fragment);
                fragments++;
            }
            return new ChunkStats(bytes, fragments);
        }
    }

    private static final class ChunkStats {
        private final long bytes;
        private final long fragments;

        private ChunkStats(long bytes, long fragments) {
            this.bytes = bytes;
            this.fragments = fragments;
        }

        long bytes() {
            return bytes;
        }

        long fragments() {
            return fragments;
        }
    }

    private static final class ProgressTracker {
        private final long totalBytes;
        private final AtomicLong processedBytes = new AtomicLong();
        private final AtomicLong processedFragments = new AtomicLong();
        private final long startNanos = System.nanoTime();
        private volatile long lastLog = startNanos;

        private ProgressTracker(long totalBytes) {
            this.totalBytes = totalBytes;
        }

        static ProgressTracker create(Path input) {
            long total = 0L;
            try {
                total = Files.size(input);
            } catch (IOException ignored) {
            }
            return new ProgressTracker(total);
        }

        void increment(ChunkStats stats) {
            if (stats == null) {
                return;
            }
            processedBytes.addAndGet(stats.bytes());
            processedFragments.addAndGet(stats.fragments());
            maybeLog();
        }

        private void maybeLog() {
            long now = System.nanoTime();
            if (TimeUnit.NANOSECONDS.toMillis(now - lastLog) < 5000) {
                return;
            }
            lastLog = now;
            long fragments = processedFragments.get();
            if (totalBytes > 0) {
                long percent = Math.min(100, (processedBytes.get() * 100) / totalBytes);
                System.out.println("[EnamineImport] " + percent + "% processed (" + fragments + " fragments)");
            } else {
                System.out.println("[EnamineImport] processed " + fragments + " fragments");
            }
        }

        void finish() {
            long fragments = processedFragments.get();
            if (fragments == 0) {
                return;
            }
            if (totalBytes > 0) {
                processedBytes.set(totalBytes);
                System.out.println("[EnamineImport] 100% processed (" + fragments + " fragments)");
            } else {
                System.out.println("[EnamineImport] processed " + fragments + " fragments");
            }
        }
    }

    private static final class ReactionRecord {
        private final String reactionId;
        private final Map<Integer, List<RawSynthon>> fragmentSets = new LinkedHashMap<>();
        private final Map<String, Map<String, String>> fragmentAttributes = new LinkedHashMap<>();

        private ReactionRecord(String reactionId) {
            this.reactionId = reactionId;
        }

        synchronized void addFragment(int synthonIdx, RawSynthon fragment) {
            fragmentSets.computeIfAbsent(synthonIdx, key -> new ArrayList<>()).add(fragment);
        }

        synchronized void addFragmentAttribute(String fragmentId, String key, String value) {
            if (fragmentId == null || fragmentId.isBlank() || key == null || key.isBlank() || value == null) {
                return;
            }
            fragmentAttributes.computeIfAbsent(fragmentId, id -> new LinkedHashMap<>()).put(key, value);
        }

        boolean isEmpty() {
            return fragmentSets.isEmpty();
        }

        int fragmentSetCount() {
            return fragmentSets.size();
        }

        void validate() throws Exception {
            Map<Integer, List<Object>> molecules = new LinkedHashMap<>();
            fragmentSets.forEach((idx, list) -> {
                List<Object> idcodes = new ArrayList<>(list.size());
                for (RawSynthon frag : list) {
                    idcodes.add(frag.getIdcode());
                }
                molecules.put(idx, idcodes);
            });
            SynthonReactionValidator.validate(molecules);
        }

        void apply(RawSynthonSpace.Builder builder) {
            fragmentSets.forEach((idx, list) -> builder.addRawFragments(reactionId, idx, list));
            fragmentAttributes.forEach((fragmentId, attributes) ->
                    attributes.forEach((key, value) -> builder.addFragmentAttribute(reactionId, fragmentId, key, value)));
        }

        String reactionId() {
            return reactionId;
        }
    }

    private static final class DescriptorContext {
        private final DescriptorHandler<long[], StereoMolecule> prototype;
        private final ThreadLocal<DescriptorHandler<long[], StereoMolecule>> handlers;
        private final int bits;
        private final String shortName;
        private final String modeLabel;

        private DescriptorContext(DescriptorHandler<long[], StereoMolecule> prototype,
                                  int bits,
                                  String shortName,
                                  String modeLabel) {
            this.prototype = prototype;
            this.handlers = ThreadLocal.withInitial(prototype::getThreadSafeCopy);
            this.bits = bits;
            this.shortName = shortName;
            this.modeLabel = modeLabel;
        }

        static DescriptorContext forMode(String mode) {
            String normalized = mode == null || mode.isBlank() ? "FragFp" : mode;
            switch (normalized.toLowerCase(Locale.ROOT)) {
                case "mode_pfp":
                case "pfp":
                    DescriptorHandlerLongFFP1024_plus pfp = new DescriptorHandlerLongFFP1024_plus("pfp");
                    return new DescriptorContext(pfp,
                            1024,
                            pfp.getInfo().shortName,
                            normalized);
                case "mode_ffp":
                case "ffp":
                    DescriptorHandlerLongFFP1024_plus ffp = new DescriptorHandlerLongFFP1024_plus("ffp");
                    return new DescriptorContext(ffp,
                            1024,
                            ffp.getInfo().shortName,
                            normalized);
                case "mode_ppcore":
                case "ppc":
                    DescriptorHandlerPPCore ppcore = new DescriptorHandlerPPCore();
                    return new DescriptorContext(ppcore,
                            2048,
                            ppcore.getInfo().shortName,
                            normalized);
                case "pathfp":
                    DescriptorHandlerLongPFP512 path = new DescriptorHandlerLongPFP512();
                    return new DescriptorContext(path,
                            512,
                            path.getInfo().shortName,
                            normalized);
                case "fragfp":
                    DescriptorHandlerLongFFP512 frag = new DescriptorHandlerLongFFP512();
                    return new DescriptorContext(frag,
                            512,
                            frag.getInfo().shortName,
                            normalized);
                default:
                    throw new IllegalArgumentException("Descriptor mode not supported for raw import: " + mode);
            }
        }

        int bits() {
            return bits;
        }

        String shortName() {
            return shortName;
        }

        String modeLabel() {
            return modeLabel;
        }
    }

    private static final class ColumnMapping {
        private final int smilesIdx;
        private final int idIdx;
        private final Integer synthonIdxColumn;
        private final int defaultSynthonIdx;
        private final Integer priceIdx;

        private ColumnMapping(int smilesIdx,
                              int idIdx,
                              Integer synthonIdxColumn,
                              int defaultSynthonIdx,
                              Integer priceIdx) {
            this.smilesIdx = smilesIdx;
            this.idIdx = idIdx;
            this.synthonIdxColumn = synthonIdxColumn;
            this.defaultSynthonIdx = defaultSynthonIdx;
            this.priceIdx = priceIdx;
        }

        static ColumnMapping fromHeader(Map<String, Integer> headerIndex, CsvDirectoryOptions options) {
            int smiles = requireColumn(headerIndex, options.smilesColumn(), "SMILES");
            int idCol = requireColumn(headerIndex, options.idColumn(), "fragment id");
            Integer synthonIdxColumn = resolveOptionalColumn(headerIndex, options.synthonSetColumn());
            Integer price = resolveOptionalColumn(headerIndex, options.priceColumn());
            return new ColumnMapping(smiles, idCol, synthonIdxColumn, options.defaultSynthonSet(), price);
        }

        int smilesIdx() {
            return smilesIdx;
        }

        int idIdx() {
            return idIdx;
        }

        int synthonIdx(List<String> parts) {
            if (synthonIdxColumn == null) {
                return defaultSynthonIdx;
            }
            String token = parts.get(synthonIdxColumn).trim();
            if (token.isEmpty()) {
                return defaultSynthonIdx;
            }
            try {
                return Integer.parseInt(token);
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException("Invalid synthon index: " + token, ex);
            }
        }

        Optional<String> priceValue(List<String> parts) {
            if (priceIdx == null) {
                return Optional.empty();
            }
            String value = parts.get(priceIdx).trim();
            if (value.isEmpty()) {
                return Optional.empty();
            }
            return Optional.of(value);
        }

        private static int requireColumn(Map<String, Integer> headerIndex, String name, String label) {
            if (name == null) {
                throw new IllegalArgumentException("Missing configuration for " + label + " column");
            }
            Integer idx = headerIndex.get(name.toLowerCase(Locale.ROOT));
            if (idx == null) {
                throw new IllegalArgumentException("Column \"" + name + "\" not found in CSV header");
            }
            return idx;
        }

        private static Integer resolveOptionalColumn(Map<String, Integer> headerIndex, String columnName) {
            if (columnName == null) {
                return null;
            }
            return headerIndex.get(columnName.toLowerCase(Locale.ROOT));
        }
    }
}
