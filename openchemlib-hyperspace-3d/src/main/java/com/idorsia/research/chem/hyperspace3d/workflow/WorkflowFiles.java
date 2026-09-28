package com.idorsia.research.chem.hyperspace3d.workflow;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import java.io.*;
import java.nio.channels.*;
import java.nio.file.*;
import java.security.*;
import java.util.*;

/** Atomic shared-filesystem publication; locks must be supported by the cluster filesystem. */
public final class WorkflowFiles {
    public static final ObjectMapper JSON = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT)
            .enable(SerializationFeature.ORDER_MAP_ENTRIES_BY_KEYS);
    private WorkflowFiles() {}

    public static String hash(Path file) throws IOException {
        try (InputStream in = Files.newInputStream(file)) {
            MessageDigest digest = digest();
            byte[] buffer = new byte[1 << 20];
            for (int n; (n = in.read(buffer)) != -1;) digest.update(buffer, 0, n);
            return HexFormat.of().formatHex(digest.digest());
        }
    }
    public static MessageDigest digest() {
        try { return MessageDigest.getInstance("SHA-256"); }
        catch (NoSuchAlgorithmException e) { throw new IllegalStateException(e); }
    }
    public static String identity(Object value) throws IOException {
        return HexFormat.of().formatHex(digest().digest(JSON.writeValueAsBytes(value)));
    }
    public static void json(Path path, Object value) throws IOException {
        Files.createDirectories(path.toAbsolutePath().getParent());
        Path tmp = Files.createTempFile(path.getParent(), ".json-", ".partial");
        try {
            JSON.writeValue(tmp.toFile(), value);
            Files.move(tmp, path, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } finally { Files.deleteIfExists(tmp); }
    }
    public static Path child(Path root, String relative) throws IOException {
        Path p = root.resolve(relative).normalize();
        if (Path.of(relative).isAbsolute() || !p.startsWith(root.normalize()) || p.equals(root))
            throw new IOException("unsafe relative path: " + relative);
        return p;
    }
    public static Map<String, String> hashes(Path root) throws IOException {
        Map<String, String> hashes = new TreeMap<>();
        try (var files = Files.walk(root)) {
            for (Path p : files.sorted().toList()) {
                if (Files.isSymbolicLink(p)) throw new IOException("symbolic links are not supported: " + p);
                if (Files.isRegularFile(p) && !p.getFileName().toString().equals("completion.json"))
                    hashes.put(root.relativize(p).toString(), hash(p));
            }
        }
        return hashes;
    }
    public static void verify(Path root, Map<String, String> expected) throws IOException {
        if (!hashes(root).equals(expected)) throw new IOException("checksum mismatch: " + root);
    }
    public static void copyTree(Path from, Path to) throws IOException {
        Files.createDirectories(to);
        try (var files = Files.walk(from)) {
            for (Path p : files.toList()) {
                if (Files.isSymbolicLink(p)) throw new IOException("symbolic links are not supported: " + p);
                Path out = to.resolve(from.relativize(p));
                if (Files.isDirectory(p)) Files.createDirectories(out);
                else Files.copy(p, out);
            }
        }
    }
    public static void deleteTree(Path path) throws IOException {
        if (!Files.exists(path, LinkOption.NOFOLLOW_LINKS)) return;
        try (var files = Files.walk(path)) {
            for (Path p : files.sorted(Comparator.reverseOrder()).toList()) Files.delete(p);
        }
    }
    public static Locked lock(Path path) throws IOException { return new Locked(path); }
    public static final class Locked implements AutoCloseable {
        private final FileChannel channel;
        private final FileLock lock;
        Locked(Path path) throws IOException {
            Files.createDirectories(path.toAbsolutePath().getParent());
            channel = FileChannel.open(path, StandardOpenOption.CREATE, StandardOpenOption.WRITE);
            try {
                lock = channel.tryLock();
                if (lock == null) throw new IOException("workflow is busy: " + path);
            } catch (IOException | RuntimeException e) { channel.close(); throw e; }
        }
        @Override public void close() throws IOException { try { lock.release(); } finally { channel.close(); } }
    }
}
