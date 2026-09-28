package com.idorsia.research.chem.hyperspace3d.workflow;

import java.net.URI;
import java.nio.file.*;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class FingerprintFileOpsTest {
    @TempDir Path temp;
    @Test void copyFallbackPreservesBytesWhenDestinationCannotHardLink() throws Exception {
        Path source = Files.writeString(temp.resolve("payload"), "immutable vector payload");
        try (var fs = FileSystems.newFileSystem(URI.create("jar:" + temp.resolve("destination.zip").toUri()), Map.of("create", "true"))) {
            Path target = fs.getPath("/payload");
            assertFalse(FingerprintFinalizer.linkOrCopy(source, target));
            assertArrayEquals(Files.readAllBytes(source), Files.readAllBytes(target));
        }
        assertEquals("immutable vector payload", Files.readString(source));
    }
}
