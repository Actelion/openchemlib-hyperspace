package com.idorsia.research.chem.hyperspace.io;

import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser3;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import static org.junit.jupiter.api.Assertions.*;

class RawSynthonSpaceImporterTest {

    @Test
    void importerCreatesRawSpaceWithMetadata() throws Exception {
        Path sample = Paths.get("..", "openchemlib-hyperspace-core", "src", "main", "resources",
                "testdata", "idorsia_toy_space_a.txt").toAbsolutePath().normalize();
        assertTrue(sample.toFile().exists(), "Sample TSV missing: " + sample);

        RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
        SynthonSpaceParser3.EnamineOptions options = SynthonSpaceParser3.EnamineOptions.builder()
                .input(sample)
                .spaceName("toy")
                .mode("FragFp")
                .threads(2)
                .maxSynthonSets(3)
                .buildSynthonSpace(true)
                .addDescriptorTag("FragFp")
                .build();
        RawSynthonSpaceImporter.Result result = importer.importEnamine(options);
        RawSynthonSpace raw = result.getRawSpace();

        assertEquals("toy", raw.getName());
        assertEquals("enamine-tsv", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
        assertEquals("FragFp", raw.getMetadata().get("parser.mode"));
        assertEquals("FragFp", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS));
        assertTrue(raw.getReactions().containsKey("benzoimidazole_b-8"));
        assertNotNull(result.getSynthonSpace());
        assertFalse(result.getSynthonSpace().getRxnIds().isEmpty());
    }

    @Test
    void importerParsesCsvDirectoryAndStoresPrices() throws Exception {
        Path tempDir = Files.createTempDirectory("synthon_csv");
        Path csv = tempDir.resolve("benzoimidazole_b-8.csv");
        String content = String.join("\n",
                "SMILES,bb1_parent_id,Price,SynthonSet",
                "[U]c([nH]nc1)c1-c(nn1)c[n]1[Np],015119586-618282196,100,0",
                "N#Cc1ccc(C(C[Np])=O)[s]1,956652961-712098415,200,1",
                "Clc(c(Cl)c1)cc([nH]2)c1nc2[U],558027358-494530800,300,2");
        Files.writeString(csv, content, StandardCharsets.UTF_8);

        try {
            RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
            SynthonSpaceParser3.CsvDirectoryOptions options = SynthonSpaceParser3.CsvDirectoryOptions.builder()
                    .directory(tempDir)
                    .spaceName("csvSpace")
                    .mode("FragFp")
                    .smilesColumn("SMILES")
                    .idColumn("bb1_parent_id")
                    .priceColumn("Price")
                    .priceAttributeKey("price.usd")
                    .synthonSetColumn("SynthonSet")
                    .threads(2)
                    .buildSynthonSpace(true)
                    .addDescriptorTag("FragFp")
                    .build();
            RawSynthonSpace raw = importer.importCsvDirectory(options).getRawSpace();
            assertEquals("csv-per-reaction", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
            assertTrue(raw.getReactions().containsKey("benzoimidazole_b-8"));
            assertEquals("FragFp", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.DESCRIPTOR_TAGS));
            Map<String, Map<String, String>> attributes = raw.getReactions()
                    .get("benzoimidazole_b-8").getFragmentAttributes();
            assertEquals("100", attributes.get("015119586-618282196").get("price.usd"));
            assertEquals(3, raw.getFragmentSets("benzoimidazole_b-8").size());
        } finally {
            Files.walk(tempDir)
                    .sorted(Comparator.reverseOrder())
                    .forEach(path -> {
                        try {
                            Files.deleteIfExists(path);
                        } catch (IOException ignored) {
                        }
                    });
        }
    }

    @Test
    void importerParsesXtalpiStyleCsvRoles() throws Exception {
        Path tempDir = Files.createTempDirectory("xtalpi_csv");
        Path csv = tempDir.resolve("VAST_synthon_2026_H1.csv");
        String content = String.join("\n",
                "SMILES,synton_id,synton_role,reaction_id",
                "[U]c([nH]nc1)c1-c(nn1)c[n]1[Np],015119586-618282196,synton_0,benzoimidazole_b-8",
                "N#Cc1ccc(C(C[Np])=O)[s]1,956652961-712098415,synton_1,benzoimidazole_b-8",
                "Clc(c(Cl)c1)cc([nH]2)c1nc2[U],558027358-494530800,synton_2,benzoimidazole_b-8");
        Files.writeString(csv, content, StandardCharsets.UTF_8);

        try {
            RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
            SynthonSpaceParser3.EnamineOptions options = SynthonSpaceParser3.EnamineOptions.builder()
                    .input(csv)
                    .spaceName("xtalpi")
                    .mode("FragFp")
                    .sourceFormat("xtalpi-csv")
                    .idColumn("synton_id")
                    .synthonSetColumn("synton_role")
                    .putMetadata("source.supplier", "Xtalpi")
                    .buildSynthonSpace(false)
                    .build();
            RawSynthonSpace raw = importer.importEnamine(options).getRawSpace();
            assertEquals("xtalpi-csv", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
            assertEquals("Xtalpi", raw.getMetadata().get("source.supplier"));
            assertEquals(3, raw.getFragmentSets("benzoimidazole_b-8").size());
        } finally {
            Files.walk(tempDir)
                    .sorted(Comparator.reverseOrder())
                    .forEach(path -> {
                        try {
                            Files.deleteIfExists(path);
                        } catch (IOException ignored) {
                        }
                    });
        }
    }

    @Test
    void importerParsesZipTableAndReactionMetadata() throws Exception {
        Path tempDir = Files.createTempDirectory("molecule_one_zip");
        Path zip = tempDir.resolve("space.zip");
        try (ZipOutputStream out = new ZipOutputStream(Files.newOutputStream(zip))) {
            out.putNextEntry(new ZipEntry("synthons.txt"));
            out.write(String.join("\n",
                    "SMILES\tsynthon_id\tsynthon#\treaction_id",
                    "[U]c([nH]nc1)c1-c(nn1)c[n]1[Np]\t015119586-618282196\t0\tbenzoimidazole_b-8",
                    "N#Cc1ccc(C(C[Np])=O)[s]1\t956652961-712098415\t1\tbenzoimidazole_b-8",
                    "Clc(c(Cl)c1)cc([nH]2)c1nc2[U]\t558027358-494530800\t2\tbenzoimidazole_b-8")
                    .getBytes(StandardCharsets.UTF_8));
            out.closeEntry();
            out.putNextEntry(new ZipEntry("reactions.txt"));
            out.write(String.join("\n",
                    "reaction_id\tcomponents\tReaction\tProduct\tR1\tR2\tR3\tR4",
                    "benzoimidazole_b-8\t3\t[*:1].[*:2].[*:3]\t[*:1]-[*:2]-[*:3]\tR1\tR2\tR3\t-")
                    .getBytes(StandardCharsets.UTF_8));
            out.closeEntry();
        }

        try {
            RawSynthonSpaceImporter importer = new RawSynthonSpaceImporter();
            SynthonSpaceParser3.EnamineOptions options = SynthonSpaceParser3.EnamineOptions.builder()
                    .input(zip)
                    .spaceName("molecule_one")
                    .mode("FragFp")
                    .sourceFormat("molecule-one-zip")
                    .zipEntry("synthons.txt")
                    .reactionZipEntry("reactions.txt")
                    .putMetadata("source.supplier", "Molecule.One")
                    .buildSynthonSpace(false)
                    .build();
            RawSynthonSpace raw = importer.importEnamine(options).getRawSpace();
            assertEquals("molecule-one-zip", raw.getMetadata().get(RawSynthonSpace.MetadataKeys.SOURCE_FORMAT));
            assertEquals("synthons.txt", raw.getMetadata().get("source.zipEntry"));
            assertEquals("reactions.txt", raw.getMetadata().get("source.reactionZipEntry"));
            assertEquals("Molecule.One", raw.getMetadata().get("source.supplier"));
            assertEquals(3, raw.getFragmentSets("benzoimidazole_b-8").size());
            assertEquals("3", raw.getReactions().get("benzoimidazole_b-8")
                    .getReactionMetadata().get("supplier.reaction.components"));
            assertEquals("[*:1].[*:2].[*:3]", raw.getReactions().get("benzoimidazole_b-8")
                    .getReactionMetadata().get("supplier.reaction.smarts"));
        } finally {
            Files.walk(tempDir)
                    .sorted(Comparator.reverseOrder())
                    .forEach(path -> {
                        try {
                            Files.deleteIfExists(path);
                        } catch (IOException ignored) {
                        }
                    });
        }
    }
}
