package com.idorsia.research.chem.hyperspace.workbench;

import com.idorsia.research.chem.hyperspace.rawspace.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class CatalogueFormatsTest {
    @TempDir Path temp;
    @Test void downsampledHeaderAndStandaloneProductManifestHaveExplicitCoverageWithoutLoadingPayloads()throws Exception {
        var builder=RawSynthonSpace.builder("Reduced public toy").putMetadata(RawSynthonSpace.MetadataKeys.SPACE_ROLE,RawSynthonSpace.MetadataKeys.SPACE_ROLE_DOWNSAMPLED);
        for(var reaction:ToySetup.raw().getReactions().entrySet())for(var set:reaction.getValue().getRawFragmentSets().entrySet())builder.addRawFragments(reaction.getKey(),set.getKey(),set.getValue().subList(0,1));
        Path reduced=temp.resolve("reduced.rawspace.gz");RawSynthonSpaceIO.write(builder.build(),reduced);
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){var inspection=app.assets.inspect(reduced,null);assertEquals(AssetService.Kind.DOWNSAMPLED_SPACE,inspection.kind());assertTrue(inspection.coverage().contains("Reduced"));assertNull(inspection.recordCount());
            assertThrows(WorkbenchException.class,()->app.assets.inspect(reduced,AssetService.Kind.RAW_SPACE));
            var asset=app.assets.register(0,"scientist",inspection,"Reduced",null,null,null,"v1","Reduced source");assertFalse((boolean)app.assets.resourceStatus(asset.id()).get("loaded"));
            var manifest=new ProductFingerprintIndexManifest();manifest.rawspaceIdentity="full";manifest.rawspaceSha256="a".repeat(64);manifest.downsampledRawspaceIdentity="reduced";manifest.downsampledRawspaceSha256="b".repeat(64);manifest.modelBundleHash="model";manifest.samplingAlgorithm="seeded";manifest.reactionWeightConfiguration=Map.of();manifest.molecularFilters=Map.of();manifest.recordCount=2;manifest.shards=List.of(new ProductFingerprintIndexShard("not-read-by-registration.fp",2,"c".repeat(64)));
            Path file=temp.resolve("fingerprints.json");WorkspaceStore.JSON.writeValue(file.toFile(),manifest);var product=app.assets.inspect(file,null);
            assertEquals(AssetService.Kind.PRODUCT_INDEX,product.kind());assertEquals(2,product.recordCount());assertTrue(product.coverage().contains("not exhaustive"));assertFalse(Files.exists(temp.resolve("not-read-by-registration.fp")));
        }
    }
    @Test void canonicalManifestFileNormalizesToBundleDirectoryAndLegacyIndexesAreRecognized()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){Path model=Path.of("../openchemlib-hyperspace-3d/model-bundles/deepspace7-v1").toRealPath();assertEquals(model.toString(),app.assets.inspect(model.resolve("manifest.json"),null).path());
            Path similarity=temp.resolve("public_FragFp_similarity3.data");try(var output=new GZIPOutputStream(Files.newOutputStream(similarity))){output.write(new byte[]{1,2,3});}
            var detection=app.assets.inspect(similarity,null);assertEquals(AssetService.Kind.SIMILARITY_INDEX,detection.kind());assertTrue(detection.capabilities().contains("legacy_similarity"));
            var asset=app.assets.register(0,"scientist",detection,"Legacy",null,null,null,"v1","Header-only test");assertEquals("registered",asset.status());assertFalse((boolean)app.assets.resourceStatus(asset.id()).get("loaded"));
            assertThrows(Exception.class,()->app.assets.load(1,"scientist",asset.id(),null));assertEquals("failed",app.assets.get(asset.id()).status());
        }
    }
}
