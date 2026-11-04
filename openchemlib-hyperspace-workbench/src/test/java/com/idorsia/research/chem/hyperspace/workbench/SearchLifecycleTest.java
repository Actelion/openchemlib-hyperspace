package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class SearchLifecycleTest {
    @TempDir Path temp;
    private AssetService.Asset register(WorkspaceStore store,AssetService assets,Path path,AssetService.Kind kind)throws Exception {
        return assets.register(store.snapshot().revision,"scientist",assets.inspect(path,kind),path.getFileName().toString(),null,null,null,"toy-v1","Synthetic redistributable toy");
    }
    private JobService.Job terminal(JobService jobs,String id)throws Exception {
        long deadline=System.nanoTime()+Duration.ofSeconds(30).toNanos();
        while(System.nanoTime()<deadline){var j=jobs.get(id);if(Set.of("completed","failed","cancelled","interrupted").contains(j.state()))return j;Thread.sleep(10);}
        fail("Job failed to reach terminal state: "+jobs.get(id));return null;
    }
    @Test void explicitLifecycleSharesAliasesProtectsActiveJobsAndResetsOnRestart()throws Exception {
        Path raw=temp.resolve("toy.rawspace.gz");RawSynthonSpaceIO.write(ToySetup.raw(),raw);String assetId;
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var assets=new AssetService(store)) {
            var one=register(store,assets,raw,AssetService.Kind.RAW_SPACE);var two=register(store,assets,raw,AssetService.Kind.RAW_SPACE);assetId=one.id();
            assertFalse((boolean)assets.resourceStatus(one.id()).get("loaded"));
            assertThrows(WorkbenchException.class,()->assets.acquire(one.id()));
            var loaded=assets.load(store.snapshot().revision,"scientist",one.id(),null);
            assertEquals("8",loaded.metadata().path("representedCombinations").asText());
            assertEquals(2,loaded.metadata().path("reactionCount").asInt());
            assets.load(store.snapshot().revision,"agent:local",two.id(),null);
            try(var a=assets.acquire(one.id());var b=assets.acquire(two.id())) {
                assertSame(a.value(),b.value());
                var exception=assertThrows(WorkbenchException.class,()->assets.unload(store.snapshot().revision,"scientist",two.id()));
                assertEquals("resource_busy",exception.code());
            }
            assets.unload(store.snapshot().revision,"scientist",two.id());assertTrue((boolean)assets.resourceStatus(one.id()).get("loaded"));
        }
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var assets=new AssetService(store)) {
            assertEquals("unloaded",assets.get(assetId).status());assertFalse((boolean)assets.resourceStatus(assetId).get("loaded"));
        }
    }
    @Test void invalidBundleIsRejectedAndFailedLoadingHasNoResource()throws Exception {
        Path bundle=temp.resolve("bundle");Files.createDirectories(bundle);
        var manifest=WorkspaceStore.JSON.readTree(Path.of("../openchemlib-hyperspace-3d/model-bundles/deepspace7-v1/manifest.json").toFile());
        ((com.fasterxml.jackson.databind.node.ObjectNode)manifest).put("atomFeatureDim",56).put("graphFeatureContract","deepspace7-v3");
        WorkspaceStore.JSON.writeValue(bundle.resolve("manifest.json").toFile(),manifest);
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var assets=new AssetService(store)) {
            assertThrows(Exception.class,()->assets.inspect(bundle,null));
            Path bad=temp.resolve("invalid.data");Files.write(bad,new byte[]{0x1f,(byte)0x8b,0,0});
            var asset=register(store,assets,bad,AssetService.Kind.SUBSTRUCTURE_INDEX);
            assertThrows(Exception.class,()->assets.load(store.snapshot().revision,"scientist",asset.id(),null));
            assertEquals("failed",assets.get(asset.id()).status());assertNotNull(assets.get(asset.id()).error());assertFalse((boolean)assets.resourceStatus(asset.id()).get("loaded"));
        }
    }
    @Test void jobsPersistFailureAndCancellationWaitsForResourceCleanup()throws Exception {
        Path raw=temp.resolve("toy.rawspace.gz");RawSynthonSpaceIO.write(ToySetup.raw(),raw);
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var assets=new AssetService(store);var jobs=new JobService(store,1,2)) {
            var asset=register(store,assets,raw,AssetService.Kind.RAW_SPACE);assets.load(store.snapshot().revision,"scientist",asset.id(),null);
            var started=new CountDownLatch(1);var release=new CountDownLatch(1);
            var job=jobs.submit(store.snapshot().revision,"agent:local","test-held-resource",null,Map.of(),List.of(assets.acquire(asset.id())),ctx->{ctx.detail("workerHandle","retained-handle");started.countDown();release.await();ctx.check();return Map.of();});
            assertTrue(started.await(5,TimeUnit.SECONDS));
            synchronized(store){jobs.cancel(store.snapshot().revision,"scientist",job.id());}
            assertEquals("cancelling",jobs.get(job.id()).state());assertEquals("retained-handle",jobs.get(job.id()).request().path("workerHandle").asText());
            assertThrows(WorkbenchException.class,()->assets.unload(store.snapshot().revision,"scientist",asset.id()));
            release.countDown();assertEquals("cancelled",terminal(jobs,job.id()).state());assertThrows(WorkbenchException.class,()->jobs.cancel(store.snapshot().revision,"scientist",job.id()));assertEquals("cancelled",jobs.get(job.id()).state());
            assertEquals(0,assets.resourceStatus(asset.id()).get("activeJobs"));
            var failed=jobs.submit(store.snapshot().revision,"scientist","test-failure",null,Map.of(),List.of(),ctx->{throw new IllegalStateException("intentional test failure");});
            var terminal=terminal(jobs,failed.id());assertEquals("failed",terminal.state());assertTrue(terminal.error().contains("intentional test failure"));
            assets.unload(store.snapshot().revision,"scientist",asset.id());
        }
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var jobs=new JobService(store,1,2)) {
            assertEquals(2,jobs.list().size());assertTrue(jobs.list().stream().allMatch(j->Set.of("cancelled","failed").contains(j.state())));
        }
    }
    @Test void realCombinatorialSearchBoundedAssemblyCollectionsAndExportSurviveReconnect()throws Exception {
        Map<String,String> toy=ToySetup.prepare(temp.resolve("toy"),1);String expansionId,collectionId,queryId;
        try(var store=new WorkspaceStore(temp.resolve("workspace"));var assets=new AssetService(store);var jobs=new JobService(store,2,4)) {
            var queryService=new QueryService(store);var results=new ResultService(store);var search=new SubstructureService(store,assets,queryService,jobs,results);
            var asset=register(store,assets,Path.of(toy.get("substructureIndex")),AssetService.Kind.SUBSTRUCTURE_INDEX);
            assets.load(store.snapshot().revision,"scientist",asset.id(),null);
            var query=queryService.create(store.snapshot().revision,"scientist","smiles","CCO",true,null,"Alcohol fragment",null);queryId=query.query().id();
            var started=search.search(store.snapshot().revision,"scientist",queryId,List.of(asset.id()),new SubstructureService.Options(2,100,true));
            var done=terminal(jobs,started.id());assertEquals("completed",done.state(),done.error());
            var families=search.families(started.runId(),0,100);assertTrue(families.total()>=2);assertTrue(families.items().stream().allMatch(f->f.positions().size()==2));
            assertEquals(0,results.getRun(started.runId()).retainedCount());
            var expansion=search.expand(store.snapshot().revision,"agent:local",started.runId(),families.items().stream().map(SubstructureService.Family::id).toList(),3,true);
            assertEquals("completed",terminal(jobs,expansion.id()).state());expansionId=expansion.runId();
            var run=results.getRun(expansionId);assertEquals(3,run.diagnostics().path("attemptedProducts").asInt());assertTrue(run.retainedCount()>0&&run.retainedCount()<=3);
            List<ResultService.Hit> rows=results.all(expansionId);assertTrue(rows.stream().allMatch(h->!h.matchedAtoms().isEmpty()&&h.synthons().size()==2&&h.assetId().equals(asset.id())));
            assertTrue(rows.stream().allMatch(h->h.provenance().path("parentRunId").asText().equals(started.runId())));
            assertEquals(1,results.page(expansionId,0,1,"weight",false,null).items().size());
            assertTrue(results.page(expansionId,0,2,null,false,new ResultService.Filters(1e6,null,null,null)).items().isEmpty());
            assertThrows(WorkbenchException.class,()->results.page(expansionId,0,201,null,false,null));
            var collection=results.saveCollection(store.snapshot().revision,"scientist",null,"Alcohol examples","Retained source occurrences",rows.stream().map(ResultService.Hit::id).toList());collectionId=collection.id();
            var comparison=results.compare(expansionId,expansionId,10);assertEquals(run.uniqueMoleculeCount(),comparison.shared());assertEquals(0,comparison.onlyLeft());
            assertEquals(Math.min(2,run.uniqueMoleculeCount()),results.diverse(expansionId,2).size());
            var export=results.export(null,collectionId);Path output=Path.of((String)export.get("directory"));
            assertEquals(rows.size()+1,Files.readAllLines(output.resolve("results.tsv")).size());
            assertEquals(rows.size(),Files.readString(output.resolve("results.sdf")).split("\\$\\$\\$\\$",-1).length-1);
            assertEquals(1,WorkspaceStore.JSON.readTree(output.resolve("provenance.json").toFile()).path("formatVersion").asInt());
            queryService.create(store.snapshot().revision,"scientist","smiles","c1ccccc1",true,null,"New editor query",null);
            assertEquals(queryId,results.getRun(expansionId).queryId());assertEquals(query.query().document(),results.getRun(expansionId).query().document());
            assets.unload(store.snapshot().revision,"scientist",asset.id());
        }
        try(var store=new WorkspaceStore(temp.resolve("workspace"))) {
            var results=new ResultService(store);assertEquals("completed",results.getRun(expansionId).state());assertEquals(queryId,results.getRun(expansionId).queryId());
            assertFalse(results.collectionHits(collectionId).isEmpty());assertTrue(results.hit(results.collectionHits(collectionId).get(0).id()).id().startsWith(expansionId));
        }
    }
}
