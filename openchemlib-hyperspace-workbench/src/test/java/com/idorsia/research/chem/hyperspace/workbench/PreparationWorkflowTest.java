package com.idorsia.research.chem.hyperspace.workbench;

import java.nio.file.*;
import java.time.Duration;
import java.util.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class PreparationWorkflowTest {
    @TempDir Path temp;
    private void configure(WorkbenchApplication app)throws Exception {
        Path root=Path.of("..").toAbsolutePath().normalize();
        app.preparation.configure(app.store.snapshot().revision,"scientist",Map.of("formatVersion",1,"cliJar",System.getProperty("hyperspace.cliJar",root.resolve("openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar").toString()),
                "fingerprintJar",root.resolve("openchemlib-hyperspace-3d/target/openchemlib-hyperspace-3d-3.1.0.jar").toString(),"fingerprintLibDirectory",root.resolve("openchemlib-hyperspace-3d/target/lib").toString(),
                "modelBundle",root.resolve("openchemlib-hyperspace-3d/model-bundles/deepspace7-v1").toString(),"compactBundle",root.resolve("openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16").toString(),"threads",2,"heap","512M"));
    }
    private JobService.Job done(WorkbenchApplication app,JobService.Job job)throws Exception {
        long deadline=System.nanoTime()+Duration.ofSeconds(60).toNanos();while(System.nanoTime()<deadline){var current=app.jobs.get(job.id());if(Set.of("completed","cancelled","failed","interrupted").contains(current.state()))return current;Thread.sleep(25);}fail("Preparation timeout: "+app.jobs.get(job.id()));return null;
    }
    @Test void realSupplierImportAndBuildRegisterAssetsWithoutLoadingThemAndReattachIdempotently()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){configure(app);Path input=temp.resolve("synthons.tsv");Files.writeString(input,"SMILES\tfragmentId\tsynthonSetIndex\treactionId\n[U]CCC\ta\t0\tpublic\n[U]CCN\tb\t0\tpublic\n[U]CCO\tx\t1\tpublic\n[U]CCC\ty\t1\tpublic\n");
            var sample=WorkspaceStore.JSON.valueToTree(app.preparation.inspect(input.toString(),null,false));assertTrue(sample.path("sampleOnly").asBoolean());
            var options=Map.<String,Object>of("input",input.toString(),"format","enamine","spaceName","Public source recipe","searchModes",List.of("substructure","similarity"),"threads",2,"heap","512M");
            var recipe=app.preparation.recipe("space",options);assertEquals(1,recipe.get("formatVersion"));assertTrue(recipe.get("resume").toString().contains("rawspace"));
            var submitted=app.preparation.submit(app.store.snapshot().revision,"scientist","space",options);var finished=done(app,submitted);if(finished.state().equals("failed"))System.err.println(app.preparation.log(submitted.id(),"worker",0,8192));assertEquals("completed",finished.state(),finished.error());
            assertEquals(3,app.store.snapshot().assets.size());assertTrue(app.store.snapshot().assets.values().stream().allMatch(a->a.path("status").asText().equals("registered")));
            var index=app.store.snapshot().assets.values().stream().map(n->WorkspaceStore.JSON.convertValue(n,AssetService.Asset.class)).filter(a->a.kind()==AssetService.Kind.SUBSTRUCTURE_INDEX).findFirst().orElseThrow();
            assertFalse((boolean)app.assets.resourceStatus(index.id()).get("loaded"));
            var log=WorkspaceStore.JSON.valueToTree(app.preparation.log(submitted.id(),"import",0,1024));assertTrue(log.toString().length()>30);
            var reattached=done(app,app.preparation.attach(app.store.snapshot().revision,"agent:mcp",submitted.id()));assertEquals("completed",reattached.state(),reattached.error());assertEquals(3,app.store.snapshot().assets.size());
            var similarity=app.store.snapshot().assets.values().stream().map(n->WorkspaceStore.JSON.convertValue(n,AssetService.Asset.class)).filter(a->a.kind()==AssetService.Kind.SIMILARITY_INDEX).findFirst().orElseThrow();
            app.assets.load(app.store.snapshot().revision,"scientist",similarity.id(),null);assertEquals("ready",app.assets.get(similarity.id()).status());
            app.assets.load(app.store.snapshot().revision,"scientist",index.id(),null);assertEquals("4",app.assets.get(index.id()).metadata().path("representedCombinations").asText());
            var query=app.queries.create(app.store.snapshot().revision,"scientist","smiles","CCO",true,null,"Query imported index",null);
            var search=done(app,app.substructure.search(app.store.snapshot().revision,"scientist",query.query().id(),List.of(index.id()),new SubstructureService.Options(2,20,true)));
            assertEquals("completed",search.state(),search.error());assertTrue(app.substructure.families(search.runId(),0,20).total()>0);
        }
    }
    @Test void preparingDownsampledRawspaceKeepsReducedCoverageInRegisteredOutputs()throws Exception {
        var builder=com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace.builder("Reduced public source").putMetadata(com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace.MetadataKeys.SPACE_ROLE,"downsampled");
        for(var reaction:ToySetup.raw().getReactions().entrySet())for(var set:reaction.getValue().getRawFragmentSets().entrySet())builder.addRawFragments(reaction.getKey(),set.getKey(),set.getValue().subList(0,1));
        Path raw=temp.resolve("reduced.rawspace.gz");com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO.write(builder.build(),raw);
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){configure(app);
            var options=Map.<String,Object>of("input",raw.toString(),"format","rawspace","spaceName","Reduced source","searchModes",List.of("substructure"),"threads",2,"heap","512M");
            var submitted=app.preparation.submit(app.store.snapshot().revision,"scientist","space",options);var finished=done(app,submitted);assertEquals("completed",finished.state(),finished.error());
            var assets=app.store.snapshot().assets.values().stream().map(n->WorkspaceStore.JSON.convertValue(n,AssetService.Asset.class)).toList();assertEquals(2,assets.size());
            assertTrue(assets.stream().anyMatch(a->a.kind()==AssetService.Kind.DOWNSAMPLED_SPACE));assertFalse(assets.stream().anyMatch(a->a.kind()==AssetService.Kind.RAW_SPACE));
            var repeated=done(app,app.preparation.attach(app.store.snapshot().revision,"scientist",submitted.id()));assertEquals("completed",repeated.state(),repeated.error());assertEquals(2,app.store.snapshot().assets.size());
        }
    }
    @Test void fingerprintPreparationIsDetachedResumableAndRegistersModelBoundCache()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){configure(app);Path input=temp.resolve("molecules.tsv");Files.writeString(input,"id\tsmiles\nhexanol\tCCCCCCO\nbenzene\tc1ccccc1\nsmall\tCCO\n");Path output=temp.resolve("fingerprints");
            var sample=WorkspaceStore.JSON.valueToTree(app.preparation.inspect(input.toString(),null,true));assertTrue(sample.path("bounded").asBoolean());
            var options=Map.<String,Object>of("input",input.toString(),"outputDirectory",output.toString(),"device","CPU","threads",2,"heap","512M","encoderBatchSize",2,"sourceRowsPerShard",2);
            var recipe=app.preparation.recipe("fingerprints",options);assertTrue(recipe.get("resume").toString().contains("committed"));
            var submitted=app.preparation.submit(app.store.snapshot().revision,"agent:mcp","fingerprints",options);var finished=done(app,submitted);if(finished.state().equals("failed"))System.err.println(app.preparation.log(submitted.id(),"worker",0,8192));assertEquals("completed",finished.state(),finished.error());
            var manifest=WorkspaceStore.JSON.readTree(output.resolve("manifest.json").toFile());assertEquals(2,manifest.path("recordCount").asInt());assertEquals(1,manifest.path("rejectedCount").asInt());assertFalse(manifest.path("modelBundleHash").asText().isBlank());
            assertEquals(1,app.store.snapshot().assets.size());var asset=WorkspaceStore.JSON.convertValue(app.store.snapshot().assets.values().iterator().next(),AssetService.Asset.class);assertEquals(AssetService.Kind.FLAT_INDEX,asset.kind());assertEquals(2,asset.recordCount());
            String before=Files.readString(output.resolve("manifest.json"));var repeat=done(app,app.preparation.submit(app.store.snapshot().revision,"scientist","fingerprints",options));assertEquals("completed",repeat.state(),repeat.error());assertEquals(before,Files.readString(output.resolve("manifest.json")));
        }
    }
    @Test void slurmExportCreatesFilesButNoOwnerComputationOrClusterSubmission()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))) {
            Path output=temp.resolve("bundle");var result=WorkspaceStore.JSON.valueToTree(app.preparation.exportSlurm(0,"scientist",Path.of("../fingerprint-slurm-profile.example.json").toAbsolutePath().toString(),output.toString()));
            assertTrue(result.path("exportOnly").asBoolean());assertFalse(result.path("submitted").asBoolean());assertFalse(result.path("monitored").asBoolean());
            assertTrue(Files.isRegularFile(output.resolve("workflow-config.json")));assertTrue(Files.isRegularFile(output.resolve("compute.sbatch")));assertTrue(app.jobs.list().isEmpty());
            assertEquals("preparation.slurm_export",app.store.snapshot().activity.get(0).operation());
        }
    }
    @Test void cancelledWorkerHandleIsObservedUntilTerminalAndCompletedShardsResume()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){configure(app);Path table=temp.resolve("large.tsv");var text=new StringBuilder("id\tsmiles\n");for(int n=0;n<2000;n++)text.append("public-").append(n).append("\tCCCCCCO\n");Files.writeString(table,text);Path output=temp.resolve("resume");
            var options=Map.<String,Object>of("input",table.toString(),"outputDirectory",output.toString(),"device","CPU","threads",1,"heap","512M","encoderBatchSize",1,"sourceRowsPerShard",10);
            var submitted=app.preparation.submit(app.store.snapshot().revision,"scientist","fingerprints",options);long deadline=System.nanoTime()+Duration.ofSeconds(30).toNanos();
            while(!Files.isRegularFile(output.resolve("shard-00000/.complete"))&&System.nanoTime()<deadline){var status=app.jobs.get(submitted.id());assertFalse(Set.of("failed","cancelled").contains(status.state()),status.error());Thread.sleep(10);}
            assertTrue(Files.isRegularFile(output.resolve("shard-00000/.complete")));String shard=Files.readString(output.resolve("shard-00000/shard.json"));
            app.jobs.cancel(app.store.snapshot().revision,"scientist",submitted.id());var cancelled=done(app,submitted);assertEquals("cancelled",cancelled.state(),cancelled.error());assertEquals(0,app.store.snapshot().assets.size());
            assertTrue(cancelled.request().hasNonNull("detachedStatus"));assertTrue(Set.of("cancelled","succeeded").contains(cancelled.request().path("detachedStatus").path("state").asText()));
            var resumed=done(app,app.preparation.resume(app.store.snapshot().revision,"scientist",submitted.id()));assertEquals("completed",resumed.state(),resumed.error());
            assertEquals(shard,Files.readString(output.resolve("shard-00000/shard.json")));assertEquals(2000,WorkspaceStore.JSON.readTree(output.resolve("manifest.json").toFile()).path("recordCount").asInt());
        }
    }
}
