package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.index.build.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import com.idorsia.research.chem.hyperspace3d.screening.*;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class LearnedWorkflowTest {
    @TempDir static Path shared;
    @TempDir Path temp;
    static Path modelPath=Path.of("../openchemlib-hyperspace-3d/model-bundles/deepspace7-v1").toAbsolutePath().normalize();
    static Path compactPath=Path.of("../openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16").toAbsolutePath().normalize();
    static DeepSpaceModelBundle model;static CompactSkelSpheresModelBundle compact;static Path cache;
    @BeforeAll static void fingerprints()throws Exception {
        model=DeepSpaceModelBundle.load(modelPath);compact=CompactSkelSpheresModelBundle.load(compactPath,model.manifest());cache=shared.resolve("cache");Path table=shared.resolve("molecules.tsv");
        Files.writeString(table,"id\tsmiles\nhexane\tCCCCCC\nhexanol\tCCCCCCO\nhexanol-copy\tOCCCCCC\nbenzene\tc1ccccc1\ncaffeine\tCn1c(=O)c2c(ncn2C)n(C)c1=O\nacid\tO=C(O)c1ccccc1\nstereo\tC[C@H](O)CCCN\n");
        var config=new MoleculeFingerprintIndexBuildConfig();config.inputs.library=table.toString();config.inputs.modelBundle=modelPath.toString();config.inputs.compactBundle=compactPath.toString();config.output.directory=cache.toString();
        config.output.sourceRowsPerShard=4;config.runtime.device="CPU";config.runtime.cpuWorkers=2;config.runtime.encoderBatchSize=4;
        var runtime=new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.CPU,0,2);
        try(var encoder=new DeepSpaceV1Encoder(runtime,model);var projector=new CompactSkelSpheresProjector(runtime,compact)) {
            new MoleculeFingerprintIndexBuilder(batch->{float[][] base=encoder.encode(batch);return new MoleculeFingerprintBatch(base,projector.project(base));},config,
                    new MoleculeFingerprintIndexBuildConfig.ResolvedPaths(table,modelPath,compactPath,cache)).withModelIdentity(model.bundleHash(),compact.bundleHash()).build();
        }
        assertEquals(7,MoleculeFingerprintIndexReader.loadManifest(cache).recordCount);
    }
    private AssetService.Asset register(WorkbenchApplication app,Path path,AssetService.Kind kind)throws Exception {
        var asset=app.assets.register(app.store.snapshot().revision,"scientist",app.assets.inspect(path,kind),path.getFileName().toString(),null,null,null,"test","Public model + synthetic molecules");return asset;
    }
    private String[] load(WorkbenchApplication app)throws Exception {
        var m=register(app,modelPath,AssetService.Kind.GRAPH_MODEL);app.assets.load(app.store.snapshot().revision,"scientist",m.id(),null);
        var c=register(app,compactPath,AssetService.Kind.COMPACT_MODEL);app.assets.load(app.store.snapshot().revision,"scientist",c.id(),m.id());
        var index=register(app,cache,AssetService.Kind.FLAT_INDEX);app.assets.load(app.store.snapshot().revision,"scientist",index.id(),null);return new String[]{m.id(),c.id(),index.id()};
    }
    private LearnedSearchService.Request request(WorkbenchApplication app,String[] ids)throws Exception {
        var query=app.queries.create(app.store.snapshot().revision,"scientist","smiles","CCCCCCO",false,null,"Hexanol",null);
        var r=new LearnedSearchService.Request();r.queryId=query.query().id();r.modelAssetId=ids[0];r.compactAssetId=ids[1];r.indexAssetId=ids[2];r.topK=3;r.shortlistTopK=6;r.batchSize=3;r.scanBatchSize=3;r.encoderBatchSize=3;return r;
    }
    private JobService.Job done(WorkbenchApplication app,JobService.Job job)throws Exception {
        long deadline=System.nanoTime()+Duration.ofSeconds(30).toNanos();while(System.nanoTime()<deadline){var j=app.jobs.get(job.id());if(Set.of("completed","failed","cancelled").contains(j.state())){assertEquals("completed",j.state(),j.error());return j;}Thread.sleep(10);}fail("Job timeout "+app.jobs.get(job.id()));return null;
    }
    @Test void nativeFlat2d3dCascadeMatchExistingRealOnnxScreenersAndReuseSessions()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace-"+UUID.randomUUID()))){String[] ids=load(app);var r=request(app,ids);assertEquals(true,app.learned.validate(r).get("compatible"));
            var graphJob=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));var hits=app.results.all(graphJob.runId());assertEquals(3,hits.size());
            assertEquals(7,app.results.getRun(graphJob.runId()).diagnostics().path("recordsScanned").asInt());assertEquals(1,app.assets.attachmentCount(ids[0]));
            var runtime=new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.CPU,0,2);var molecule=LearnedSearchService.concrete(app.queries.revision(r.queryId).document());
            float[] q;
            try(var encoder=new DeepSpaceV1Encoder(runtime,model);var compare=new DeepSpaceV1Comparator(runtime,model);var source=MoleculeFingerprintDataSources.open(cache)) {
                q=encoder.encode(new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(molecule)))[0];
                var expected=new MoleculePheSAScreener(compare,model.manifest(),ScreeningObjective.direct("phesa_total")).screen(source,q,3,3,2);
                for(int i=0;i<hits.size();i++){assertEquals(expected.hits().get(i).molecule().moleculeId(),hits.get(i).provenance().path("sourceMoleculeId").asText());assertEquals(expected.hits().get(i).objectiveScore(),hits.get(i).scores().get("objective_score"),1e-7);}
            }
            r.batchSize=7;var repeat=done(app,app.learned.search(app.store.snapshot().revision,"agent:mcp",r));assertEquals(1,app.assets.attachmentCount(ids[0]));
            assertEquals(hits.stream().map(h->h.provenance().path("sourceMoleculeId").asText()).toList(),app.results.all(repeat.runId()).stream().map(h->h.provenance().path("sourceMoleculeId").asText()).toList());
            r.mode=LearnedSearchService.Mode.SKELSPHERES_2D;var two=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));var exact=app.results.all(two.runId());
            assertEquals(1.0,exact.get(0).scores().get("exact_skelspheres"),1e-7);assertEquals(exact.get(0).moleculeId(),exact.get(1).moleculeId());assertNotEquals(exact.get(0).id(),exact.get(1).id());
            assertTrue(app.results.getRun(two.runId()).scope().contains("shortlist"));assertEquals(1,app.assets.attachmentCount(ids[1]));
            r.mode=LearnedSearchService.Mode.SKELSPHERES_THEN_PHESA;var cascade=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));
            assertEquals(7,app.results.getRun(cascade.runId()).diagnostics().path("recordsScanned").asInt());assertEquals(6,app.results.getRun(cascade.runId()).diagnostics().path("comparatorCandidates").asInt());
            assertTrue(app.results.all(cascade.runId()).stream().allMatch(h->h.provenance().has("compactRank")));
            LearnedSearchService.GraphSessions nativeSessions;
            try(var lease=app.assets.acquire(ids[0])) {nativeSessions=app.assets.attachment(ids[0],"graph:CPU:0:2",LearnedSearchService.GraphSessions.class,()->{throw new AssertionError("Expected reused sessions");});}
            var alias=register(app,modelPath,AssetService.Kind.GRAPH_MODEL);app.assets.load(app.store.snapshot().revision,"scientist",alias.id(),null);assertEquals(1,app.assets.attachmentCount(alias.id()));
            app.assets.unload(app.store.snapshot().revision,"scientist",ids[0]);assertEquals(0,app.assets.attachmentCount(ids[0]));
            assertEquals(1,nativeSessions.encoder.encode(new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(molecule))).length);
            app.assets.unload(app.store.snapshot().revision,"scientist",alias.id());
            assertThrows(RuntimeException.class,()->nativeSessions.encoder.encode(new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(molecule))));
        }
    }
    @Test void boundedRankingRetainsOccurrencesAndRecordsRejectedGraphs()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace-"+UUID.randomUUID()))){String[] ids=load(app);var r=request(app,ids);var query=app.queries.revision(r.queryId);
            var source=app.results.create(app.store.snapshot().revision,"scientist","bounded_expansion",query,List.of(ids[2]),Map.of("expansionAttemptCap",4),"Test retained collection");
            List<ResultService.Hit> hits=new ArrayList<>();for(int n=0;n<3;n++){var m=new StereoMolecule();new SmilesParser().parse(m,n==2?"CCO":"CCCCCCO");hits.add(ResultService.makeHit(source.id(),"source-"+n,ids[2],"family-public","reaction-public",List.of(),m,null,false,Map.of(),Map.of("sourceTuple",n)));}
            app.results.complete(source.id(),hits,Map.of("attemptedProducts",3));r.indexAssetId=null;
            var ranked=done(app,app.learned.rank(app.store.snapshot().revision,"agent:mcp",r,source.id(),null,10));var actual=app.results.all(ranked.runId());
            assertEquals(2,actual.size());assertEquals(1,actual.stream().map(ResultService.Hit::moleculeId).distinct().count());assertEquals(2,actual.stream().map(ResultService.Hit::id).distinct().count());
            var diag=app.results.getRun(ranked.runId()).diagnostics();assertEquals(3,diag.path("attemptedCandidates").asInt());assertEquals(1,diag.path("acceptedCandidateGraphs").asInt());assertEquals(1,diag.path("rejectedCandidateOccurrences").asInt());
            assertEquals(1,diag.path("rejectionReasons").path("TOO_FEW_HEAVY_ATOMS").asInt());assertTrue(actual.stream().allMatch(h->h.reaction().equals("reaction-public")&&h.provenance().has("sourceOccurrenceId")));
            var toy=ToySetup.prepare(temp.resolve("assembly-to-ranking"),1);var index=register(app,Path.of(toy.get("substructureIndex")),AssetService.Kind.SUBSTRUCTURE_INDEX);app.assets.load(app.store.snapshot().revision,"scientist",index.id(),null);
            var fragment=app.queries.create(app.store.snapshot().revision,"scientist","smiles","CCO",true,null,"Assembly fragment",null).query();
            var search=done(app,app.substructure.search(app.store.snapshot().revision,"scientist",fragment.id(),List.of(index.id()),new SubstructureService.Options(2,20,true)));
            var families=app.substructure.families(search.runId(),0,20);var expansion=done(app,app.substructure.expand(app.store.snapshot().revision,"scientist",search.runId(),families.items().stream().map(SubstructureService.Family::id).toList(),4,true));
            var realRank=done(app,app.learned.rank(app.store.snapshot().revision,"scientist",r,expansion.runId(),null,4));var products=app.results.all(realRank.runId());
            assertFalse(products.isEmpty());assertTrue(products.stream().allMatch(h->h.synthons().size()==2&&h.provenance().path("sourceRunId").asText().equals(expansion.runId())));
            assertTrue(app.results.getRun(realRank.runId()).settings().path("sourceRuns").path(expansion.runId()).has("expansionSettings"));
        }
    }
    @Test void rejectsQueryFeaturesAndIncompatibleCacheBeforeAdmittingAJob()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace-"+UUID.randomUUID()))){String[] ids=load(app);var r=request(app,ids);var original=app.queries.inspect(r.queryId);
            var changed=app.queries.edit(app.store.snapshot().revision,"agent:mcp",r.queryId,List.of(new QueryChemistry.Edit(QueryChemistry.Operation.ATOM_ELEMENTS,List.of(original.atoms().get(0).ref()),List.of(6,7),false,null,null,null,null,null)),"Query alternatives");
            r.queryId=changed.query().id();assertThrows(WorkbenchException.class,()->app.learned.search(app.store.snapshot().revision,"scientist",r));assertTrue(app.jobs.list().isEmpty(),app.jobs.list().toString());
            r.queryId=app.queries.create(app.store.snapshot().revision,"scientist","smiles","CCCCCCO",false,null,"Concrete",null).query().id();
            Path bad=temp.resolve("bad-cache-"+UUID.randomUUID());Files.createDirectories(bad);
            try(var paths=Files.walk(cache)){for(Path file:paths.toList()){Path target=bad.resolve(cache.relativize(file));if(Files.isDirectory(file))Files.createDirectories(target);else Files.copy(file,target);}}
            var manifest=WorkspaceStore.JSON.readTree(cache.resolve("manifest.json").toFile());((com.fasterxml.jackson.databind.node.ObjectNode)manifest).put("modelBundleHash","wrong-content");WorkspaceStore.JSON.writeValue(bad.resolve("manifest.json").toFile(),manifest);
            var asset=register(app,bad,AssetService.Kind.FLAT_INDEX);app.assets.load(app.store.snapshot().revision,"scientist",asset.id(),null);r.indexAssetId=asset.id();
            assertThrows(IllegalArgumentException.class,()->app.learned.search(app.store.snapshot().revision,"scientist",r));assertTrue(app.jobs.list().isEmpty(),app.jobs.list().toString());assertEquals(0,app.assets.resourceStatus(asset.id()).get("activeJobs"));
        }
    }
    @Test void cpuDistributionRejectsCudaWithoutFallbackAndReleasesLeases()throws Exception {
        Assumptions.assumeFalse(ai.onnxruntime.OrtEnvironment.getAvailableProviders().contains(ai.onnxruntime.OrtProvider.CUDA),"CPU-distribution negative test");
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){String[] ids=load(app);var r=request(app,ids);r.device="CUDA";
            var submitted=app.learned.search(app.store.snapshot().revision,"scientist",r);long deadline=System.nanoTime()+Duration.ofSeconds(20).toNanos();
            JobService.Job status;do{status=app.jobs.get(submitted.id());if(status.state().equals("failed"))break;Thread.sleep(10);}while(System.nanoTime()<deadline);
            assertEquals("failed",status.state());assertTrue(status.error().contains("CPU fallback is disabled"));assertEquals("CUDA",app.results.getRun(submitted.runId()).settings().path("device").asText());
            assertEquals(0,app.assets.resourceStatus(ids[0]).get("activeJobs"));assertEquals(0,app.assets.attachmentCount(ids[0]));
        }
    }
    @Test void pythonSupplierFormatMatchesNativeRealLearnedScores()throws Exception {
        Path python=PythonCacheFixture.copy(cache,temp.resolve("python"),model,compact);
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){String[] ids=load(app);var r=request(app,ids);
            var source=register(app,python,AssetService.Kind.PYTHON_CACHE);app.assets.load(app.store.snapshot().revision,"scientist",source.id(),null);
            for(var mode:List.of(LearnedSearchService.Mode.PHESA_3D,LearnedSearchService.Mode.SKELSPHERES_2D,LearnedSearchService.Mode.SKELSPHERES_THEN_PHESA)) {
                r.mode=mode;r.indexAssetId=ids[2];var nativeJob=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));var expected=app.results.all(nativeJob.runId());
                r.indexAssetId=source.id();var pythonJob=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));var actual=app.results.all(pythonJob.runId());
                assertEquals(expected.stream().map(ResultService.Hit::moleculeId).toList(),actual.stream().map(ResultService.Hit::moleculeId).toList());
                for(int n=0;n<actual.size();n++)assertEquals(expected.get(n).scores().get("objective_score"),actual.get(n).scores().get("objective_score"),1e-7);
            }
        }
    }
    @Test void suppliedPoseUsesExternalVerifiedBundleAndRecordsDirectedSemantics()throws Exception {
        Path posePath=Path.of(System.getProperty("hyperspace.poseBundle","../../pose_deployment/shared-candidate-seed17-verified")).toAbsolutePath().normalize();
        Assumptions.assumeTrue(Files.isRegularFile(posePath.resolve("manifest.json")),"External pose bundle not configured");
        try(var app=new WorkbenchApplication(temp.resolve("workspace-"+UUID.randomUUID()))){String[] ids=load(app);var r=request(app,ids);var p=register(app,posePath,AssetService.Kind.POSE_MODEL);app.assets.load(app.store.snapshot().revision,"scientist",p.id(),ids[0]);
            var args=WorkspaceStore.JSON.valueToTree(Map.of("path",Path.of("../openchemlib-hyperspace-3d/src/test/resources/pose-golden/query-0.mol").toAbsolutePath().toString(),"name","Supplied pose"));
            var query=(QueryService.PublishResult)app.invoke("query.import_pose",args,app.store.snapshot().revision,"scientist");r.queryId=query.query().id();r.mode=LearnedSearchService.Mode.POSE_PHESA_3D;r.poseAssetId=p.id();r.posePpWeight=.25;
            var job=done(app,app.learned.search(app.store.snapshot().revision,"scientist",r));var run=app.results.getRun(job.runId());assertEquals("learned_pose_to_ensemble",run.diagnostics().path("scoreSemantics").asText());
            assertEquals(.25,run.diagnostics().path("posePpWeight").asDouble());assertEquals(3,app.results.all(run.id()).size());
            for(var h:app.results.all(run.id()))assertEquals(.25*h.scores().get("phesa_pharmacophore")+.75*h.scores().get("phesa_shape"),h.scores().get("phesa_total"),1e-6);
            assertEquals(1,app.assets.attachmentCount(p.id()));
        }
    }
}
