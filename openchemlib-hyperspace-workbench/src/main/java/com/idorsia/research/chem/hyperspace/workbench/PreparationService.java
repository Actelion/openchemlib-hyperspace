package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import com.idorsia.research.chem.hyperspace.mcp.DetachedPreparation;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.CancellationException;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Reuses existing detached workers; only this owner registers output assets and mutates the visible catalogue. */
public final class PreparationService {
    private static final Set<String> TERMINAL=Set.of("succeeded","failed","cancelled","interrupted");
    private final WorkspaceStore store;private final AssetService assets;private final JobService jobs;
    public PreparationService(WorkspaceStore store,AssetService assets,JobService jobs){this.store=store;this.assets=assets;this.jobs=jobs;}
    public Map<String,Object> environment() {
        var configured=store.snapshot().preparationRuntime;if(configured!=null&&!configured.isNull())return map(configured);
        Path root=Path.of("").toAbsolutePath();var defaults=new LinkedHashMap<String,Object>();defaults.put("formatVersion",1);
        defaults.put("cliJar",root.resolve("openchemlib-hyperspace-cli/target/openchemlib-hyperspace-cli.jar").toString());
        defaults.put("fingerprintJar",root.resolve("openchemlib-hyperspace-3d/target/openchemlib-hyperspace-3d-3.1.0.jar").toString());
        defaults.put("fingerprintLibDirectory",root.resolve("openchemlib-hyperspace-3d/target/lib").toString());
        defaults.put("modelBundle",root.resolve("openchemlib-hyperspace-3d/model-bundles/deepspace7-v1").toString());
        defaults.put("compactBundle",root.resolve("openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16").toString());
        defaults.put("javaExecutable",Path.of(System.getProperty("java.home"),"bin","java").toString());defaults.put("threads",2);defaults.put("heap","4G");return defaults;
    }
    public Object configure(long revision,String actor,Map<String,Object> settings)throws Exception {
        require(settings.keySet().stream().allMatch(Set.of("formatVersion","cliJar","fingerprintJar","fingerprintLibDirectory","modelBundle","compactBundle","javaExecutable","threads","heap")::contains),"Unknown preparation runtime setting");
        var runtime=new LinkedHashMap<>(environment());runtime.putAll(settings);require(((Number)runtime.get("formatVersion")).intValue()==1,"Unsupported preparation configuration version");
        driver(runtime);return store.mutate(revision,actor,"preparation.configure",s->{s.preparationRuntime=WorkspaceStore.JSON.valueToTree(runtime);return runtime;});
    }
    public Object inspect(String input,String zipEntry,boolean molecules)throws Exception {return DetachedPreparation.inspect(input,zipEntry,molecules);}
    public Map<String,Object> recipe(String kind,Map<String,Object> options)throws Exception {
        return recipe(kind,options,environment());
    }
    private Map<String,Object> recipe(String kind,Map<String,Object> options,Map<String,Object> runtime)throws Exception {
        var normalized=driver(runtime).plan(kind,options);
        if(kind.equals("space"))require(Files.isRegularFile(Path.of(runtime.get("cliJar").toString())),"CLI JAR missing; build openchemlib-hyperspace-cli and configure its absolute path");
        else {
            require(Files.isRegularFile(Path.of(runtime.get("fingerprintJar").toString())),"Fingerprint JAR missing; package the 3D thin distribution first");
            var graph=DeepSpaceModelBundle.load(Path.of(runtime.get("modelBundle").toString()));CompactSkelSpheresModelBundle.load(Path.of(runtime.get("compactBundle").toString()),graph.manifest());
        }
        var out=new LinkedHashMap<String,Object>();out.put("formatVersion",1);out.put("kind",kind);out.put("options",normalized);out.put("runtime",runtime);
        out.put("produces",kind.equals("space")?List.of("rawspace","substructure index","optional legacy similarity index"):List.of("native dual-vector flat cache","accepted/rejected counts","model/checkpoint-bound manifest"));
        out.put("coverage",kind.equals("space")?"Combinatorial combinations defined by imported synthon sets; not enumerated product fingerprints.":"Accepted complete molecules in the input table, excluding recorded featurization/filter rejections.");
        out.put("requiredDependencies",kind.equals("space")?"Existing Java CLI import/build commands; Linux setsid; explicit child heap and CPU threads.":"Thin 3D JAR/lib, canonical graph and compatible compact model, ONNX Runtime CPU or explicit CUDA provider/native libraries, Linux setsid.");
        out.put("resume",kind.equals("space")?"Completed import/rawspace artifacts remain. Rebuild from a completed rawspace to reuse import; failed index stages restart in a new job.":"Reuse the same outputDirectory to resume committed source-order fingerprint shards; model and build identities are validated by the existing builder.");
        out.put("execution","Detached local worker; browser disconnect does not stop it. Outputs are registered by the workspace owner after validation.");return out;
    }
    public JobService.Job submit(long revision,String actor,String kind,Map<String,Object> options)throws Exception {
        return start(revision,actor,recipe(kind,options),null);
    }
    private JobService.Job start(long revision,String actor,Map<String,Object> recipe,String previousJob)throws Exception {
        String kind=recipe.get("kind").toString();var runtime=cast(recipe.get("runtime"));var normalized=cast(recipe.get("options"));
        var request=new LinkedHashMap<String,Object>();request.put("kind",kind);request.put("options",normalized);request.put("runtime",runtime);request.put("recipe",recipe);
        if(previousJob!=null)request.put("resumedFromJobId",previousJob);
        return jobs.submit(revision,actor,"prepare_"+kind,null,request,List.of(),ctx->{
            var driver=driver(runtime);ctx.stage("Starting detached "+kind+" worker");var accepted=driver.submit(kind,normalized);
            require(Boolean.TRUE.equals(accepted.get("accepted")),"Detached preparation workspace is busy: "+accepted.get("reason"));String child=accepted.get("jobId").toString();
            ctx.detail("detachedJobId",child);ctx.log("Detached worker "+child+"; resource limits and recipe are retained in the job request.");return monitor(ctx,actor,driver,child,kind,normalized);
        });
    }
    public JobService.Job resume(long revision,String actor,String ownerJob)throws Exception {
        var previous=jobs.get(ownerJob);require(!Set.of("queued","running","cancelling").contains(previous.state()),"Owner preparation is still active; monitor or cancel it first");
        String kind=previous.request().path("kind").asText();require(Set.of("space","fingerprints").contains(kind),"Not a resumable preparation job");
        var runtime=map(previous.request().path("runtime"));var options=map(previous.request().path("options"));
        if(previous.request().hasNonNull("detachedJobId")) {
            var status=driver(runtime).status(previous.request().path("detachedJobId").asText());
            if(Set.of("pending","running","succeeded").contains(status.get("state")))return attach(revision,actor,ownerJob);
            var artifacts=cast(status.getOrDefault("artifacts",Map.of()));
            if(kind.equals("space")&&artifacts.containsKey("rawspace")){options.put("input",artifacts.get("rawspace"));options.put("format","rawspace");options.put("importOptions",Map.of());}
        }
        return start(revision,actor,recipe(kind,options,runtime),ownerJob);
    }
    public JobService.Job attach(long revision,String actor,String previousJobId)throws Exception {
        var previous=jobs.get(previousJobId);require(previous.request().hasNonNull("detachedJobId"),"This owner job has no recorded detached worker handle");
        String child=previous.request().path("detachedJobId").asText(),kind=previous.request().path("kind").asText();var runtime=map(previous.request().path("runtime"));var options=map(previous.request().path("options"));
        require(jobs.list().stream().noneMatch(j->Set.of("queued","running","cancelling").contains(j.state())&&j.request().path("detachedJobId").asText().equals(child)),"This detached worker is already monitored by an active owner job");
        var driver=driver(runtime);var status=driver.status(child); // Revalidate the specific process/status; never assume a state file means work is live.
        require(Set.of("pending","running","succeeded").contains(status.get("state")),"Worker is terminal without successful outputs; submit a new recipe (reuse fingerprint outputDirectory or completed rawspace)");
        return jobs.submit(revision,actor,"reattach_preparation",null,Map.of("detachedJobId",child,"previousJobId",previousJobId,"kind",kind,"runtime",runtime,"options",options),List.of(),ctx->monitor(ctx,actor,driver,child,kind,options));
    }
    public Object log(String ownerJob,String stream,long offset,int limit)throws Exception {
        require(offset>=0&&limit>=1&&limit<=65536,"Log offset must be nonnegative; limit 1–65536 bytes");var job=jobs.get(ownerJob);
        require(job.request().hasNonNull("detachedJobId"),"No detached preparation log for this job");
        return driver(map(job.request().path("runtime"))).log(job.request().path("detachedJobId").asText(),stream,offset,limit);
    }
    private Object monitor(JobService.Context ctx,String actor,DetachedPreparation driver,String child,String kind,Map<String,Object> options)throws Exception {
        String previous="";boolean cancelSent=false;
        while(true){var status=driver.status(child);String state=status.get("state").toString();
            if(ctx.cancelled()&&!cancelSent&&!TERMINAL.contains(state)){driver.cancel(child);cancelSent=true;}
            String stage=Objects.toString(status.getOrDefault("stage",state));if(!ctx.cancelled()&&!stage.equals(previous)){ctx.stage(stage);previous=stage;}
            if(TERMINAL.contains(state)) {
                ctx.detail("detachedStatus",status);
                if(ctx.cancelled()||state.equals("cancelled"))throw new CancellationException("Detached worker cancelled; completed artifacts remain");
                require(state.equals("succeeded"),"Detached worker "+state+": "+status.getOrDefault("error","Inspect stage logs"));
                var registered=register(actor,child,kind,options,cast(status.getOrDefault("artifacts",Map.of())));
                return Map.of("detachedJobId",child,"status",status,"registeredAssetIds",registered,"scope","Successful output assets only; source coverage is recorded separately.");
            }
            Thread.sleep(250);
        }
    }
    private List<String> register(String actor,String child,String kind,Map<String,Object> options,Map<String,Object> artifacts)throws Exception {
        List<String> ids=new ArrayList<>();String library=null,snapshot=null;
        for(var entry:List.of(Map.entry("rawspace",AssetService.Kind.RAW_SPACE),Map.entry("substructure",AssetService.Kind.SUBSTRUCTURE_INDEX),Map.entry("similarity",AssetService.Kind.SIMILARITY_INDEX),Map.entry("moleculeIndex",AssetService.Kind.FLAT_INDEX))) {
            if(!artifacts.containsKey(entry.getKey()))continue;
            var detected=assets.inspect(Path.of(artifacts.get(entry.getKey()).toString()),entry.getValue()==AssetService.Kind.RAW_SPACE?null:entry.getValue());
            var existing=store.snapshot().assets.values().stream().map(n->WorkspaceStore.JSON.convertValue(n,AssetService.Asset.class)).filter(a->a.metadata().path("preparedByDetachedJobId").asText().equals(child)&&a.kind()==detected.kind()).findFirst();
            if(existing.isPresent()){var asset=existing.get();ids.add(asset.id());library=asset.libraryId();snapshot=asset.sourceSnapshotId();continue;}
            var metadata=detected.metadata().deepCopy();
            ((com.fasterxml.jackson.databind.node.ObjectNode)metadata).put("preparedByDetachedJobId",child);
            var detection=new AssetService.Detection(detected.path(),detected.kind(),detected.capabilities(),detected.recordCount(),detected.coverage(),metadata,detected.diagnostics());
            String name=Objects.toString(options.getOrDefault("spaceName",Path.of(options.get("input").toString()).getFileName().toString()))+" · "+entry.getKey();
            synchronized(store){var asset=assets.register(store.snapshot().revision,actor,detection,name,library,null,snapshot,"prepared-v1",WorkspaceStore.JSON.writeValueAsString(Map.of("kind",kind,"input",options.get("input"),"detachedJobId",child,"options",options)));
                ids.add(asset.id());library=asset.libraryId();snapshot=asset.sourceSnapshotId();}
        }
        require(!ids.isEmpty(),"Successful worker did not publish recognized prepared assets");return ids;
    }
    public Object exportSlurm(long revision,String actor,String profile,String output)throws Exception {
        Object bundle=DetachedPreparation.exportSlurm(Path.of(profile),Path.of(output));store.mutate(revision,actor,"preparation.slurm_export",s->null);
        return Map.of("exportOnly",true,"submitted",false,"monitored",false,"bundle",bundle);
    }
    private DetachedPreparation driver(Map<String,Object> runtime)throws Exception {return new DetachedPreparation(store.directory().resolve("detached-preparation"),runtime);}
    private static Map<String,Object> map(JsonNode value){return WorkspaceStore.JSON.convertValue(value,new TypeReference<>(){});}
    @SuppressWarnings("unchecked") private static Map<String,Object> cast(Object value){return new LinkedHashMap<>((Map<String,Object>)value);}
}
