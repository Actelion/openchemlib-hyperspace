package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.*;
import com.fasterxml.jackson.core.type.TypeReference;
import com.idorsia.research.chem.hyperspace3d.feature.PoseQueryInput;
import java.nio.file.*;
import java.util.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Shared domain operations for HTTP and MCP. No transport owns chemistry, jobs or state. */
public final class WorkbenchApplication implements AutoCloseable {
    public final WorkspaceStore store;public final AssetService assets;public final QueryService queries;
    public final JobService jobs;public final ResultService results;public final SubstructureService substructure;public final LearnedSearchService learned;public final PreparationService preparation;
    public record Operation(String name,String description,boolean writes,Map<String,Object> arguments){}
    public WorkbenchApplication(Path workspace)throws Exception {
        store=new WorkspaceStore(workspace);
        try{assets=new AssetService(store);queries=new QueryService(store);jobs=new JobService(store,2,8);results=new ResultService(store);
            substructure=new SubstructureService(store,assets,queries,jobs,results);learned=new LearnedSearchService(store,assets,queries,jobs,results);preparation=new PreparationService(store,assets,jobs);
        }catch(Exception e){store.close();throw e;}
    }
    public Object invoke(String operation,JsonNode a,long revision,String actor)throws Exception {
        require(a!=null&&a.isObject(),"Arguments must be an object");
        return switch(operation) {
            case "workspace.state" -> state();
            case "workspace.history" -> history(text(a,"kind"),integer(a,"offset",0),integer(a,"limit",20));
            case "workspace.rename" -> store.mutate(revision,actor,"workspace.rename",s->{String name=text(a,"name");require(name.length()<=200,"Workspace name exceeds 200 characters");s.name=name;return Map.of("workspaceId",s.workspaceId,"name",name);});
            case "query.history" -> queries.history(integer(a,"offset",0),integer(a,"limit",50));
            case "query.branch" -> queries.branch(revision,actor,text(a,"queryId"),optional(a,"name"));
            case "run.branch_query" -> queries.branch(revision,actor,results.getRun(text(a,"runId")).queryId(),optional(a,"name"));
            case "workspace.navigate" -> store.mutate(revision,actor,"workspace.navigate",s->{String view=text(a,"view");require(Set.of("libraries","screening","saved").contains(view),"Unknown supported view");
                String entity=optional(a,"entityId");if(entity!=null)require(s.runs.containsKey(entity)||s.collections.containsKey(entity),"Unknown run or collection");s.view=view;s.viewEntityId=entity;return Map.of("view",view);});
            case "workspace.undo" -> {queries.undo(revision,actor);yield state();}
            case "query.create" -> queries.create(revision,actor,optional(a,"format"),text(a,"structure"),a.hasNonNull("fragment")?a.get("fragment").asBoolean():null,optional(a,"coordinates"),optional(a,"name"),optional(a,"branchFrom"));
            case "query.inspect" -> queries.inspect(text(a,"queryId"));
            case "query.edit" -> queries.edit(revision,actor,text(a,"queryId"),WorkspaceStore.JSON.convertValue(a.get("edits"),new TypeReference<List<QueryChemistry.Edit>>(){}),optional(a,"name"));
            case "query.edit_selection" -> queries.editSelection(revision,actor,QueryChemistry.Operation.valueOf(text(a,"operation")),a.has("elements")?WorkspaceStore.JSON.convertValue(a.get("elements"),new TypeReference<List<Integer>>(){}):null,
                    optional(a,"state"),a.hasNonNull("minAtoms")?a.get("minAtoms").asInt():null,a.hasNonNull("maxAtoms")?a.get("maxAtoms").asInt():null,a.hasNonNull("enabled")?a.get("enabled").asBoolean():null,a.hasNonNull("types")?strings(a,"types"):null,a.path("exclude").asBoolean(false));
            case "query.test" -> QueryChemistry.test(queries.revision(text(a,"queryId")).document(),strings(a,"examples"));
            case "query.from_hit" -> {var hit=results.hit(text(a,"hitId"));yield queries.create(revision,actor,"idcode",hit.document().idcode(),a.path("fragment").asBoolean(true),hit.document().coordinates(),optional(a,"name"),store.snapshot().currentQueryId);}
            case "query.import_pose" -> importPose(revision,actor,a);
            case "query.upload" -> {Object output=upload(a);store.mutate(revision,actor,"query.upload",s->null);yield output;}
            case "selection.change" -> {var selection=WorkspaceStore.JSON.treeToValue(a,WorkspaceStore.Selection.class);results.validateHitIds(selection.hits(),10000);yield queries.select(revision,actor,selection);}
            case "asset.inspect" -> assets.inspect(Path.of(text(a,"path")),a.hasNonNull("kind")?AssetService.Kind.valueOf(a.get("kind").asText()):null);
            case "asset.register" -> assets.register(revision,actor,assets.inspect(Path.of(text(a,"path")),a.hasNonNull("kind")?AssetService.Kind.valueOf(a.get("kind").asText()):null),text(a,"name"),optional(a,"libraryId"),optional(a,"libraryName"),optional(a,"sourceSnapshotId"),optional(a,"sourceVersion"),optional(a,"provenance"));
            case "asset.load" -> loadAsset(revision,actor,text(a,"assetId"),optional(a,"graphModelId"));
            case "asset.unload" -> assets.unload(revision,actor,text(a,"assetId"));
            case "search.substructure" -> substructure.search(revision,actor,text(a,"queryId"),strings(a,"assetIds"),new SubstructureService.Options(integer(a,"threads",2),integer(a,"familyLimit",100),a.path("alignedDepictions").asBoolean(true)));
            case "preparation.environment" -> preparation.environment();
            case "preparation.configure" -> preparation.configure(revision,actor,WorkspaceStore.JSON.convertValue(a.get("runtime"),new TypeReference<Map<String,Object>>(){}));
            case "preparation.inspect" -> preparation.inspect(text(a,"input"),optional(a,"zipEntry"),a.path("molecules").asBoolean());
            case "preparation.recipe" -> preparation.recipe(text(a,"kind"),WorkspaceStore.JSON.convertValue(a.get("options"),new TypeReference<Map<String,Object>>(){}));
            case "preparation.submit" -> preparation.submit(revision,actor,text(a,"kind"),WorkspaceStore.JSON.convertValue(a.get("options"),new TypeReference<Map<String,Object>>(){}));
            case "preparation.resume" -> preparation.resume(revision,actor,text(a,"jobId"));
            case "preparation.attach" -> preparation.attach(revision,actor,text(a,"jobId"));
            case "preparation.log" -> preparation.log(text(a,"jobId"),a.path("stream").asText("worker"),a.path("offset").asLong(0),integer(a,"limit",8192));
            case "preparation.export_slurm" -> preparation.exportSlurm(revision,actor,text(a,"profilePath"),text(a,"outputDirectory"));
            case "search.learned" -> learned.search(revision,actor,WorkspaceStore.JSON.treeToValue(a.get("request"),LearnedSearchService.Request.class));
            case "search.validate" -> learned.validate(WorkspaceStore.JSON.treeToValue(a.get("request"),LearnedSearchService.Request.class));
            case "search.rank_candidates" -> learned.rank(revision,actor,WorkspaceStore.JSON.treeToValue(a.get("request"),LearnedSearchService.Request.class),optional(a,"sourceRunId"),optional(a,"collectionId"),integer(a,"maxCandidates",2000));
            case "search.expand" -> substructure.expand(revision,actor,text(a,"runId"),strings(a,"familyIds"),integer(a,"maxProducts",100),a.path("alignedDepictions").asBoolean(true));
            case "run.families" -> substructure.families(text(a,"runId"),integer(a,"offset",0),integer(a,"limit",50));
            case "run.results" -> results.page(text(a,"runId"),integer(a,"offset",0),integer(a,"limit",50),optional(a,"sort"),a.path("descending").asBoolean(false),a.hasNonNull("filters")?WorkspaceStore.JSON.treeToValue(a.get("filters"),ResultService.Filters.class):null);
            case "run.inspect" -> results.getRun(text(a,"runId"));
            case "hit.inspect" -> results.hit(text(a,"hitId"));
            case "run.compare" -> results.compare(text(a,"leftRunId"),text(a,"rightRunId"),integer(a,"sampleLimit",50));
            case "run.diverse" -> results.diverse(text(a,"runId"),integer(a,"count",20));
            case "collection.save" -> results.saveCollection(revision,actor,optional(a,"collectionId"),text(a,"name"),optional(a,"notes"),a.hasNonNull("hitIds")?strings(a,"hitIds"):null);
            case "collection.results" -> {var rows=results.collectionHits(text(a,"collectionId"));int offset=integer(a,"offset",0),limit=integer(a,"limit",50);require(offset>=0&&offset<=ResultService.MAX_RETAINED&&limit>=1&&limit<=200,"Page limit must be 1–200");yield new ResultService.Page<>(rows.subList(Math.min(offset,rows.size()),Math.min(offset+limit,rows.size())),offset,limit,rows.size(),"Named collection; source occurrences retained.");}
            case "results.export" -> {Object output=results.export(optional(a,"runId"),optional(a,"collectionId"));store.mutate(revision,actor,"results.export",s->null);yield output;}
            case "job.status" -> jobs.get(text(a,"jobId"));
            case "job.cancel" -> jobs.cancel(revision,actor,text(a,"jobId"));
            case "toy.prepare" -> toy(revision,actor);
            default -> throw new WorkbenchException("operation_not_found",404,"Unknown application operation: "+operation);
        };
    }
    private Object loadAsset(long revision,String actor,String id,String graph)throws Exception {
        // Validation and resource construction run in an owner job, not an HTTP request thread.
        assets.get(id);return jobs.submit(revision,actor,"load_asset",null,Map.of("assetId",id),List.of(),ctx->{
            ctx.stage("Loading registered resource");return assets.loadFromJob(actor,id,graph);
        });
    }
    private Object toy(long revision,String actor) {
        return jobs.submit(revision,actor,"prepare_toy",null,Map.of("formatVersion",1),List.of(),ctx->{
            ctx.stage("Building redistributable synthetic index");var paths=ToySetup.prepare(store.directory().resolve("toy"),1);ctx.check();
            synchronized(store){var raw=assets.register(store.snapshot().revision,actor,assets.inspect(Path.of(paths.get("rawspace")),AssetService.Kind.RAW_SPACE),"Toy · raw synthons",null,"Synthetic alkyl library",null,"toy-v1","Generated freely redistributable structures");
                var index=assets.register(store.snapshot().revision,actor,assets.inspect(Path.of(paths.get("substructureIndex")),AssetService.Kind.SUBSTRUCTURE_INDEX),"Toy · substructure",raw.libraryId(),null,raw.sourceSnapshotId(),null,null);
                return Map.of("paths",paths,"rawspaceAssetId",raw.id(),"indexAssetId",index.id());}
        });
    }
    private Object upload(JsonNode args)throws Exception {
        String content=text(args,"content");require(content.length()<=1_000_000,"Query upload exceeds 1 MB");String format=text(args,"format");require(Set.of("mol","sdf").contains(format),"Upload format is mol or sdf");
        Path dir=store.directory().resolve("query-uploads");Files.createDirectories(dir);Path file=dir.resolve(UUID.randomUUID()+"."+format);Files.writeString(file,content);return Map.of("path",file.toString());
    }
    private Object importPose(long revision,String actor,JsonNode args)throws Exception {
        Path input=Path.of(text(args,"path")).toRealPath();require(Files.size(input)<=16*1024*1024,"Query input exceeds 16 MiB; select a smaller supplied-pose file");
        var pose=PoseQueryInput.read(input,args.hasNonNull("record")?args.get("record").asInt():null);
        Path directory=store.directory().resolve("query-poses");Files.createDirectories(directory);Path saved=directory.resolve(pose.provenance().get("inputSha256")+".sdf");
        if(!Files.exists(saved))Files.copy(input,saved);var provenance=new LinkedHashMap<>(pose.provenance());provenance.put("originalPath",input.toString());provenance.put("input",saved.toString());
        return queries.attachPose(revision,actor,QueryChemistry.canonical(pose.molecule()),WorkspaceStore.JSON.valueToTree(provenance),optional(args,"name"));
    }
    public Map<String,Object> state() {
        var s=store.snapshot();Map<String,Object> result=new LinkedHashMap<>();
        result.put("schemaVersion",s.schemaVersion);result.put("workspaceId",s.workspaceId);result.put("name",s.name);result.put("revision",s.revision);result.put("sequence",s.sequence);
        result.put("currentQuery",s.currentQueryId==null?null:s.queries.get(s.currentQueryId));result.put("inspection",s.currentQueryId==null?null:queries.inspect(s.currentQueryId));result.put("selection",s.selection);
        result.put("libraries",s.libraries);result.put("sourceSnapshots",s.sourceSnapshots);result.put("assets",s.assets);result.put("resourceStatus",s.assets.keySet().stream().collect(java.util.stream.Collectors.toMap(id->id,assets::resourceStatus)));
        result.put("runCount",s.runs.size());result.put("collectionCount",s.collections.size());result.put("jobCount",s.jobs.size());
        result.put("runs",visible(s.runs,s.viewEntityId).stream().map(WorkbenchApplication::runSummary).toList());result.put("collections",visible(s.collections,s.viewEntityId).stream().map(WorkbenchApplication::collectionSummary).toList());result.put("jobs",last(s.jobs.values(),100).stream().map(WorkbenchApplication::jobSummary).toList());result.put("activity",last(s.activity,100));
        result.put("queryHistory",last(s.queries.values().stream().map(q->Map.of("id",q.id(),"name",q.name(),"createdAt",q.createdAt())).toList(),100));
        result.put("view",s.view);result.put("viewEntityId",s.viewEntityId);result.put("undoAvailable",!s.undo.isEmpty());result.put("identityPolicy",ResultService.IDENTITY_POLICY);return result;
    }
    /** The selected saved entity stays available even when it is older than the recent summary window. */
    private static List<JsonNode> visible(Map<String,JsonNode> values,String selected){
        var rows=new ArrayList<>(last(values.values(),100));JsonNode chosen=values.get(selected);
        if(chosen!=null&&rows.stream().noneMatch(n->n.path("id").asText().equals(selected)))rows.add(chosen);return rows;
    }
    public ResultService.Page<JsonNode> history(String kind,int offset,int limit){
        require(offset>=0&&limit>=1&&limit<=100,"History offset must be nonnegative and page limit 1–100");var s=store.snapshot();
        List<JsonNode> values=switch(kind){
            case "runs" -> s.runs.values().stream().map(WorkbenchApplication::runSummary).toList();
            case "collections" -> s.collections.values().stream().map(WorkbenchApplication::collectionSummary).toList();
            case "jobs" -> s.jobs.values().stream().map(WorkbenchApplication::jobSummary).toList();
            default -> throw new WorkbenchException("invalid_argument",400,"History kind must be runs, collections or jobs");
        };
        var rows=new ArrayList<>(values);Collections.reverse(rows);int start=Math.min(offset,rows.size()),end=(int)Math.min((long)start+limit,rows.size());
        return new ResultService.Page<>(List.copyOf(rows.subList(start,end)),offset,limit,rows.size(),"Persistent workspace history, newest first; results and collection members are paged separately.");
    }
    private static JsonNode runSummary(JsonNode node) {
        var result=node.deepCopy();if(result.path("settings") instanceof com.fasterxml.jackson.databind.node.ObjectNode settings){
            if(settings.has("sourceOccurrenceIds"))settings.put("sourceOccurrenceCount",settings.path("sourceOccurrenceIds").size());settings.remove(List.of("sourceOccurrenceIds","sourceRuns"));
        }return result;
    }
    private static JsonNode collectionSummary(JsonNode node){var result=(com.fasterxml.jackson.databind.node.ObjectNode)node.deepCopy();result.put("hitCount",result.path("hitIds").size());result.remove("hitIds");return result;}
    private static JsonNode jobSummary(JsonNode node){var result=(com.fasterxml.jackson.databind.node.ObjectNode)node.deepCopy();if(result.path("request") instanceof com.fasterxml.jackson.databind.node.ObjectNode request)request.remove(List.of("sourceOccurrenceIds","sourceRuns"));
        var lines=result.putArray("log");var log=node.path("log");for(int i=Math.max(0,log.size()-10);i<log.size();i++)lines.add(log.get(i));return result;}
    private static <T> List<T> last(java.util.Collection<T> values,int limit){var list=new ArrayList<>(values);return List.copyOf(list.subList(Math.max(0,list.size()-limit),list.size()));}
    public static List<Operation> operations(){return List.of(
        op("workspace.state","Attach/discover this owner workspace, current query, last 100 jobs/runs/activity, asset capabilities and shared selection.",false),
        op("workspace.history","Page all saved run, collection or job summaries newest first; kind runs/collections/jobs, limit 1–100. Large results remain separate.",false,"kind:string","offset:integer?","limit:integer?"),
        op("workspace.rename","Name this persistent workspace without changing its stable ID or filesystem location.",true,"name:string"),
        op("query.history","Page immutable query revisions newest first (limit 1–100), including parent references and original pose provenance.",false,"offset:integer?","limit:integer?"),
        op("query.branch","Publish a new visible branch from a previous query revision, retaining its graph/pose. Original queries/runs remain immutable.",true,"queryId:string","name:string?"),
        op("run.branch_query","Publish a new branch from a saved run's original query, preserving original graph/pose settings and leaving the run unchanged.",true,"runId:string","name:string?"),
        op("workspace.navigate","Open libraries, screening or saved; optional entityId selects a retained run/collection.",true,"view:string","entityId:string?"),
        op("workspace.undo","Undo a query or meaningful selection edit; immutable run computations are unaffected.",true),
        op("query.create","Create/publish a canonical OCL query revision. IDCode preserves atom/bond features; fragment true enables advanced queries.",true,"format:string?","structure:string","fragment:boolean?","coordinates:string?","name:string?","branchFrom:string?"),
        op("query.inspect","Read atoms/bonds and revision-scoped refs and bridge values. Never invent IDCodes.",false,"queryId:string"),
        op("query.edit","Publish typed edits: operation, refs, elements/exclude/state/types/enabled/minAtoms/maxAtoms. References must belong to the visible revision.",true,"queryId:string","edits:array","name:string?"),
        op("query.edit_selection","Apply an ATOM_*/BOND_* typed operation to scientist-selected refs. Bridge min is 0–15, max-min 0–15; 0..0 clears the encoding and is rejected.",true,"operation:string","elements:array?","state:string?","minAtoms:integer?","maxAtoms:integer?","enabled:boolean?","types:array?","exclude:boolean?"),
        op("query.test","Test current advanced OCL query against at most 50 concrete SMILES examples.",false,"queryId:string","examples:array"),
        op("query.from_hit","Publish a new query from one retained hit occurrence.",true,"hitId:string","fragment:boolean?","name:string?"),
        op("query.import_pose","Attach supplied MOL/SDF coordinates. Explicit zero-based record required for multi-record input; graph edits detach this pose.",true,"path:string","record:integer?","name:string?"),
        op("query.upload","Save a small MOL/SDF upload (1 MB maximum) for explicit pose record selection.",true,"content:string","format:string"),
        op("selection.change","Share current revision's selected atoms/bonds and retained occurrence IDs; stale refs are rejected.",true,"queryRevisionId:string?","atoms:array","bonds:array","hits:array"),
        op("asset.inspect","Inspect local prepared-asset path/manifest without loading large indexes or scanning fingerprint payloads.",false,"path:string","kind:string?"),
        op("asset.register","Register a prepared local asset, optionally linking an explicit library/source snapshot. Registration does not load it.",true,"path:string","kind:string?","name:string","libraryId:string?","libraryName:string?","sourceSnapshotId:string?","sourceVersion:string?","provenance:string?"),
        op("asset.load","Explicitly load an asset in an owner job; compact/pose bundles require graphModelId. No silent device substitution.",true,"assetId:string","graphModelId:string?"),
        op("asset.unload","Unload an asset; active or queued jobs holding leases prevent unloading.",true,"assetId:string"),
        op("search.substructure","Start real headless combinatorial search. Retains bounded families, with overlapping represented counts, not expanded products.",true,"queryId:string","assetIds:array","threads:integer?","familyLimit:integer?","alignedDepictions:boolean?"),
        op("preparation.environment","Read explicit detached-worker Java/CLI/fingerprint/model paths and resource defaults.",false),
        op("preparation.configure","Persist runtime paths and default threads/child heap. Runtime object has formatVersion=1, cliJar, fingerprintJar, fingerprintLibDirectory, modelBundle, compactBundle, javaExecutable, threads, heap.",true,"runtime:object"),
        op("preparation.inspect","Inspect bounded supplier input/ZIP entry or enumerated molecule-table sample. No large browser upload required; samples are data, not instructions.",false,"input:string","zipEntry:string?","molecules:boolean?"),
        op("preparation.recipe","Validate and propose existing-worker recipe: kind space (options input,format,spaceName,searchModes,importOptions,threads,heap) or fingerprints (options input,outputDirectory,smilesColumn,idColumn,device,threads,heap,encoderBatchSize,cudaDeviceId,sourceRowsPerShard). Returns produced coverage, dependencies, runtime and resume behavior.",false,"kind:string","options:object"),
        op("preparation.submit","Start the validated detached recipe and monitor it in the owner job drawer. Successful rawspace/index/cache outputs auto-register with explicit provenance; no desktop launch.",true,"kind:string","options:object"),
        op("preparation.resume","Recover a terminal owner preparation using its pinned runtime/options. Reattach if the recorded worker is live/succeeded; otherwise restart fingerprint shards at the same output path or reuse a completed rawspace import. Never blindly restart a live process.",true,"jobId:string"),
        op("preparation.attach","Reattach an interrupted owner preparation to its recorded live/completed detached worker; verify status first. Idempotent registration of completed outputs.",true,"jobId:string"),
        op("preparation.log","Read up to 65536 bytes of worker/import/build/fingerprint log for an owner preparation job.",false,"jobId:string","stream:string?","offset:integer?","limit:integer?"),
        op("preparation.export_slurm","Export existing Slurm fingerprint bundle only. No submission, SSH or cluster monitoring.",true,"profilePath:string","outputDirectory:string"),
        op("search.learned","Screen a loaded native/Python flat cache with a loaded canonical graph model: request has queryId, indexAssetId, modelAssetId, mode (SKELSPHERES_2D/PHESA_3D/SKELSPHERES_THEN_PHESA/POSE_PHESA_3D), optional compactAssetId/poseAssetId, objective, topK, shortlistTopK, exactSkelSpheres, device CPU/CUDA, intraOpThreads, batchSize, prefetchDepth, posePpWeight. Strict model/cache validation before scanning.",true,"request:object"),
        op("search.validate","Validate query and loaded model/cache provenance without scanning payloads or opening ONNX sessions. Same typed request as search.learned.",false,"request:object"),
        op("search.rank_candidates","Freeze a bounded retained run/collection prefix, featurize/encode unique graphs, and rank every accepted source occurrence with graph or pose PheSA. Preserves synthon provenance and records rejection/expansion counts; no exhaustive hit-family ranking claim.",true,"request:object","sourceRunId:string?","collectionId:string?","maxCandidates:integer?"),
        op("search.expand","Assemble a bounded deterministic prefix of selected retained families; maxProducts is a global attempts cap, not exhaustive ranking.",true,"runId:string","familyIds:array","maxProducts:integer?","alignedDepictions:boolean?"),
        op("run.families","Page retained combinatorial families; maximum 200 rows per page.",false,"runId:string","offset:integer?","limit:integer?"),
        op("run.results","Page/filter/sort retained result occurrences only; filtering top-K does not imply globally optimal filtered top-K.",false,"runId:string","offset:integer?","limit:integer?","sort:string?","descending:boolean?","filters:object?"),
        op("run.inspect","Inspect immutable query/settings, coverage, counts and diagnostics for a saved run.",false,"runId:string"),
        op("hit.inspect","Inspect one occurrence, graph, matching atoms, synthons, raw score components and source provenance.",false,"hitId:string"),
        op("run.compare","Compare canonical molecular identities between two retained runs; never average incompatible method scores.",false,"leftRunId:string","rightRunId:string","sampleLimit:integer?"),
        op("run.diverse","Select up to 500 diverse unique molecular identities using OCL's FFP512 DiversitySelector.",false,"runId:string","count:integer?"),
        op("collection.save","Create/update named collection and notes with explicit occurrence IDs; preserves source occurrences.",true,"collectionId:string?","name:string","notes:string?","hitIds:array?"),
        op("collection.results","Page a named collection, at most 200 source occurrences per page.",false,"collectionId:string","offset:integer?","limit:integer?"),
        op("results.export","Export one run or collection to SDF, TSV, Markdown and JSON provenance inside this workspace.",true,"runId:string?","collectionId:string?"),
        op("job.status","Read persistent status and last 100 bounded log lines; browser disconnect does not cancel work.",false,"jobId:string"),
        op("job.cancel","Request cooperative cancellation; resources stay leased until all underlying work returns.",true,"jobId:string"),
        op("toy.prepare","Build and register redistributable synthetic rawspace and a real combinatorial substructure index.",true));}
    private static Operation op(String name,String description,boolean writes,String... fields) {
        Map<String,Object> properties=new LinkedHashMap<>();List<String> required=new ArrayList<>();for(String field:fields){String[] parts=field.split(":");boolean optional=parts[1].endsWith("?");String type=parts[1].replace("?","");
            var schema=new LinkedHashMap<String,Object>();schema.put("type",optional?List.of(type,"null"):type);if(type.equals("array"))schema.put("items",Map.of());properties.put(parts[0],schema);if(!optional)required.add(parts[0]);}
        return new Operation(name,description,writes,Map.of("type","object","properties",properties,"required",required,"additionalProperties",false));
    }
    private static String text(JsonNode node,String field){String value=optional(node,field);require(value!=null&&!value.isBlank(),"Missing argument: "+field);return value;}
    private static String optional(JsonNode node,String field){return node.hasNonNull(field)?node.get(field).asText():null;}
    private static int integer(JsonNode node,String field,int fallback){if(!node.hasNonNull(field))return fallback;require(node.get(field).canConvertToInt(),"Invalid integer: "+field);return node.get(field).asInt();}
    private static List<String> strings(JsonNode node,String field){require(node.path(field).isArray(),"Expected array: "+field);List<String> result=new ArrayList<>();node.get(field).forEach(n->{require(n.isTextual(),"Expected string entries: "+field);result.add(n.asText());});return result;}
    @Override public void close()throws Exception{jobs.close();try{assets.close();}finally{store.close();}}
}
