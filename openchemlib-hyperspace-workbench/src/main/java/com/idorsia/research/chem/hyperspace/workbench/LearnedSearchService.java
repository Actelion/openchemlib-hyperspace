package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import com.idorsia.research.chem.hyperspace3d.screening.*;
import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Existing format-neutral screeners using owner-leased caches and shared native model sessions. */
public final class LearnedSearchService {
    public enum Mode {SKELSPHERES_2D,PHESA_3D,SKELSPHERES_THEN_PHESA,POSE_PHESA_3D}
    public static final class Request {
        public int formatVersion=1;public String queryId,indexAssetId,modelAssetId,compactAssetId,poseAssetId;
        public Mode mode=Mode.PHESA_3D;public String objective="phesa_total";public Double compositePpWeight;
        public double posePpWeight=0.5;public boolean exactSkelSpheres=true;
        public int encoderBatchSize=128;public int topK=100,shortlistTopK=2000,batchSize=4096,scanBatchSize=65536,prefetchDepth=2;
        public String device="CPU";public int cudaDeviceId,intraOpThreads=2;
        public void validate(){require(formatVersion==1&&mode!=null,"Unsupported learned search request");require(queryId!=null&&modelAssetId!=null,"Query and canonical model asset required");
            require(Set.of("phesa_total","phesa_shape","phesa_pharmacophore").contains(objective),"Choose total, shape or pharmacophore objective");
            require(topK>=1&&topK<=ResultService.MAX_RETAINED,"Retained topK must be 1–20000");require(shortlistTopK>=topK&&shortlistTopK<=200000,"Shortlist must be >=topK and <=200000");
            require(encoderBatchSize>=1&&encoderBatchSize<=1024,"Encoder batch size must be 1–1024");
            require(batchSize>=1&&batchSize<=65536&&scanBatchSize>=1&&scanBatchSize<=1000000&&prefetchDepth>=1&&prefetchDepth<=8,"Invalid bounded batch/prefetch settings");
            require(Set.of("CPU","CUDA").contains(device)&&cudaDeviceId>=0&&intraOpThreads>=1&&intraOpThreads<=64,"Select explicit CPU/CUDA runtime and 1–64 threads");
            require(Double.isFinite(posePpWeight)&&posePpWeight>=0&&posePpWeight<=1,"Pose PP weight must be in [0,1]");
            require(compositePpWeight==null||(Double.isFinite(compositePpWeight)&&compositePpWeight>=0&&compositePpWeight<=1),"Composite PP weight must be in [0,1]");
            require(mode!=Mode.POSE_PHESA_3D||compositePpWeight==null,"Pose totals use the pose comparator's configured weight; do not add an independent composite");
        }
        public boolean compact(){return mode==Mode.SKELSPHERES_2D||mode==Mode.SKELSPHERES_THEN_PHESA;}
        public boolean pose(){return mode==Mode.POSE_PHESA_3D;}
        public String runtimeKey(){return device+":"+cudaDeviceId+":"+intraOpThreads;}
    }
    static final class GraphSessions implements AutoCloseable {
        final DeepSpaceV1Encoder encoder;final DeepSpaceV1Comparator comparator;
        GraphSessions(DeepSpaceOnnxEnvironment runtime,DeepSpaceModelBundle bundle){encoder=new DeepSpaceV1Encoder(runtime,bundle);
            try{comparator=new DeepSpaceV1Comparator(runtime,bundle);}catch(RuntimeException e){encoder.close();throw e;}}
        @Override public void close(){try{comparator.close();}finally{encoder.close();}}
    }
    private record Inputs(WorkspaceStore.QueryRevision query,DeepSpaceModelBundle model,CompactSkelSpheresModelBundle compact,
                          SharedCandidatePoseBundle pose,PoseQueryInput supplied,StereoMolecule graph,MoleculeFingerprintDataSource source,List<AssetService.Lease> leases){}
    private final WorkspaceStore store;private final AssetService assets;private final QueryService queries;private final JobService jobs;private final ResultService results;
    public LearnedSearchService(WorkspaceStore store,AssetService assets,QueryService queries,JobService jobs,ResultService results){this.store=store;this.assets=assets;this.queries=queries;this.jobs=jobs;this.results=results;}
    public JobService.Job search(long revision,String actor,Request submitted)throws Exception {
        Request request=WorkspaceStore.JSON.convertValue(submitted,Request.class);request.validate();require(request.indexAssetId!=null,"Select a flat fingerprint cache");Inputs input=acquire(request,true);
        try{synchronized(store){String scope=switch(request.mode){
            case SKELSPHERES_2D->request.exactSkelSpheres?"All compact vectors scanned; exact OCL SkelSpheres ranking only within the learned top-"+request.shortlistTopK+" shortlist.":"Exhaustive learned 16D SkelSpheres scan over accepted cache rows.";
            case SKELSPHERES_THEN_PHESA->"All compact vectors scanned; learned PheSA ranking only within the 16D top-"+request.shortlistTopK+" shortlist.";
            case PHESA_3D->"Exhaustive learned graph-to-ensemble PheSA ranking over accepted 128D cache rows.";
            case POSE_PHESA_3D->"Exhaustive learned directed supplied-pose → graph-only ensemble predictions over accepted cache rows. No candidate conformers or computed 3D alignment.";};
            var run=results.create(revision,actor,request.mode.name().toLowerCase(Locale.ROOT),input.query(),assetIds(request),request,scope);
            return jobs.submit(store.snapshot().revision,actor,"learned_screen",run.id(),request,input.leases(),ctx->screen(ctx,run,request,input));
        }}catch(Exception e){close(input.leases());throw e;}
    }
    public Map<String,Object> validate(Request request)throws Exception {
        request.validate();Inputs i=acquire(request,request.indexAssetId!=null);
        try {var out=new LinkedHashMap<String,Object>();out.put("compatible",true);out.put("modelBundleHash",i.model().bundleHash());
            out.put("scoreSemantics",request.pose()?"learned_pose_to_ensemble":"graph_only");if(i.source()!=null)out.put("acceptedCacheRows",i.source().recordCount());
            return out;
        }finally{close(i.leases());}
    }
    /** Rank a frozen bounded list of already assembled products; provenance/occurrences are never collapsed. */
    public JobService.Job rank(long revision,String actor,Request submitted,String sourceRun,String collection,int maximum)throws Exception {
        Request request=WorkspaceStore.JSON.convertValue(submitted,Request.class);request.validate();require(request.mode==Mode.PHESA_3D||request.mode==Mode.POSE_PHESA_3D,"Bounded ranking supports graph/pose PheSA");
        require((sourceRun==null)!=(collection==null),"Choose exactly one source run or collection");require(maximum>=1&&maximum<=ResultService.MAX_RETAINED,"Candidate cap must be 1–20000");
        List<ResultService.Hit> all=sourceRun!=null?results.all(sourceRun):results.collectionHits(collection);
        List<ResultService.Hit> selected=List.copyOf(all.subList(0,Math.min(maximum,all.size())));require(!selected.isEmpty(),"No concrete candidate occurrences to rank; expand families first");
        Inputs input=acquire(request,false);try{synchronized(store){
            var settings=new LinkedHashMap<String,Object>();settings.put("request",request);settings.put("sourceRunId",sourceRun);settings.put("collectionId",collection);settings.put("maxCandidates",maximum);
            settings.put("sourceOccurrenceIds",selected.stream().map(ResultService.Hit::id).toList());settings.put("sourceAvailableCount",all.size());
            Map<String,Object> sourceRuns=new LinkedHashMap<>();for(var hit:selected)if(!sourceRuns.containsKey(hit.runId())) {
                var source=results.getRun(hit.runId());var provenance=new LinkedHashMap<String,Object>();provenance.put("id",source.id());provenance.put("mode",source.mode());provenance.put("queryId",source.queryId());provenance.put("assetIds",source.assetIds());provenance.put("scope",source.scope());provenance.put("diagnostics",source.diagnostics());
                if(source.mode().equals("bounded_expansion"))provenance.put("expansionSettings",source.settings());
                sourceRuns.put(hit.runId(),provenance);
            }
            settings.put("sourceRuns",sourceRuns);
            settings.put("policy","Frozen retained-source order prefix; one embedding per unique accepted canonical graph; each source occurrence remains separate.");
            Set<String> ids=new LinkedHashSet<>(assetIds(request));selected.forEach(h->ids.add(h.assetId()));
            var run=results.create(revision,actor,"bounded_learned_ranking",input.query(),List.copyOf(ids),settings,
                    "Learned "+(request.pose()?"directed pose → ensemble":"graph → ensemble")+" ranking of at most "+maximum+" retained candidate occurrences. Not exhaustive over represented synthon combinations.");
            return jobs.submit(store.snapshot().revision,actor,"bounded_learned_ranking",run.id(),Map.of("runId",run.id(),"maxCandidates",maximum,"request",request),input.leases(),ctx->rank(ctx,run,request,input,selected,all.size()));
        }}catch(Exception e){close(input.leases());throw e;}
    }
    private Object rank(JobService.Context ctx,ResultService.Run run,Request request,Inputs input,List<ResultService.Hit> source,int available)throws Exception {
        ctx.stage("Encoding query for bounded ranking");float[] query=query(request,input);var sessions=graphSessions(request,input);
        var comparator=request.pose()?(EmbeddingComparator)poseSessions(request,input):sessions.comparator;
        var manifest=request.pose()?input.pose().scoreManifest():input.model().manifest();
        Map<String,float[]> scores=new HashMap<>();Map<String,String> rejected=new HashMap<>();Map<String,ResultService.Hit> unique=new LinkedHashMap<>();source.forEach(h->unique.putIfAbsent(h.moleculeId(),h));
        List<Map.Entry<String,ResultService.Hit>> molecules=new ArrayList<>(unique.entrySet());var featurizer=new OCLDeepSpaceFeaturizer();var builder=new DeepSpaceTensorBatchBuilder(featurizer);long preparation=0,encoding=0,comparison=0;
        ctx.stage("Featurizing/encoding bounded candidates on "+request.device);
        for(int start=0;start<molecules.size();start+=request.encoderBatchSize) {
            ctx.check();List<FeaturizationResult> features=new ArrayList<>();List<String> identities=new ArrayList<>();long t=System.nanoTime();
            for(var entry:molecules.subList(start,Math.min(start+request.encoderBatchSize,molecules.size()))) {
                var f=featurizer.featurize(QueryChemistry.parse(entry.getValue().document()));
                if(!f.accepted())rejected.put(entry.getKey(),f.rejectionReason());else{features.add(f);identities.add(entry.getKey());}
            }
            if(features.isEmpty()){preparation+=System.nanoTime()-t;continue;}var batch=builder.buildFromFeatures(features);preparation+=System.nanoTime()-t;
            ctx.check();t=System.nanoTime();float[][] embeddings=sessions.encoder.encode(batch);encoding+=System.nanoTime()-t;
            ctx.check();t=System.nanoTime();float[][] predicted=comparator.compare(query,embeddings);comparison+=System.nanoTime()-t;
            for(int n=0;n<identities.size();n++)scores.put(identities.get(n),predicted[n]);
        }
        List<ResultService.Hit> hits=new ArrayList<>();for(var old:source) {
            ctx.check();float[] predicted=scores.get(old.moleculeId());if(predicted==null)continue;Map<String,Double> values=new LinkedHashMap<>();
            for(int n=0;n<manifest.availableTargets.size();n++)values.put(manifest.availableTargets.get(n),(double)predicted[n]);values.put("objective_score",objective(request).score(predicted,manifest));
            var provenance=WorkspaceStore.JSON.createObjectNode();provenance.put("sourceOccurrenceId",old.id());provenance.put("sourceRunId",old.runId());provenance.set("sourceProvenance",old.provenance());provenance.set("sourceScores",WorkspaceStore.JSON.valueToTree(old.scores()));
            provenance.put("modelBundleHash",input.model().bundleHash());provenance.put("scoreSemantics",request.pose()?"learned_pose_to_ensemble":"graph_only");
            hits.add(new ResultService.Hit(run.id()+"/hit-"+ResultService.hash(old.id()),run.id(),old.moleculeId(),old.document(),old.smiles(),old.assetId(),old.familyId(),old.reaction(),old.synthons(),values,old.properties(),old.matchedAtoms(),provenance));
        }
        hits.sort(Comparator.comparingDouble((ResultService.Hit h)->h.scores().get("objective_score")).reversed().thenComparing(ResultService.Hit::id));
        var diagnostics=new LinkedHashMap<String,Object>();diagnostics.put("availableSourceOccurrences",available);diagnostics.put("attemptedCandidates",source.size());diagnostics.put("uniqueCandidateGraphs",unique.size());
        diagnostics.put("objective",request.compositePpWeight==null?request.objective:"explicit_shape_pp_composite");diagnostics.put("acceptedCandidateGraphs",scores.size());diagnostics.put("rejectedCandidateGraphs",rejected.size());diagnostics.put("scoredCandidateOccurrences",hits.size());diagnostics.put("rejectedCandidateOccurrences",source.stream().filter(h->rejected.containsKey(h.moleculeId())).count());
        diagnostics.put("rejectionReasons",rejected.values().stream().collect(java.util.stream.Collectors.groupingBy(v->v,java.util.stream.Collectors.counting())));
        diagnostics.put("modelBundleHash",input.model().bundleHash());diagnostics.put("scoreSemantics",request.pose()?"learned_pose_to_ensemble":"graph_only");diagnostics.put("timings",Map.of("cpuPreparationNanos",preparation,"encodingNanos",encoding,"comparatorNanos",comparison));
        if(input.pose()!=null){diagnostics.put("poseBundleHash",input.pose().hash());diagnostics.put("posePpWeight",request.posePpWeight);}
        ctx.check();results.complete(run.id(),hits.subList(0,Math.min(request.topK,hits.size())),diagnostics);return diagnostics;
    }
    private Inputs acquire(Request r,boolean cache)throws Exception {
        List<AssetService.Lease> leases=new ArrayList<>();try {
            var query=queries.revision(r.queryId);var model=resource(r.modelAssetId,AssetService.Kind.GRAPH_MODEL,DeepSpaceModelBundle.class,leases);
            CompactSkelSpheresModelBundle compact=null;SharedCandidatePoseBundle pose=null;PoseQueryInput supplied=null;StereoMolecule graph=null;
            if(r.compact()){compact=resource(r.compactAssetId,AssetService.Kind.COMPACT_MODEL,CompactSkelSpheresModelBundle.class,leases);compact.manifest().validateSource(model.manifest());}
            if(r.pose()) {
                pose=resource(r.poseAssetId,AssetService.Kind.POSE_MODEL,SharedCandidatePoseBundle.class,leases);
                // Validate against this request's actual graph encoder, even if the pose asset was loaded against a different model alias.
                require(pose.manifest().path("requiredCandidateEncoderSha256").asText().equals(model.manifest().encoderSha256)
                        &&pose.manifest().path("canonicalPredictorSha256").asText().equals(model.manifest().predictorCheckpointSha256)
                        &&pose.manifest().path("foundationCheckpointSha256").asText().equals(model.manifest().foundationCheckpointSha256),"Pose and canonical model checkpoint provenance differs");
                require(query.pose()!=null&&!query.pose().isNull(),"Attach one validated supplied MOL/SDF conformation before pose screening");
                Path path=Path.of(query.pose().path("input").asText());
                require(SharedCandidatePoseBundle.sha256(path).equals(query.pose().path("inputSha256").asText()),"Stored query pose checksum changed");
                supplied=PoseQueryInput.read(path,query.pose().path("record").asInt());
                require(QueryChemistry.canonical(supplied.molecule()).idcode().equals(query.document().idcode()),"Stored pose graph differs from this query revision");
            } else graph=concrete(query.document());
            MoleculeFingerprintDataSource source=null;if(cache){
                var kind=assets.get(r.indexAssetId).kind();require(kind==AssetService.Kind.FLAT_INDEX||kind==AssetService.Kind.PYTHON_CACHE,"Select a native/Python flat-molecule fingerprint cache");
                var lease=assets.acquire(r.indexAssetId);leases.add(lease);require(lease.value() instanceof MoleculeFingerprintDataSource,"Unsupported loaded cache resource");
                source=(MoleculeFingerprintDataSource)lease.value();source.validateCompatibility(model,compact,r.compact());
            }
            return new Inputs(query,model,compact,pose,supplied,graph,source,leases);
        }catch(Exception e){close(leases);throw e;}
    }
    private <T>T resource(String id,AssetService.Kind kind,Class<T> type,List<AssetService.Lease> leases){require(id!=null,"Missing required "+kind+" asset; register and explicitly load its model bundle");require(assets.get(id).kind()==kind,"Wrong asset kind; expected "+kind);var lease=assets.acquire(id);leases.add(lease);require(type.isInstance(lease.value()),"Wrong loaded resource type");return type.cast(lease.value());}
    static StereoMolecule concrete(QueryChemistry.Document document){
        StereoMolecule m=QueryChemistry.parse(document);for(int a=0;a<m.getAllAtoms();a++)require(m.getAtomQueryFeatures(a)==0&&m.getAtomList(a)==null,"Learned scoring requires a concrete graph without atom query features");
        for(int b=0;b<m.getAllBonds();b++)require(m.getBondQueryFeatures(b)==0,"Learned scoring requires a concrete graph without bond query features");
        m.setFragment(false);var result=new OCLDeepSpaceFeaturizer().featurize(m);require(result.accepted(),"Query rejected by Deepspace featurization: "+result.rejectionReason());return m;
    }
    private DeepSpaceOnnxEnvironment runtime(Request r){return new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.valueOf(r.device),r.cudaDeviceId,r.intraOpThreads);}
    private GraphSessions graphSessions(Request r,Inputs i)throws Exception{return assets.attachment(r.modelAssetId,"graph:"+r.runtimeKey(),GraphSessions.class,()->new GraphSessions(runtime(r),i.model()));}
    private SharedCandidatePoseRuntime poseSessions(Request r,Inputs i)throws Exception{return assets.attachment(r.poseAssetId,"pose:"+r.runtimeKey()+":"+r.posePpWeight,SharedCandidatePoseRuntime.class,()->new SharedCandidatePoseRuntime(runtime(r),i.pose(),r.posePpWeight));}
    private float[] query(Request r,Inputs i)throws Exception {
        if(r.pose())return poseSessions(r,i).encode(i.supplied().batch(),i.supplied().distances())[0];
        return graphSessions(r,i).encoder.encode(new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(i.graph())))[0];
    }
    private ScreeningObjective objective(Request r){return r.compositePpWeight==null?ScreeningObjective.direct(r.objective):ScreeningObjective.phesaComposite(1-r.compositePpWeight,r.compositePpWeight);}
    private Object screen(JobService.Context ctx,ResultService.Run run,Request r,Inputs input)throws Exception {
        results.state(run.id(),"running",Map.of());ctx.stage("Encoding query on "+r.device);long started=System.nanoTime();float[] q=query(r,input);long queryNanos=System.nanoTime()-started;
        ctx.check();var source=guard(input.source(),ctx);List<ResultService.Hit> rows=new ArrayList<>();long scanned,compared=0;Map<String,Object> timings=new LinkedHashMap<>();
        if(r.mode==Mode.SKELSPHERES_2D||r.mode==Mode.SKELSPHERES_THEN_PHESA){
            var projector=assets.attachment(r.compactAssetId,"projection:"+r.runtimeKey(),CompactSkelSpheresProjector.class,()->new CompactSkelSpheresProjector(runtime(r),input.compact()));
            float[] q16=projector.project(q);ctx.stage("Scanning learned 16D SkelSpheres vectors");var screener=new MoleculeSkelSpheresScreener(new CompactSkelSpheresScorer(input.compact().manifest()));
            if(r.mode==Mode.SKELSPHERES_2D){
                var compact=screener.screen(source,q16,r.exactSkelSpheres?r.shortlistTopK:r.topK,r.scanBatchSize,r.prefetchDepth);scanned=compact.recordsScanned();
                List<MoleculeSkelSpheresHit> hits=compact.hits();ctx.check();if(r.exactSkelSpheres){ctx.stage("Exact SkelSpheres reranking of learned shortlist");hits=new MoleculeExactSkelSpheresReranker().rerank(input.graph(),hits);}
                for(var h:hits.stream().limit(r.topK).toList()){ctx.check();Map<String,Double> scores=new LinkedHashMap<>();scores.put("predicted_skelspheres",h.predictedSimilarity());scores.put("compact_dot_product",h.dotProduct());
                    if(h.exactSimilarity()!=null)scores.put("exact_skelspheres",h.exactSimilarity());scores.put("objective_score",r.exactSkelSpheres?h.exactSimilarity():h.predictedSimilarity());rows.add(hit(run,r,h.reference(),h.molecule(),scores,Map.of("exactRerankScope",r.exactSkelSpheres?"learned_shortlist":"disabled")));}
                timings.put("compactScanNanos",compact.scanNanos());timings.put("metadataResolutionNanos",compact.metadataResolutionNanos());
            }else{
                var compact=screener.scanUnresolved(source,q16,r.shortlistTopK,r.scanBatchSize,r.prefetchDepth);scanned=compact.recordsScanned();ctx.stage("PheSA reranking of 16D shortlist");
                var result=new MoleculePheSAScreener(guard(graphSessions(r,input).comparator,ctx),input.model().manifest(),objective(r)).rerank(source,q,compact.hits(),r.topK,r.batchSize);
                compared=result.recordsScanned();render3d(rows,result,run,r,input.model().manifest());timings.put("compactScanNanos",compact.scanNanos());timings.put("comparatorNanos",result.comparatorNanos());
            }
        }else{
            var manifest=r.pose()?input.pose().scoreManifest():input.model().manifest();EmbeddingComparator comparator=r.pose()?poseSessions(r,input):graphSessions(r,input).comparator;
            ctx.stage(r.pose()?"Screening directed pose → ensemble scores":"Screening learned ensemble PheSA scores");
            var result=new MoleculePheSAScreener(guard(comparator,ctx),manifest,objective(r)).screen(source,q,r.topK,r.batchSize,r.prefetchDepth);
            scanned=result.recordsScanned();compared=scanned;render3d(rows,result,run,r,manifest);timings.put("comparatorNanos",result.comparatorNanos());timings.put("vectorPreparationNanos",result.vectorPreparationNanos());timings.put("metadataResolutionNanos",result.metadataResolutionNanos());
        }
        ctx.check();timings.put("queryEncodingNanos",queryNanos);timings.put("totalNanos",System.nanoTime()-started);
        var diagnostics=new LinkedHashMap<String,Object>();diagnostics.put("recordsScanned",scanned);diagnostics.put("libraryRecords",source.recordCount());diagnostics.put("comparatorCandidates",compared);
        diagnostics.put("scope",run.scope());diagnostics.put("modelBundleHash",input.model().bundleHash());diagnostics.put("scoreSemantics",r.pose()?"learned_pose_to_ensemble":"graph_only");
        diagnostics.put("objective",r.mode==Mode.SKELSPHERES_2D?(r.exactSkelSpheres?"exact_skelspheres":"predicted_skelspheres"):(r.compositePpWeight==null?r.objective:"explicit_shape_pp_composite"));diagnostics.put("timings",timings);diagnostics.put("runtime",Map.of("device",r.device,"cudaDeviceId",r.cudaDeviceId,"intraOpThreads",r.intraOpThreads));
        if(input.compact()!=null)diagnostics.put("compactBundleHash",input.compact().bundleHash());if(input.pose()!=null){diagnostics.put("poseBundleHash",input.pose().hash());diagnostics.put("queryPose",input.query().pose());diagnostics.put("posePpWeight",r.posePpWeight);}
        results.complete(run.id(),rows,diagnostics);return diagnostics;
    }
    private void render3d(List<ResultService.Hit> rows,MoleculePheSAScreeningResult result,ResultService.Run run,Request request,DeepSpaceModelManifest manifest)throws Exception {
        for(var h:result.hits()){var scores=new LinkedHashMap<String,Double>();for(int n=0;n<manifest.availableTargets.size();n++)scores.put(manifest.availableTargets.get(n),(double)h.modelScores()[n]);scores.put("objective_score",h.objectiveScore());if(h.predictedSkelSpheres()!=null)scores.put("predicted_skelspheres",h.predictedSkelSpheres());
            var extra=new LinkedHashMap<String,Object>();if(h.compactRank()!=null)extra.put("compactRank",h.compactRank());rows.add(hit(run,request,h.reference(),h.molecule(),scores,extra));}
    }
    private ResultService.Hit hit(ResultService.Run run,Request r,MoleculeVectorReference reference,MoleculeFingerprintMetadata metadata,Map<String,Double> scores,Map<String,Object> extra)throws Exception {
        var molecule=new StereoMolecule();new SmilesParser().parse(molecule,metadata.canonicalSmiles()==null||metadata.canonicalSmiles().isBlank()?metadata.smiles():metadata.canonicalSmiles());
        var provenance=new LinkedHashMap<String,Object>(extra);provenance.put("cacheReference",reference);provenance.put("sourceRow",metadata.sourceRow());provenance.put("sourceMoleculeId",metadata.moleculeId());provenance.put("sourceSnapshotId",assets.get(r.indexAssetId).sourceSnapshotId());
        provenance.put("scoreSemantics",r.pose()?"learned_pose_to_ensemble":"graph_only");provenance.put("modelAssetId",r.modelAssetId);
        return ResultService.makeHit(run.id(),r.indexAssetId+":"+reference.shardIndex()+":"+reference.localRow(),r.indexAssetId,null,null,List.of(),molecule,null,false,scores,provenance);
    }
    private static List<String> assetIds(Request r){List<String> ids=new ArrayList<>();if(r.indexAssetId!=null)ids.add(r.indexAssetId);ids.add(r.modelAssetId);if(r.compact())ids.add(r.compactAssetId);if(r.pose())ids.add(r.poseAssetId);return ids;}
    private static void close(List<AssetService.Lease> leases){for(var lease:leases)lease.close();}
    private static EmbeddingComparator guard(EmbeddingComparator delegate,JobService.Context ctx){return new EmbeddingComparator(){
        public float[][] compare(float[] q,float[][] c){ctx.check();var result=delegate.compare(q,c);ctx.check();return result;}
        public float[][] compareFlat(float[] q,float[] c,int count){ctx.check();var result=delegate.compareFlat(q,c,count);ctx.check();return result;}
    };}
    private static MoleculeFingerprintDataSource guard(MoleculeFingerprintDataSource delegate,JobService.Context ctx){return new MoleculeFingerprintDataSource(){
        public String artifactType(){return delegate.artifactType();}public long recordCount(){return delegate.recordCount();}public int shardCount(){return delegate.shardCount();}public MoleculeFingerprintProvenance provenance(){return delegate.provenance();}
        public void validateCompatibility(DeepSpaceModelBundle model,CompactSkelSpheresModelBundle compact,boolean required){delegate.validateCompatibility(model,compact,required);}
        public Map<MoleculeVectorReference,MoleculeFingerprintMetadata> resolveMetadata(Collection<MoleculeVectorReference> refs)throws IOException{ctx.check();var values=delegate.resolveMetadata(refs);ctx.check();return values;}
        public MoleculeFingerprintShardReader openShard(int shard,MoleculeFingerprintColumn column)throws IOException {ctx.check();var reader=delegate.openShard(shard,column);return new MoleculeFingerprintShardReader(){
            public MoleculeVectorBatch readBatch(int maximum)throws IOException{ctx.check();var batch=reader.readBatch(maximum);ctx.check();return batch;}
            public float[] readVector(long row)throws IOException{ctx.check();return reader.readVector(row);}public long recordCount(){return reader.recordCount();}public void close()throws IOException{reader.close();}
        };}
    };}
}
