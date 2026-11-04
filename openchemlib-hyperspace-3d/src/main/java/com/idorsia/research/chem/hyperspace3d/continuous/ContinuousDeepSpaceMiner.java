package com.idorsia.research.chem.hyperspace3d.continuous;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.TorsionDB;
import com.fasterxml.jackson.databind.*;
import com.idorsia.research.chem.hyperspace.SynthonAssembler;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductTupleSampler;
import com.idorsia.research.chem.hyperspace3d.model.*;
import com.idorsia.research.chem.hyperspace3d.screening.ScoredProductFingerprint;
import java.io.*;
import java.nio.channels.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;

/** Bounded CPU preparation -> one batch inference consumer, checkpointed at consumed batches. */
public final class ContinuousDeepSpaceMiner {
    public static final ObjectMapper JSON=new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);
    private static final Comparator<ScoredProductFingerprint> BEST=Comparator
            .comparingDouble(ScoredProductFingerprint::objective).reversed()
            .thenComparing(h->h.tuple().reactionId()).thenComparing(h->h.tuple().synthonIds().toString())
            .thenComparing(h->h.tuple().synthonOrdinals().toString());
    private final RawSynthonSpace full,sampling;
    private final ProductEmbeddingBatchEncoder encoder;
    private final EmbeddingComparator comparator;
    private final DeepSpaceModelManifest scores;
    private final float[] query;
    private final ContinuousDeepSpaceConfig config;
    private final String identity;
    private final Map<String,Object> provenance;
    private final AtomicBoolean stopped=new AtomicBoolean();
    private final CountDownLatch finished=new CountDownLatch(1);
    private final ThreadLocal<IDCodeParser> parsers=ThreadLocal.withInitial(IDCodeParser::new);
    private ProductFingerprintIndexWriter writer;
    private long shardRecords;

    public ContinuousDeepSpaceMiner(RawSynthonSpace full,RawSynthonSpace sampling,
            ProductEmbeddingBatchEncoder encoder,EmbeddingComparator comparator,DeepSpaceModelManifest scores,
            float[] query,ContinuousDeepSpaceConfig config,String identity,Map<String,Object> provenance) {
        config.validate();this.full=full;this.sampling=sampling;this.encoder=encoder;this.comparator=comparator;
        this.scores=scores;this.query=query.clone();this.config=config;this.identity=identity;this.provenance=provenance;
        if(query.length!=128)throw new IllegalArgumentException("query must be 128D");
        scores.targetIndex(config.run.target);
    }
    public void requestStop(){stopped.set(true);}
    public boolean isStopRequested(){return stopped.get();}
    public String runIdentity(){return identity;}
    public String scoreSemantics(){return config.pose()?"learned_pose_to_ensemble":"learned_graph_to_ensemble";}
    public float component(ScoredProductFingerprint hit,String target){return hit.components()[scores.targetIndex(target)];}
    public void awaitStop()throws InterruptedException{finished.await(30,TimeUnit.SECONDS);}

    public State run(Path directory)throws Exception {
        return run(directory,s->{});
    }
    @FunctionalInterface public interface Finisher {void finish(State state)throws Exception;}
    public State run(Path directory,Finisher finisher)throws Exception {
        Path dir=directory.toAbsolutePath().normalize();Files.createDirectories(dir);
        try(var channel=FileChannel.open(dir.resolve(".miner.lock"),StandardOpenOption.CREATE,StandardOpenOption.WRITE);
            var lock=channel.tryLock()) {
            if(lock==null)throw new IOException("continuous-miner output is locked");
            State state=runLocked(dir);if(state.complete&&!stopped.get())finisher.finish(state);return state;
        } finally {finished.countDown();}
    }
    private State runLocked(Path dir)throws Exception {
        var sampler=new ProductTupleSampler(full,sampling,config.run.seed,config.reactionWeighting.toWeighting());
        Path checkpoint=dir.resolve("checkpoint.json");
        State state;
        if(Files.exists(checkpoint)) {
            if(!config.output.resume)throw new IOException("checkpoint exists but resume is disabled");
            state=JSON.readValue(checkpoint.toFile(),State.class);
            if(state.formatVersion!=1||!identity.equals(state.identity))throw new IOException("checkpoint provenance/configuration differs");
            for(var shard:state.shards) {
                Path p=dir.resolve(shard.path()).normalize();
                if(!p.startsWith(dir)||!SharedCandidatePoseBundle.sha256(p).equals(shard.sha256()))throw new IOException("checkpoint fingerprint shard mismatch");
            }
            sampler.restore(state.sampler);
        } else {
            try(var files=Files.list(dir)) {
                if(files.anyMatch(p->!p.getFileName().toString().equals(".miner.lock")))throw new IOException("output directory is nonempty without a checkpoint");
            }
            state=new State();state.identity=identity;state.sampler=sampler.snapshot();
            atomicJson(checkpoint,state);
        }
        var global=new TreeSet<ScoredProductFingerprint>(BEST);global.addAll(state.global);
        var reactions=new TreeMap<String,TreeSet<ScoredProductFingerprint>>();
        state.byReaction.forEach((id,hits)->{var set=new TreeSet<ScoredProductFingerprint>(BEST);set.addAll(hits);reactions.put(id,set);});
        if(state.complete){publish(dir,state);return state;}
        state.stopReason="running";atomicJson(checkpoint,state);publish(dir,state);
        long started=System.nanoTime(),previousElapsed=state.elapsedNanos,lastCheckpoint=started;
        var queue=new ArrayBlockingQueue<Batch>(config.runtime.queueCapacity);
        var producer=Executors.newSingleThreadExecutor();
        long remainingProposals=config.run.maxProposals==0?Long.MAX_VALUE:Math.max(0,config.run.maxProposals-state.proposals);
        Future<?> producing=producer.submit(()->produce(sampler,queue,remainingProposals));
        try {
            while(true) {
                state.elapsedNanos=previousElapsed+System.nanoTime()-started;
                String reason=limit(state);
                if(reason!=null){state.stopReason=reason;state.complete=!reason.equals("interrupted");break;}
                Batch batch=queue.poll(200,TimeUnit.MILLISECONDS);
                if(batch==null) {
                    if(producing.isDone()){producing.get();state.stopReason="space_exhausted";state.complete=true;break;}
                    continue;
                }
                if(batch.failure!=null)throw new IllegalStateException("CPU preparation failed",batch.failure);
                state.proposals+=batch.proposals;state.assemblyNanos+=batch.assemblyNanos;state.featurizationNanos+=batch.featureNanos;
                batch.rejections.forEach((k,v)->state.rejections.merge(k,v,Long::sum));
                List<Prepared> accepted=batch.accepted;
                if(config.run.maxAccepted>0) {
                    int keep=(int)Math.min(accepted.size(),config.run.maxAccepted-state.accepted);
                    state.discardedAtAcceptedLimit+=accepted.size()-keep;accepted=accepted.subList(0,keep);
                }
                if(!accepted.isEmpty()) {
                    long stage=System.nanoTime();
                    var tensors=new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).buildFromFeatures(accepted.stream().map(Prepared::features).toList());
                    state.packingNanos+=System.nanoTime()-stage;stage=System.nanoTime();
                    float[][] embeddings=encoder.encode(tensors);state.encodingNanos+=System.nanoTime()-stage;stage=System.nanoTime();
                    if(embeddings.length!=accepted.size())throw new IllegalStateException("encoder returned wrong batch size");
                    float[][] values=comparator.compare(query,embeddings);state.comparisonNanos+=System.nanoTime()-stage;
                    if(values.length!=accepted.size())throw new IllegalStateException("comparator returned wrong batch size");
                    for(int i=0;i<values.length;i++) {
                        if(values[i].length!=scores.availableTargets.size())throw new IllegalStateException("comparator target width mismatch");
                        for(float v:values[i])if(!Float.isFinite(v))throw new IllegalStateException("nonfinite comparator output");
                        var hit=new ScoredProductFingerprint(accepted.get(i).tuple,config.objective().score(values[i],scores),values[i]);
                        retain(global,hit,config.output.globalTopK);
                        retain(reactions.computeIfAbsent(hit.tuple().reactionId(),id->new TreeSet<>(BEST)),hit,config.output.perReactionTopK);
                        if(config.output.persistFingerprints) {
                            if(writer==null){Files.createDirectories(dir.resolve("fingerprints"));writer=new ProductFingerprintIndexWriter(partial(dir,state.shards.size()));}
                            writer.write(new ProductFingerprintRecord(hit.tuple(),embeddings[i]));shardRecords++;
                        }
                        state.accepted++;
                    }
                }
                // The producer may be ahead. Only the cursor belonging to this consumed batch is durable.
                state.sampler=batch.snapshot;
                if(System.nanoTime()-lastCheckpoint>=config.output.checkpointSeconds*1_000_000_000L) {
                    state.elapsedNanos=previousElapsed+System.nanoTime()-started;
                    checkpoint(dir,state,global,reactions);lastCheckpoint=System.nanoTime();
                    System.out.printf("Continuous DeepSpace: %,d proposals, %,d scored, best %.5f%n",state.proposals,state.accepted,global.isEmpty()?0:global.first().objective());
                }
            }
            checkpoint(dir,state,global,reactions);
            return state;
        } finally {
            producing.cancel(true);producer.shutdownNow();
            if(writer!=null){writer.close();writer=null;shardRecords=0;}
        }
    }
    private String limit(State s) {
        if(stopped.get()||Thread.currentThread().isInterrupted())return "interrupted";
        if(config.run.maxAccepted>0&&s.accepted>=config.run.maxAccepted)return "accepted_limit";
        if(config.run.maxProposals>0&&s.proposals>=config.run.maxProposals)return "proposal_limit";
        if(config.run.maxRuntimeSeconds>0&&s.elapsedNanos/1e9>=config.run.maxRuntimeSeconds)return "time_limit";
        return null;
    }
    private void produce(ProductTupleSampler sampler,BlockingQueue<Batch> queue,long proposalBudget) {
        var workers=new ThreadPoolExecutor(config.runtime.cpuWorkers,config.runtime.cpuWorkers,0,TimeUnit.MILLISECONDS,
                new ArrayBlockingQueue<Runnable>(config.runtime.encoderBatchSize),new ThreadPoolExecutor.CallerRunsPolicy());
        try {
            long proposed=0;
            while(!stopped.get()&&!Thread.currentThread().isInterrupted()&&proposed<proposalBudget) {
                List<Callable<Prepared>> tasks=new ArrayList<>();
                for(int i=0;i<config.runtime.encoderBatchSize&&proposed<proposalBudget;i++) {
                    var sample=sampler.nextWeighted();if(sample==null)break;
                    proposed++;tasks.add(()->prepare(sample));
                }
                if(tasks.isEmpty())return;
                List<Prepared> accepted=new ArrayList<>();Map<String,Long> rejections=new TreeMap<>();long assembly=0,feature=0;
                for(var future:workers.invokeAll(tasks)) {
                    var p=future.get();assembly+=p.assemblyNanos;feature+=p.featureNanos;
                    if(p.reason==null)accepted.add(p);else rejections.merge(p.reason,1L,Long::sum);
                }
                queue.put(new Batch(accepted,tasks.size(),assembly,feature,rejections,sampler.snapshot(),null));
            }
        } catch(InterruptedException e){Thread.currentThread().interrupt();}
        catch(Exception e){try{queue.put(new Batch(List.of(),0,0,0,Map.of(),null,e));}catch(InterruptedException interrupted){Thread.currentThread().interrupt();}}
        finally{workers.shutdownNow();}
    }
    private Prepared prepare(ProductTupleSampler.Sample sample) {
        long started=System.nanoTime(),assembly=0;
        try {
            List<StereoMolecule> parts=new ArrayList<>();
            for(var synthon:sample.synthons()){StereoMolecule m=new StereoMolecule();parsers.get().parse(m,synthon.getIdcode());m.ensureHelperArrays(Molecule.cHelperCIP);parts.add(m);}
            StereoMolecule molecule=SynthonAssembler.assembleSynthons_faster(parts);molecule.ensureHelperArrays(Molecule.cHelperCIP);
            assembly=System.nanoTime()-started;started=System.nanoTime();
            String invalid=null;
            if(molecule.getAtoms()<config.filters.minHeavyAtoms)invalid="TOO_FEW_HEAVY_ATOMS";
            if(molecule.getAtoms()>config.filters.maxHeavyAtoms)invalid="TOO_MANY_HEAVY_ATOMS";
            if(config.filters.maxRotatableBonds>0) {
                boolean[] rotatable=new boolean[molecule.getBonds()];TorsionDB.findRotatableBonds(molecule,true,rotatable);int count=0;for(boolean value:rotatable)if(value)count++;
                if(count>config.filters.maxRotatableBonds)invalid="TOO_MANY_ROTATABLE_BONDS";
            }
            FeaturizationResult result=invalid==null?new OCLDeepSpaceFeaturizer().featurize(molecule):null;
            return new Prepared(sample.tuple(),result,invalid!=null?invalid:result.rejectionReason(),assembly,System.nanoTime()-started);
        } catch(RuntimeException e){return new Prepared(sample.tuple(),null,"ASSEMBLY_OR_FEATURE_FAILED:"+e.getClass().getSimpleName(),assembly==0?System.nanoTime()-started:assembly,0);}
    }
    public StereoMolecule assemble(ProductTuple tuple) {
        var reaction=full.getReactions().get(tuple.reactionId());if(reaction==null)throw new IllegalArgumentException("unknown reaction");
        var positions=new ArrayList<>(reaction.getRawFragmentSets().keySet());Collections.sort(positions);
        if(positions.size()!=tuple.synthonIds().size())throw new IllegalArgumentException("tuple position mismatch");
        List<StereoMolecule> parts=new ArrayList<>();
        for(int i=0;i<positions.size();i++) {
            var synthon=reaction.getRawFragmentSets().get(positions.get(i)).get(tuple.synthonOrdinals().get(i));
            if(!synthon.getFragmentId().equals(tuple.synthonIds().get(i)))throw new IllegalArgumentException("tuple ID/ordinal mismatch");
            var m=new StereoMolecule();parsers.get().parse(m,synthon.getIdcode());m.ensureHelperArrays(Molecule.cHelperCIP);parts.add(m);
        }
        return SynthonAssembler.assembleSynthons_faster(parts);
    }
    private static void retain(TreeSet<ScoredProductFingerprint> set,ScoredProductFingerprint hit,int capacity){set.add(hit);if(set.size()>capacity)set.pollLast();}
    private static Path partial(Path dir,int shard){return dir.resolve(String.format("fingerprints/part-%06d.bin.partial",shard));}
    private void checkpoint(Path dir,State s,TreeSet<ScoredProductFingerprint> global,Map<String,TreeSet<ScoredProductFingerprint>> reactions)throws Exception {
        if(writer!=null) {
            writer.close();writer=null;
            Path p=partial(dir,s.shards.size()),target=p.resolveSibling(p.getFileName().toString().replace(".partial",""));
            Files.move(p,target,StandardCopyOption.REPLACE_EXISTING,StandardCopyOption.ATOMIC_MOVE);
            s.shards.add(new ProductFingerprintIndexShard(dir.relativize(target).toString(),shardRecords,SharedCandidatePoseBundle.sha256(target)));shardRecords=0;
        }
        s.global=new ArrayList<>(global);s.byReaction=new TreeMap<>();reactions.forEach((id,hits)->s.byReaction.put(id,new ArrayList<>(hits)));
        atomicJson(dir.resolve("checkpoint.json"),s);publish(dir,s);
    }
    private void publish(Path dir,State s)throws Exception {
        Map<String,Object> run=new LinkedHashMap<>();run.put("artifactType","hyperspace3d-continuous-learned-screen");run.put("formatVersion",1);
        run.put("identity",identity);run.put("provenance",provenance);run.put("configuration",config);run.put("targetNames",scores.availableTargets);
        run.put("scoreSemantics",scoreSemantics());
        run.put("samplingAlgorithm",ProductTupleSampler.ALGORITHM);run.put("state",s);
        atomicJson(dir.resolve("run.json"),run);
        writeHits(dir.resolve("hits.tsv"),s.global,"global");
        writeHits(dir.resolve("reaction-hits.tsv"),s.byReaction.values().stream().flatMap(List::stream).toList(),"reaction_elite");
        if(config.output.persistFingerprints) {
            var m=new ProductFingerprintIndexManifest();
            m.rawspaceIdentity=Objects.toString(provenance.get("rawspaceIdentity"));m.rawspaceSha256=Objects.toString(provenance.get("rawspaceSha256"));
            m.downsampledRawspaceIdentity=Objects.toString(provenance.get("downsampledRawspaceIdentity"));m.downsampledRawspaceSha256=Objects.toString(provenance.get("downsampledRawspaceSha256"));
            m.modelBundleHash=Objects.toString(provenance.get("modelBundleHash"));m.samplingAlgorithm=ProductTupleSampler.ALGORITHM;m.samplingSeed=config.run.seed;
            m.reactionWeightConfiguration=JSON.convertValue(config.reactionWeighting,new com.fasterxml.jackson.core.type.TypeReference<Map<String,Object>>(){});
            m.molecularFilters=JSON.convertValue(config.filters,new com.fasterxml.jackson.core.type.TypeReference<Map<String,Object>>(){});
            m.recordCount=s.accepted;m.shards=s.shards;m.buildConfigurationHash=identity;m.validate();atomicJson(dir.resolve("fingerprints.json"),m);
        }
        Path p=dir.resolve("summary.md.partial");
        Files.writeString(p,"# Continuous DeepSpace screening\n\nScored products: "+s.accepted+"; proposals: "+s.proposals+"; state: "+s.stopReason+".\n\n"
                +"Scores are learned "+(config.pose()?"pose→ensemble":"graph→ensemble")+" predictions. Exact PheSA, when enabled, is reported separately in exact-hits.tsv.\n\n"
                +"Elapsed: "+s.elapsedNanos/1e9+" s; summed CPU assembly: "+s.assemblyNanos/1e9+" s; summed CPU features: "+s.featurizationNanos/1e9
                +" s; tensor packing: "+s.packingNanos/1e9+" s; encoding: "+s.encodingNanos/1e9+" s; comparison: "+s.comparisonNanos/1e9+" s.\n");
        Files.move(p,dir.resolve("summary.md"),StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
    }
    private void writeHits(Path path,List<ScoredProductFingerprint> hits,String source)throws Exception {
        Path p=path.resolveSibling(path.getFileName()+".partial");
        try(var out=Files.newBufferedWriter(p)) {
            out.write("format_version\trank\tsource\treaction_id\tsynthon_ids_json\tsynthon_ordinals_json\tobjective\tpredicted_total\tpredicted_shape\tpredicted_pp\n");int rank=0;
            for(var h:hits)out.write("1\t"+(++rank)+"\t"+source+"\t"+clean(h.tuple().reactionId())+"\t"+compact(h.tuple().synthonIds())+"\t"+compact(h.tuple().synthonOrdinals())+"\t"+h.objective()+"\t"
                    +h.components()[scores.targetIndex("phesa_total")]+"\t"+h.components()[scores.targetIndex("phesa_shape")]+"\t"+h.components()[scores.targetIndex("phesa_pharmacophore")]+"\n");
        }
        Files.move(p,path,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
    }
    public static String compact(Object value)throws Exception{return JSON.writer().without(SerializationFeature.INDENT_OUTPUT).writeValueAsString(value);}
    public static String clean(String s){return s.replace('\t',' ').replace('\n',' ').replace('\r',' ');}
    public static void atomicJson(Path path,Object value)throws Exception {
        Path p=path.resolveSibling(path.getFileName()+".partial");JSON.writeValue(p.toFile(),value);
        Files.move(p,path,StandardCopyOption.REPLACE_EXISTING,StandardCopyOption.ATOMIC_MOVE);
    }
    private record Prepared(ProductTuple tuple,FeaturizationResult features,String reason,long assemblyNanos,long featureNanos){}
    private record Batch(List<Prepared> accepted,long proposals,long assemblyNanos,long featureNanos,Map<String,Long> rejections,ProductTupleSampler.Snapshot snapshot,Throwable failure){}
    public static final class State {
        public int formatVersion=1;public String identity,stopReason="running";public boolean complete;
        public long proposals,accepted,discardedAtAcceptedLimit,elapsedNanos,assemblyNanos,featurizationNanos,packingNanos,encodingNanos,comparisonNanos;
        public ProductTupleSampler.Snapshot sampler;
        public Map<String,Long> rejections=new TreeMap<>();
        public List<ScoredProductFingerprint> global=new ArrayList<>();
        public Map<String,List<ScoredProductFingerprint>> byReaction=new TreeMap<>();
        public List<ProductFingerprintIndexShard> shards=new ArrayList<>();
    }
}
