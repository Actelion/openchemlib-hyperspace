package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace.rawspace.*;
import com.idorsia.research.chem.hyperspace3d.continuous.*;
import com.idorsia.research.chem.hyperspace3d.cli.ContinuousDeepSpaceScreeningCLI;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.atomic.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class ContinuousDeepSpaceMinerTest {
    @TempDir Path temp;
    static RawSynthonSpace space()throws Exception {
        var builder=RawSynthonSpace.builder("toy-continuous");
        for(String id:List.of("r","s")) {
            builder.addRawFragments(id,0,List.of(synthon(id,0,"a","[U]CCC"),synthon(id,0,"b","[U]CCCC")));
            builder.addRawFragments(id,1,List.of(synthon(id,1,"x","[U]CCC"),synthon(id,1,"y","[U]CCO")));
        }
        return builder.build();
    }
    static RawSynthon synthon(String r,int position,String id,String smiles)throws Exception {
        var m=new StereoMolecule();new SmilesParser().parse(m,smiles);return RawSynthon.fromMolecule(r,position,id,m);
    }
    static ContinuousDeepSpaceConfig config() {
        var c=new ContinuousDeepSpaceConfig();c.inputs.rawFull="toy";c.inputs.modelBundle="model";c.output.directory="out";
        c.query.smiles="CCCCCC";c.runtime.encoderBatchSize=2;c.runtime.cpuWorkers=2;c.runtime.queueCapacity=2;
        c.output.globalTopK=3;c.output.perReactionTopK=2;c.output.checkpointSeconds=0;c.run.maxAccepted=8;c.run.maxRuntimeSeconds=0;c.exact.topK=2;
        return c;
    }
    static ProductEmbeddingBatchEncoder encoder(){return batch->{
        float[][] result=new float[batch.batchSize()][128];
        for(int i=0;i<result.length;i++)for(int a=0;a<32;a++)if(batch.atomMask()[i*32+a])result[i][0]+=.1f;
        return result;
    };}
    static EmbeddingComparator comparator(){return(q,candidates)->{
        float[][] scores=new float[candidates.length][3];for(int i=0;i<candidates.length;i++)Arrays.fill(scores[i],candidates[i][0]);return scores;
    };}
    ContinuousDeepSpaceMiner miner(ContinuousDeepSpaceConfig c,ProductEmbeddingBatchEncoder e,String identity)throws Exception {
        var model=new DeepSpaceModelManifest();model.availableTargets=List.of("phesa_total","phesa_shape","phesa_pharmacophore");
        return new ContinuousDeepSpaceMiner(space(),space(),e,comparator(),model,new float[128],c,identity,Map.of("rawspaceIdentity","toy",
                "rawspaceSha256","a".repeat(64),"downsampledRawspaceIdentity","toy","downsampledRawspaceSha256","a".repeat(64),"modelBundleHash","model"));
    }
    @Test void finiteRunProducesBoundedDeterministicElitesAndNoDefaultVectorCache()throws Exception {
        var c=config();Path out=temp.resolve("one");var a=miner(c,encoder(),"same").run(out);
        c.runtime.cpuWorkers=1;var b=miner(c,encoder(),"same").run(temp.resolve("two"));
        assertEquals(8,a.accepted);assertEquals("accepted_limit",a.stopReason);assertTrue(a.complete);
        assertEquals(3,a.global.size());assertTrue(a.byReaction.values().stream().allMatch(v->v.size()<=2));
        assertEquals(a.global.stream().map(h->h.tuple()).toList(),b.global.stream().map(h->h.tuple()).toList());
        assertFalse(Files.exists(out.resolve("fingerprints")));assertTrue(Files.isRegularFile(out.resolve("run.json")));
        assertEquals(8,miner(c,encoder(),"same").run(out).accepted);
        assertThrows(java.io.IOException.class,()->miner(c,encoder(),"changed").run(out));
    }
    @Test void recoveryReplaysUnconsumedPreparedBatchesAndRetainsCommittedShards()throws Exception {
        var c=config();c.output.persistFingerprints=true;c.run.maxAccepted=5;
        var calls=new AtomicInteger();ProductEmbeddingBatchEncoder broken=batch->{if(calls.incrementAndGet()==2)throw new IllegalStateException("injected failure");return encoder().encode(batch);};
        Path out=temp.resolve("resume");assertThrows(IllegalStateException.class,()->miner(c,broken,"same").run(out));
        var saved=ContinuousDeepSpaceMiner.JSON.readValue(out.resolve("checkpoint.json").toFile(),ContinuousDeepSpaceMiner.State.class);
        assertEquals(2,saved.accepted);
        var resumed=miner(c,encoder(),"same").run(out);var reference=miner(c,encoder(),"same").run(temp.resolve("reference"));
        assertEquals(5,resumed.accepted);assertEquals(reference.global.stream().map(h->h.tuple()).toList(),resumed.global.stream().map(h->h.tuple()).toList());
        var manifest=ContinuousDeepSpaceMiner.JSON.readValue(out.resolve("fingerprints.json").toFile(),ProductFingerprintIndexManifest.class);manifest.validate();
        Set<ProductTuple> tuples=new HashSet<>();
        for(var shard:manifest.shards)try(var in=new ProductFingerprintIndexReader(out.resolve(shard.path()))){for(var row:in.readBatch(100))assertTrue(tuples.add(row.tuple()));}
        assertEquals(5,tuples.size());
        Files.write(out.resolve(manifest.shards.get(0).path()),new byte[]{1});
        assertThrows(java.io.IOException.class,()->miner(c,encoder(),"same").run(out));
    }
    @Test void rejectedProductsExhaustSpaceWithoutInferenceAndProposalLimitIsExact()throws Exception {
        var c=config();c.filters.minHeavyAtoms=32;c.run.maxAccepted=0;
        var result=miner(c,b->{throw new AssertionError("all products should be rejected");},"reject").run(temp.resolve("reject"));
        assertEquals(0,result.accepted);assertEquals(8,result.proposals);assertEquals("space_exhausted",result.stopReason);
        c=config();c.run.maxProposals=3;c.run.maxAccepted=0;
        var limited=miner(c,encoder(),"limit").run(temp.resolve("limit"));assertEquals(3,limited.proposals);assertEquals("proposal_limit",limited.stopReason);
    }
    @Test void stopRequestsCheckpointAndResumeWhereasTimeLimitCompletes()throws Exception {
        var c=config();var ref=new AtomicReference<ContinuousDeepSpaceMiner>();
        var m=miner(c,b->{ref.get().requestStop();return encoder().encode(b);},"stop");ref.set(m);
        var stopped=m.run(temp.resolve("stop"));assertFalse(stopped.complete);assertEquals("interrupted",stopped.stopReason);
        assertEquals(8,miner(c,encoder(),"stop").run(temp.resolve("stop")).accepted);
        c=config();c.run.maxRuntimeSeconds=1e-9;
        var timed=miner(c,encoder(),"time").run(temp.resolve("time"));assertEquals("time_limit",timed.stopReason);assertTrue(timed.complete);
    }
    @Test void realOnnxCliAndExactFinishingPass()throws Exception {
        Path raw=temp.resolve("toy.rawspace.gz");RawSynthonSpaceIO.write(space(),raw);
        var c=config();c.inputs.rawFull=raw.toString();c.inputs.modelBundle=Path.of("model-bundles/deepspace7-v1").toAbsolutePath().toString();
        c.output.directory=temp.resolve("real").toString();c.output.globalTopK=8;c.run.maxAccepted=4;c.query.smiles="c1ccccc1";
        c.runtime.device=System.getProperty("continuous.device","CPU");
        c.exact.enabled=true;c.exact.topK=2;c.exact.threads=2;c.exact.maxConformers=1;c.exact.minScore=1;
        Path file=temp.resolve("config.json");ContinuousDeepSpaceMiner.JSON.writeValue(file.toFile(),c);
        ContinuousDeepSpaceScreeningCLI.main(new String[]{"--config",file.toString()});
        assertTrue(Files.isRegularFile(temp.resolve("real/hits.sdf")));
        var report=ContinuousDeepSpaceMiner.JSON.readTree(temp.resolve("real/exact-report.json").toFile());
        assertTrue(report.path("successful").asInt()>0);assertTrue(report.path("selected").asInt()<=2);
        assertEquals(0,report.path("passedFilter").asInt());assertEquals(0,Files.size(temp.resolve("real/exact-hits.sdf")));
        var rows=Files.readAllLines(temp.resolve("real/exact-hits.tsv")).stream().skip(1).map(row->row.split("\t")[5]).toList();
        assertEquals(rows.size(),new HashSet<>(rows).size());
        ContinuousDeepSpaceScreeningCLI.main(new String[]{"--config",file.toString()});
        assertEquals(4,ContinuousDeepSpaceMiner.JSON.readTree(temp.resolve("real/checkpoint.json").toFile()).path("accepted").asLong());
        assertEquals(report.path("elapsedNanos"),ContinuousDeepSpaceMiner.JSON.readTree(temp.resolve("real/exact-report.json").toFile()).path("elapsedNanos"));
    }
    @Test void realPoseCliPreservesQueryAndProducesDirectedScores()throws Exception {
        Assumptions.assumeTrue(System.getProperty("pose.deployment")!=null);
        Path raw=temp.resolve("pose.rawspace.gz");RawSynthonSpaceIO.write(space(),raw);
        var c=config();c.inputs.rawFull=raw.toString();c.inputs.modelBundle=Path.of("model-bundles/deepspace7-v1").toAbsolutePath().toString();
        c.inputs.poseBundle=Path.of(System.getProperty("pose.deployment"),"shared-candidate-seed17-verified").toString();
        c.query.smiles=null;c.query.file=Path.of("src/test/resources/pose-golden/query-0.mol").toAbsolutePath().toString();
        c.output.directory=temp.resolve("pose").toString();c.run.maxAccepted=2;
        c.runtime.device=System.getProperty("continuous.device","CPU");
        Path file=temp.resolve("pose-config.json");ContinuousDeepSpaceMiner.JSON.writeValue(file.toFile(),c);
        ContinuousDeepSpaceScreeningCLI.main(new String[]{file.toString()});
        var run=ContinuousDeepSpaceMiner.JSON.readTree(temp.resolve("pose/run.json").toFile());
        assertEquals("learned_pose_to_ensemble",run.path("scoreSemantics").asText());
        assertEquals("input_mol_axes",run.path("provenance").path("queryPose").path("coordinateFrame").asText());
        // A second JVM must compute the same identity despite randomized Map iteration order.
        var command=new ArrayList<String>(List.of(Path.of(System.getProperty("java.home"),"bin/java").toString(),"-cp",
                System.getProperty("java.class.path"),ContinuousDeepSpaceScreeningCLI.class.getName(),file.toString()));
        var process=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(temp.resolve("restart.log").toFile()).start();
        assertTrue(process.waitFor(30,java.util.concurrent.TimeUnit.SECONDS));assertEquals(0,process.exitValue(),Files.readString(temp.resolve("restart.log")));
    }
    @Test void refusesConcurrentWritersAndNonemptyUncheckpointedOutput()throws Exception {
        Path out=Files.createDirectories(temp.resolve("locked"));
        try(var channel=java.nio.channels.FileChannel.open(out.resolve(".miner.lock"),StandardOpenOption.CREATE,StandardOpenOption.WRITE);
            var lock=channel.lock()) {
            assertThrows(java.nio.channels.OverlappingFileLockException.class,()->miner(config(),encoder(),"lock").run(out));
        }
        Files.writeString(out.resolve("user.txt"),"user data");
        assertThrows(java.io.IOException.class,()->miner(config(),encoder(),"lock").run(out));
        assertEquals("user data",Files.readString(out.resolve("user.txt")));
    }
    @Test void fixedPoseExactDescriptorPreservesPlanarHeavyAtomDistances()throws Exception {
        var pose=com.idorsia.research.chem.hyperspace3d.feature.PoseQueryInput.read(Path.of("src/test/resources/pose-golden/query-2.mol"),null);
        var descriptor=ContinuousPheSAFinisher.fixedPoseDescriptor(pose.molecule());assertNotNull(descriptor);assertEquals(1,descriptor.getVolumes().size());
        var actual=descriptor.getConformer(0);actual.ensureHelperArrays(Molecule.cHelperCIP);
        for(int a=0;a<6;a++)for(int b=0;b<6;b++)assertEquals(pose.distances()[a*32+b],actual.getAtomCoordinates(a).distance(actual.getAtomCoordinates(b)),1e-5);
    }
}
