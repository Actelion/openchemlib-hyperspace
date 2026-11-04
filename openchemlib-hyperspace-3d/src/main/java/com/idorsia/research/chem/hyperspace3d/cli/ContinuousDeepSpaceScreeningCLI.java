package com.idorsia.research.chem.hyperspace3d.cli;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace.rawspace.*;
import com.idorsia.research.chem.hyperspace3d.continuous.*;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuilder;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.*;

/** Continuous random enumeration and batched learned PheSA screening. */
public final class ContinuousDeepSpaceScreeningCLI {
    public static void main(String[] args)throws Exception {
        if(args.length==1&&args[0].equals("--help")){System.out.println("Usage: ContinuousDeepSpaceScreeningCLI --config mining.json");return;}
        Path file=args.length==1?Path.of(args[0]):args.length==2&&args[0].equals("--config")?Path.of(args[1]):null;
        if(file==null)throw new IllegalArgumentException("Usage: ContinuousDeepSpaceScreeningCLI --config mining.json");
        var c=ContinuousDeepSpaceMiner.JSON.readValue(file.toFile(),ContinuousDeepSpaceConfig.class);c.validate();
        Path fullPath=ContinuousDeepSpaceConfig.resolve(file,c.inputs.rawFull);
        Path samplePath=c.inputs.rawDownsampled==null||c.inputs.rawDownsampled.isBlank()?fullPath:ContinuousDeepSpaceConfig.resolve(file,c.inputs.rawDownsampled);
        Path modelPath=ContinuousDeepSpaceConfig.resolve(file,c.inputs.modelBundle),output=ContinuousDeepSpaceConfig.resolve(file,c.output.directory);
        var full=RawSynthonSpaceIO.read(fullPath);var sample=samplePath.equals(fullPath)?full:RawSynthonSpaceIO.read(samplePath);
        var model=DeepSpaceModelBundle.load(modelPath);
        if(!c.pose()&&Math.abs(c.run.ppWeight-model.manifest().phesaPpWeight)>1e-12)
            throw new IllegalArgumentException("graph-mode ppWeight must match the model manifest's trained total head weight");
        var pose=c.pose()?SharedCandidatePoseBundle.load(ContinuousDeepSpaceConfig.resolve(file,c.inputs.poseBundle),model):null;
        var poseQuery=c.pose()?PoseQueryInput.read(ContinuousDeepSpaceConfig.resolve(file,c.query.file),c.query.record):null;
        StereoMolecule query;
        if(poseQuery!=null)query=poseQuery.molecule();else{query=new StereoMolecule();new SmilesParser().parse(query,c.query.smiles);}
        Map<String,Object> provenance=new LinkedHashMap<>();
        provenance.put("rawspaceIdentity",ProductFingerprintIndexBuilder.identity(full));provenance.put("rawspaceSha256",SharedCandidatePoseBundle.sha256(fullPath));
        provenance.put("downsampledRawspaceIdentity",ProductFingerprintIndexBuilder.identity(sample));provenance.put("downsampledRawspaceSha256",SharedCandidatePoseBundle.sha256(samplePath));
        provenance.put("modelBundleHash",model.bundleHash());provenance.put("queryIdcode",new Canonizer(query).getIDCode());
        provenance.put("modelContract",model.manifest());
        if(pose!=null){provenance.put("poseBundleHash",pose.hash());provenance.put("queryPose",poseQuery.provenance());}
        String identity=HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(
                ContinuousDeepSpaceMiner.JSON.writer().with(com.fasterxml.jackson.databind.SerializationFeature.ORDER_MAP_ENTRIES_BY_KEYS)
                        .writeValueAsBytes(List.of(provenance,c))));
        var runtime=new DeepSpaceOnnxEnvironment(DeepSpaceOnnxEnvironment.Device.valueOf(c.runtime.device),c.runtime.cudaDeviceId,c.runtime.intraOpThreads);
        try(var encoder=new DeepSpaceV1Encoder(runtime,model);
            var graph=pose==null?new DeepSpaceV1Comparator(runtime,model):null;
            var directed=pose==null?null:new SharedCandidatePoseRuntime(runtime,pose,c.run.ppWeight)) {
            long started=System.nanoTime();
            float[] embedding=pose==null?encoder.encode(new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(query)))[0]
                    :directed.encode(poseQuery.batch(),poseQuery.distances())[0];
            provenance.put("queryEncodingNanos",System.nanoTime()-started);
            var miner=new ContinuousDeepSpaceMiner(full,sample,encoder,pose==null?graph:directed,pose==null?model.manifest():pose.scoreManifest(),embedding,c,identity,provenance);
            Thread hook=new Thread(()->{miner.requestStop();try{miner.awaitStop();}catch(InterruptedException e){Thread.currentThread().interrupt();}},"continuous-miner-shutdown");
            java.lang.Runtime.getRuntime().addShutdownHook(hook);
            try {
                var result=miner.run(output,state->ContinuousPheSAFinisher.finish(output,miner,state,c,query));
                System.out.printf("Continuous DeepSpace stopped: %s; %,d scored / %,d proposals; results: %s%n",result.stopReason,result.accepted,result.proposals,output);
            } finally {try{java.lang.Runtime.getRuntime().removeShutdownHook(hook);}catch(IllegalStateException shuttingDown){/* JVM shutdown hook is already running. */}}
        }
    }
}
