package com.idorsia.research.chem.hyperspace.mcp;

import java.io.File;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/** Typed pose operations, launched through the configured native-runtime distribution. */
final class PoseSearchJobs {
    static List<String> command(ServerConfig c,String main,List<String> args,String heap) {
        if(c.fingerprintJar()==null || c.fingerprintLibDirectory()==null)
            throw new IllegalArgumentException("Configure fingerprintJar and fingerprintLibDirectory for the 3D runtime");
        var command=new ArrayList<>(List.of(c.javaExecutable(),"-Xmx"+ServerConfig.validHeap(heap),"-cp",
                c.fingerprintJar()+File.pathSeparator+c.fingerprintLibDirectory().resolve("*"),main));
        command.addAll(args);return command;
    }
    static Map<String,Object> inspect(ServerConfig c,Map<String,Object> args)throws Exception {
        var parameters=new ArrayList<>(List.of(Path.of(JsonFiles.string(args,"queryFile")).toRealPath().toString()));
        if(args.containsKey("record"))parameters.add(""+JsonFiles.integer(args,"record",0,0,Integer.MAX_VALUE));
        Path log=Files.createTempFile(c.workspace(),"pose-inspect-",".json");
        Process p=null;
        try {
            p=new ProcessBuilder(command(c,"com.idorsia.research.chem.hyperspace3d.cli.InspectPoseQueryCLI",parameters,"512M"))
                    .redirectError(ProcessBuilder.Redirect.INHERIT).redirectOutput(log.toFile()).start();
            if(!p.waitFor(30,TimeUnit.SECONDS))throw new IllegalStateException("query inspection timed out");
            if(p.exitValue()!=0)throw new IllegalArgumentException("query inspection failed; verify MOL/SDF coordinates and record selection");
            if(Files.size(log)>65536)throw new IllegalStateException("query summary exceeds 64 KiB");
            return JsonFiles.read(log);
        } finally {if(p!=null && p.isAlive())Jobs.terminate(p.toHandle());Files.deleteIfExists(log);}
    }
    static Map<String,Object> normalize(ServerConfig c,Map<String,Object> args)throws Exception {
        var result=new LinkedHashMap<String,Object>();
        result.put("jobKind","pose-search");
        for(String name:List.of("queryFile","index","modelBundle","poseBundle")) {
            Path path=Path.of(JsonFiles.string(args,name)).toRealPath();result.put(name,path.toString());
        }
        if(args.containsKey("record"))result.put("record",JsonFiles.integer(args,"record",0,0,Integer.MAX_VALUE));
        result.put("identifier",args.getOrDefault("identifier","pose-query"));
        String device=Objects.toString(args.getOrDefault("device","CPU")).toUpperCase(Locale.ROOT);
        if(!Set.of("CPU","CUDA").contains(device))throw new IllegalArgumentException("device must be CPU or CUDA");
        result.put("device",device);
        String target=Objects.toString(args.getOrDefault("objective","phesa_total"));
        if(!Set.of("phesa_total","phesa_shape","phesa_pharmacophore").contains(target))throw new IllegalArgumentException("unsupported pose objective");
        result.put("objective",target);
        double weight=((Number)args.getOrDefault("ppWeight",.5)).doubleValue();
        if(!Double.isFinite(weight)||weight<0||weight>1)throw new IllegalArgumentException("ppWeight must be in [0,1]");
        result.put("ppWeight",weight);
        result.put("topK",JsonFiles.integer(args,"topK",1000,1,1_000_000));
        result.put("batchSize",JsonFiles.integer(args,"batchSize",32768,1,1_000_000));
        result.put("prefetchDepth",JsonFiles.integer(args,"prefetchDepth",2,1,32));
        result.put("cudaDeviceId",JsonFiles.integer(args,"cudaDeviceId",0,0,1024));
        result.put("intraOpThreads",JsonFiles.integer(args,"intraOpThreads",1,1,4096));
        result.put("heap",ServerConfig.validHeap(Objects.toString(args.getOrDefault("heap",c.heap()))));
        return result;
    }
    static Map<String,Object> config(Map<String,Object> r,Path output) {
        var query=new LinkedHashMap<String,Object>();query.put("file",r.get("queryFile"));
        query.put("identifier",r.get("identifier"));if(r.containsKey("record"))query.put("record",r.get("record"));
        return Map.of("formatVersion",1,"inputs",Map.of("index",r.get("index"),"modelBundle",r.get("modelBundle"),"poseBundle",r.get("poseBundle")),
                "query",query,"screening",Map.of("mode","pose_phesa_3d","posePpWeight",r.get("ppWeight"),
                        "resultTopK",r.get("topK"),"cascadeShortlistTopK",r.get("topK"),"objective",Map.of("type","direct","target",r.get("objective"))),
                "runtime",Map.of("device",r.get("device"),"cudaDeviceId",r.get("cudaDeviceId"),"comparatorBatchSize",r.get("batchSize"),"prefetchDepth",r.get("prefetchDepth"),"intraOpThreads",r.get("intraOpThreads")),
                "output",Map.of("hitsTsv",output.resolve("hits.tsv").toString(),"hitsSdf",output.resolve("hits.sdf").toString(),
                        "summaryMarkdown",output.resolve("summary.md").toString(),"runManifest",output.resolve("run.json").toString(),
                        "reportTopK",Math.min(25,((Number)r.get("topK")).intValue())));
    }
}
