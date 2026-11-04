package com.idorsia.research.chem.hyperspace.mcp;

import java.nio.file.*;
import java.util.*;

/** Public, chemistry-free adapter to existing detached import/build/fingerprint jobs. No workspace catalogue ownership. */
public final class DetachedPreparation {
    private final ServerConfig config;private final Jobs jobs;
    public DetachedPreparation(Path workerWorkspace,Map<String,Object> runtime)throws Exception {
        Path cli=Path.of(Objects.toString(runtime.get("cliJar"),"missing-cli.jar")).toAbsolutePath().normalize();
        config=new ServerConfig(workerWorkspace.toAbsolutePath(),cli,Objects.toString(runtime.getOrDefault("javaExecutable",Path.of(System.getProperty("java.home"),"bin","java"))),
                JsonFiles.integer(runtime,"threads",2,1,4096),ServerConfig.validHeap(Objects.toString(runtime.getOrDefault("heap","4G"))),"4G",
                path(runtime,"fingerprintJar"),path(runtime,"fingerprintLibDirectory"),path(runtime,"modelBundle"),path(runtime,"compactBundle"));
        jobs=new Jobs(config,System.getProperty("java.class.path"));
    }
    private static Path path(Map<String,Object> args,String key){return args.get(key)==null||args.get(key).toString().isBlank()?null:Path.of(args.get(key).toString()).toAbsolutePath().normalize();}
    public static Object inspect(String input,String zipEntry,boolean molecules)throws Exception{return molecules?FingerprintJobs.inspect(Path.of(input)):InputInspector.inspect(Path.of(input),zipEntry);}
    public Map<String,Object> plan(String kind,Map<String,Object> input)throws Exception {
        if(kind.equals("fingerprints"))return FingerprintJobs.normalize(input,config);
        if(!kind.equals("space"))throw new IllegalArgumentException("Preparation kind is space or fingerprints");
        if(Boolean.TRUE.equals(input.get("launchGui")))throw new IllegalArgumentException("Workbench preparation does not launch a desktop GUI");
        return Preparation.normalize(input,config);
    }
    public Map<String,Object> submit(String kind,Map<String,Object> input)throws Exception{return kind.equals("fingerprints")?jobs.submitFingerprint(input):jobs.submit(input);}
    public Map<String,Object> status(String id)throws Exception{return jobs.status(id);}
    public Map<String,Object> cancel(String id)throws Exception{return jobs.cancel(id);}
    public Object log(String id,String stream,long offset,int limit)throws Exception{return jobs.log(id,stream,offset,limit);}
    public static Object exportSlurm(Path profile,Path output)throws Exception{return SlurmFingerprintBundle.export(profile,output);}
}
