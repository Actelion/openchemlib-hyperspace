package com.idorsia.research.chem.hyperspace3d.continuous;

import com.idorsia.research.chem.hyperspace3d.index.build.ProductFingerprintIndexBuildConfig;
import com.idorsia.research.chem.hyperspace3d.screening.ScreeningObjective;
import java.nio.file.Path;
import java.util.Set;

/** Versioned configuration for query-specific streaming combinatorial mining. */
public final class ContinuousDeepSpaceConfig {
    public int formatVersion=1;
    public Inputs inputs=new Inputs();
    public Query query=new Query();
    public Run run=new Run();
    public Output output=new Output();
    public Runtime runtime=new Runtime();
    public Exact exact=new Exact();
    public ProductFingerprintIndexBuildConfig.Filters filters=new ProductFingerprintIndexBuildConfig.Filters();
    public ProductFingerprintIndexBuildConfig.ReactionWeighting reactionWeighting=new ProductFingerprintIndexBuildConfig.ReactionWeighting();
    public void validate() {
        require(formatVersion==1,"unsupported formatVersion");
        require(text(inputs.rawFull)&&text(inputs.modelBundle)&&text(output.directory),"rawFull, modelBundle and output directory are required");
        require(text(query.identifier) && (text(query.smiles)^text(query.file)),"supply exactly one query smiles or MOL/SDF file, plus identifier");
        require(!pose() || text(inputs.poseBundle),"poseBundle is required for coordinate queries");
        require(query.record==null || (pose()&&query.record>=0),"record is zero based and requires a query file");
        require(run.maxAccepted>=0&&run.maxProposals>=0&&Double.isFinite(run.maxRuntimeSeconds)&&run.maxRuntimeSeconds>=0,"invalid run limits");
        require(output.globalTopK>0&&output.perReactionTopK>0&&output.checkpointSeconds>=0,"invalid elite/checkpoint limits");
        require(runtime.cpuWorkers>0&&runtime.encoderBatchSize>0&&runtime.encoderBatchSize<=65536
                &&runtime.queueCapacity>0&&runtime.intraOpThreads>=0&&runtime.cudaDeviceId>=0,"invalid runtime resources");
        require(Set.of("CPU","CUDA").contains(runtime.device),"device must be CPU or CUDA");
        require(filters.minHeavyAtoms>=6&&filters.maxHeavyAtoms<=32&&filters.maxHeavyAtoms>=filters.minHeavyAtoms&&filters.maxRotatableBonds>=0,"invalid molecular filters");
        require(Set.of("phesa_total","phesa_shape","phesa_pharmacophore").contains(run.target),"unsupported ranking target");
        require(Double.isFinite(run.ppWeight)&&run.ppWeight>=0&&run.ppWeight<=1,"ppWeight must be in [0,1]");
        require(exact.topK>0&&exact.maxConformers>0&&exact.threads>0,"invalid exact-PheSA resources/shortlist");
        require(Double.isFinite(exact.minScore)&&exact.minScore>=0&&exact.minScore<=1,"exact.minScore must be in [0,1]");
        reactionWeighting.validate();
    }
    public boolean pose(){return text(query.file);}
    public ScreeningObjective objective(){return ScreeningObjective.direct(run.target);}
    public static Path resolve(Path config,String value){Path p=Path.of(value);return (p.isAbsolute()?p:config.toAbsolutePath().getParent().resolve(p)).normalize();}
    private static boolean text(String s){return s!=null&&!s.isBlank();}
    private static void require(boolean b,String s){if(!b)throw new IllegalArgumentException(s);}
    public static final class Inputs {public String rawFull,rawDownsampled,modelBundle,poseBundle;}
    public static final class Query {public String identifier="query",smiles,file;public Integer record;}
    public static final class Run {public long seed=17,maxAccepted=0,maxProposals=0;public double maxRuntimeSeconds=3600,ppWeight=.5;public String target="phesa_total";}
    public static final class Output {public String directory;public boolean resume=true,persistFingerprints=false;public int globalTopK=10000,perReactionTopK=100,checkpointSeconds=30;}
    public static final class Runtime {public String device="CPU";public int cudaDeviceId=0,cpuWorkers=8,encoderBatchSize=512,queueCapacity=2,intraOpThreads=1;}
    public static final class Exact {public boolean enabled=false;public int topK=1000,maxConformers=8,threads=4;public double minScore=0;}
}
