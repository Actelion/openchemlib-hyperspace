package com.idorsia.research.chem.hyperspace3d.continuous;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.*;
import com.actelion.research.chem.phesa.*;
import com.idorsia.research.chem.hyperspace3d.mining.*;
import com.idorsia.research.chem.hyperspace3d.screening.ScoredProductFingerprint;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;

/** Optional final exact-label pass; never feeds back into the sampler or neural scorer. */
public final class ContinuousPheSAFinisher {
    private record Result(ScoredProductFingerprint hit,String idcode,String smiles,ExactPheSALabel exact){}
    public static void finish(Path output,ContinuousDeepSpaceMiner miner,ContinuousDeepSpaceMiner.State state,
            ContinuousDeepSpaceConfig config,StereoMolecule query)throws Exception {
        long started=System.nanoTime();
        writeSdf(output.resolve("hits.sdf"),state.global,miner,null);
        if(!config.exact.enabled)return;
        Path completedReport=output.resolve("exact-report.json");
        if(Files.isRegularFile(completedReport)&&Files.isRegularFile(output.resolve("exact-hits.tsv"))&&Files.isRegularFile(output.resolve("exact-hits.sdf"))) {
            var previous=ContinuousDeepSpaceMiner.JSON.readTree(completedReport.toFile());
            if(previous.path("complete").asBoolean()&&state.identity.equals(previous.path("runIdentity").asText()))return;
        }
        var base=new ExactPheSALabeler(config.exact.maxConformers,config.run.ppWeight);
        PheSAMolecule queryDescriptor=config.pose()?fixedPoseDescriptor(query):base.describe(query);
        String encodedQuery=base.encode(queryDescriptor);
        if(encodedQuery==null)throw new IllegalStateException("exact-PheSA query descriptor failed; learned results remain available");
        var handlers=ThreadLocal.withInitial(()->new ExactPheSALabeler(config.exact.maxConformers,config.run.ppWeight));
        var queryDescriptors=ThreadLocal.withInitial(()->handlers.get().decode(encodedQuery));
        Set<String> seen=new HashSet<>();List<Callable<Result>> tasks=new ArrayList<>();
        for(var hit:state.global) {
            if(tasks.size()>=config.exact.topK||miner.isStopRequested())break;
            StereoMolecule molecule=miner.assemble(hit.tuple());String idcode=new Canonizer(molecule).getIDCode();
            if(!seen.add(idcode))continue;
            String smiles=new IsomericSmilesCreator(molecule).getSmiles();
            tasks.add(()-> {
                if(miner.isStopRequested())throw new InterruptedException("exact finishing interrupted");
                var handler=handlers.get();return new Result(hit,idcode,smiles,
                        handler.score(queryDescriptors.get(),handler.describe(molecule)));
            });
        }
        var workers=Executors.newFixedThreadPool(config.exact.threads);
        List<Result> results=new ArrayList<>();
        try {
            for(var future:workers.invokeAll(tasks))results.add(future.get());
        } finally{workers.shutdownNow();}
        results.sort(Comparator.<Result>comparingDouble(r->r.exact.status()==0?exactObjective(r.exact,config.run.target):Double.NEGATIVE_INFINITY)
                .reversed().thenComparing(Result::idcode));
        Path path=output.resolve("exact-hits.tsv"),partial=output.resolve("exact-hits.tsv.partial");
        try(var out=Files.newBufferedWriter(partial)) {
            out.write("format_version\trank\treaction_id\tsynthon_ids_json\tsynthon_ordinals_json\tcanonical_smiles\tlearned_objective\texact_status\texact_total\texact_shape\texact_pp\tpasses_filter\tpredicted_total\tpredicted_shape\tpredicted_pp\n");int rank=0;
            for(var r:results)out.write("1\t"+(++rank)+"\t"+ContinuousDeepSpaceMiner.clean(r.hit.tuple().reactionId())+"\t"
                    +ContinuousDeepSpaceMiner.compact(r.hit.tuple().synthonIds())+"\t"+ContinuousDeepSpaceMiner.compact(r.hit.tuple().synthonOrdinals())
                    +"\t"+r.smiles+"\t"+r.hit.objective()+"\t"+r.exact.status()+"\t"+number(r.exact.total())+"\t"+number(r.exact.shape())+"\t"+number(r.exact.pharmacophore())+"\t"+passes(r.exact,config)
                    +"\t"+miner.component(r.hit,"phesa_total")+"\t"+miner.component(r.hit,"phesa_shape")+"\t"+miner.component(r.hit,"phesa_pharmacophore")+"\n");
        }
        Files.move(partial,path,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
        var exactByTuple=new HashMap<com.idorsia.research.chem.hyperspace3d.index.ProductTuple,ExactPheSALabel>();
        for(var r:results)exactByTuple.put(r.hit.tuple(),r.exact);
        writeSdf(output.resolve("exact-hits.sdf"),results.stream().filter(r->passes(r.exact,config)).map(Result::hit).toList(),miner,exactByTuple);
        ContinuousDeepSpaceMiner.atomicJson(output.resolve("exact-report.json"),Map.of("artifactType","hyperspace3d-continuous-exact-phesa",
                "formatVersion",1,"runIdentity",state.identity,"queryConvention",config.pose()?"fixed_supplied_pose":"generated_query_ensemble",
                "ppWeight",config.run.ppWeight,"maxConformers",config.exact.maxConformers,"selected",results.size(),
                "successful",results.stream().filter(r->r.exact.status()==0).count(),"elapsedNanos",System.nanoTime()-started));
        var report=ContinuousDeepSpaceMiner.JSON.readTree(output.resolve("exact-report.json").toFile());
        ((com.fasterxml.jackson.databind.node.ObjectNode)report).put("minScore",config.exact.minScore)
                .put("target",config.run.target).put("passedFilter",results.stream().filter(r->passes(r.exact,config)).count()).put("complete",true);
        if(miner.isStopRequested())throw new InterruptedException("exact finishing interrupted; rerun to complete exports");
        ContinuousDeepSpaceMiner.atomicJson(output.resolve("exact-report.json"),report);
    }
    public static PheSAMolecule fixedPoseDescriptor(StereoMolecule supplied) {
        StereoMolecule copy=new StereoMolecule(supplied);copy.ensureHelperArrays(Molecule.cHelperCIP);int n=copy.getAtoms();
        double[][] coordinates=new double[n][3];for(int a=0;a<n;a++)coordinates[a]=new double[]{copy.getAtomX(a),copy.getAtomY(a),copy.getAtomZ(a)};
        new HydrogenAssembler(copy).addImplicitHydrogens();
        for(int a=0;a<n;a++)if(copy.getAtomX(a)!=coordinates[a][0]||copy.getAtomY(a)!=coordinates[a][1]||copy.getAtomZ(a)!=coordinates[a][2])
            throw new IllegalStateException("hydrogen preparation changed supplied query pose");
        var set=new ConformerSet();set.add(new Conformer(copy));
        // The ConformerSet overload preserves even planar supplied poses.
        return new DescriptorHandlerShape().createDescriptor(set);
    }
    private static double exactObjective(ExactPheSALabel score,String target){return switch(target){case "phesa_shape"->score.shape();case "phesa_pharmacophore"->score.pharmacophore();default->score.total();};}
    private static boolean passes(ExactPheSALabel score,ContinuousDeepSpaceConfig c){return score.status()==0&&exactObjective(score,c.run.target)>=c.exact.minScore;}
    private static String number(float v){return Float.isFinite(v)?Float.toString(v):"";}
    private static void writeSdf(Path path,List<ScoredProductFingerprint> hits,ContinuousDeepSpaceMiner miner,
            Map<com.idorsia.research.chem.hyperspace3d.index.ProductTuple,ExactPheSALabel> exact)throws Exception {
        Path partial=path.resolveSibling(path.getFileName()+".partial");
        try(var out=Files.newBufferedWriter(partial)) {
            int rank=0;for(var hit:hits) {
                if(miner.isStopRequested())break;
                StereoMolecule molecule=miner.assemble(hit.tuple());new com.actelion.research.chem.coords.CoordinateInventor().invent(molecule);
                molecule.setName(hit.tuple().reactionId());new MolfileCreator(molecule).writeMolfile(out);
                property(out,"RANK",Integer.toString(++rank));property(out,"REACTION_ID",hit.tuple().reactionId());
                property(out,"FORMAT_VERSION","1");property(out,"RUN_IDENTITY",miner.runIdentity());property(out,"SCORE_SEMANTICS",miner.scoreSemantics());
                property(out,"SYNTHON_IDS",ContinuousDeepSpaceMiner.compact(hit.tuple().synthonIds()));
                property(out,"SYNTHON_ORDINALS",ContinuousDeepSpaceMiner.compact(hit.tuple().synthonOrdinals()));
                property(out,"LEARNED_OBJECTIVE",Double.toString(hit.objective()));
                property(out,"PREDICTED_TOTAL",Float.toString(miner.component(hit,"phesa_total")));
                property(out,"PREDICTED_SHAPE",Float.toString(miner.component(hit,"phesa_shape")));
                property(out,"PREDICTED_PP",Float.toString(miner.component(hit,"phesa_pharmacophore")));
                if(exact!=null){var score=exact.get(hit.tuple());property(out,"EXACT_TOTAL",number(score.total()));property(out,"EXACT_SHAPE",number(score.shape()));property(out,"EXACT_PP",number(score.pharmacophore()));}
                property(out,"STRUCTURE_SEMANTICS","2D product depiction; no aligned conformer exported");out.write("$$$$\n");
            }
        }
        Files.move(partial,path,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
    }
    private static void property(java.io.Writer out,String name,String value)throws Exception{out.write("> <"+name+">\n"+ContinuousDeepSpaceMiner.clean(value)+"\n\n");}
}
