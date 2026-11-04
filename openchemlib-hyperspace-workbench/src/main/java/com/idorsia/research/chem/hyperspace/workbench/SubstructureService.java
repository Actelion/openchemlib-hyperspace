package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.fasterxml.jackson.core.type.TypeReference;
import com.idorsia.research.chem.hyperspace.*;
import java.io.*;
import java.math.BigInteger;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.atomic.*;
import java.util.zip.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Headless real combinatorial search and explicitly bounded, deterministic expansion. */
public final class SubstructureService {
    public record Options(int threads,int familyLimit,boolean alignedDepictions) {
        public void validate(){require(threads>=1&&threads<=64,"Search threads must be 1–64");require(familyLimit>=1&&familyLimit<=1000,"Family limit must be 1–1000");}
    }
    public record Position(int position,int count,List<ResultService.Synthon> sample){}
    public record Family(String id,String runId,String assetId,String reaction,String representedCombinations,List<Position> positions){}
    public record FamilyData(Family summary,List<List<ResultService.Synthon>> prefix){}
    private final WorkspaceStore store;private final AssetService assets;private final QueryService queries;private final JobService jobs;private final ResultService results;
    public SubstructureService(WorkspaceStore store,AssetService assets,QueryService queries,JobService jobs,ResultService results){this.store=store;this.assets=assets;this.queries=queries;this.jobs=jobs;this.results=results;}
    public JobService.Job search(long expected,String actor,String queryId,List<String> assetIds,Options options)throws Exception {
        options.validate();require(assetIds!=null&&!assetIds.isEmpty()&&assetIds.size()<=8,"Select 1–8 combinatorial indexes");
        require(assetIds.stream().distinct().count()==assetIds.size(),"Duplicate search assets");var query=queries.revision(queryId);
        var molecule=QueryChemistry.parse(query.document());molecule.setFragment(true);
        require(molecule.getAtoms()<=64,"Combinatorial query is limited to 64 atoms in this prototype");
        SubstructureSearchHelper.expandBridgedSearches(molecule); // Validate bridge limit before admitting a job.
        List<AssetService.Lease> leases=new ArrayList<>();
        try {
            for(String id:assetIds) {
                require(assets.get(id).kind()==AssetService.Kind.SUBSTRUCTURE_INDEX,"Select a substructure index; rawspaces must be prepared first");
                var lease=assets.acquire(id);leases.add(lease);require(lease.value() instanceof SynthonSpace,"Asset is not a SynthonSpace");
            }
            synchronized(store) {
                var run=results.create(expected,actor,"substructure",query,assetIds,options,
                        "Combinatorial family search with the engine's per-pattern candidate guard; at most "+options.familyLimit()+" retained families. Family counts can overlap and are not a deduplicated space total."
                        +(QueryChemistry.inspect(queryId,query.document()).bonds().stream().anyMatch(QueryChemistry.Bond::bridge)?" Bridge families are expanded-path candidates; concrete products are validated with OCL's original bridge query after assembly.":""));
                return jobs.submit(store.snapshot().revision,actor,"substructure",run.id(),Map.of("runId",run.id(),"queryId",queryId,"assets",assetIds,"options",options),leases,
                        ctx->executeSearch(ctx,run,leases,options));
            }
        }catch(Exception e){for(var lease:leases)lease.close();throw e;}
    }
    private Object executeSearch(JobService.Context ctx,ResultService.Run run,List<AssetService.Lease> leases,Options options)throws Exception {
        results.state(run.id(),"running",Map.of());List<Family> summaries=new ArrayList<>();Set<String> signatures=new HashSet<>();
        AtomicInteger discovered=new AtomicInteger();AtomicBoolean truncated=new AtomicBoolean();AtomicReference<Exception> failure=new AtomicReference<>();
        Path directory=results.directory(run.id());Files.createDirectories(directory.resolve("families"));
        for(int assetIndex=0;assetIndex<leases.size();assetIndex++) {
            ctx.check();String asset=run.assetIds().get(assetIndex);SynthonSpace space=(SynthonSpace)leases.get(assetIndex).value();
            ctx.stage("Searching "+assets.get(asset).name());
            var molecule=QueryChemistry.parse(run.query().document());molecule.setFragment(true);
            SubstructureSearchHelper.run_substructure_search_streaming_01_withBridgedBondsExpansion(space,new CachedDescriptorProvider(space.getDescriptorHandlerShortName()),molecule,
                    options.threads(),true,false,new SynthonSpace.CombinatorialHitReceiver(){
                        @Override public synchronized void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hits,List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> statistics) {
                            if(ctx.cancelled()||failure.get()!=null)return;
                            for(var hit:hits)try {
                                discovered.incrementAndGet();if(summaries.size()>=options.familyLimit()){truncated.set(true);continue;}
                                FamilyData data=family(run.id(),asset,hit);if(!signatures.add(data.summary().id()))continue;
                                try(var out=new GZIPOutputStream(Files.newOutputStream(directory.resolve("families").resolve(data.summary().id()+".json.gz")))){WorkspaceStore.JSON.writeValue(out,data);}
                                summaries.add(data.summary());
                            }catch(Exception e){failure.compareAndSet(null,e);return;}
                        }
                    });
            ctx.check();if(failure.get()!=null)throw failure.get();
        }
        WorkspaceStore.JSON.writeValue(directory.resolve("families.json").toFile(),summaries);
        String represented=summaries.stream().map(f->new BigInteger(f.representedCombinations())).reduce(BigInteger.ZERO,BigInteger::add).toString();
        var diagnostics=Map.of("retainedFamilies",summaries.size(),"discoveredFamilies",discovered.get(),"familyLimitReached",truncated.get(),
                "summedRepresentedCombinations",represented,"countsMayOverlap",true,"productsExpanded",0,"engineCandidateGuard",1000);
        results.complete(run.id(),List.of(),diagnostics);ctx.log("Retained "+summaries.size()+" families; expansion is an explicit bounded operation.");return diagnostics;
    }
    private FamilyData family(String run,String asset,SynthonSpace.CombinatorialHit hit) {
        var entries=hit.hit_fragments.entrySet().stream().sorted(Comparator.comparingInt(e->e.getKey().frag)).toList();
        List<Position> positions=new ArrayList<>();List<List<ResultService.Synthon>> prefixes=new ArrayList<>();BigInteger count=BigInteger.ONE;
        // Sorting defines a reproducible prefix. Only the first MAX_RETAINED synthons per position are needed for any allowed expansion prefix.
        List<Object> signature=new ArrayList<>();signature.add(asset);signature.add(hit.rxn);
        for(var entry:entries) {
            List<SynthonSpace.FragId> sorted=entry.getValue().stream().sorted(Comparator.comparing((SynthonSpace.FragId f)->Objects.toString(f.fragment_id,"")).thenComparing(f->f.idcode)).toList();
            int position=entry.getKey().frag;List<ResultService.Synthon> prefix=sorted.stream().limit(ResultService.MAX_RETAINED).map(f->new ResultService.Synthon(position,f.fragment_id,f.idcode)).toList();
            count=count.multiply(BigInteger.valueOf(sorted.size()));positions.add(new Position(position,sorted.size(),prefix.stream().limit(3).toList()));prefixes.add(prefix);
            signature.add(sorted.stream().map(f->List.of(Objects.toString(f.fragment_id,""),f.idcode)).toList());
        }
        require(!positions.isEmpty(),"Empty combinatorial family");
        String id="family-"+ResultService.hash(WorkspaceStore.JSON.valueToTree(signature).toString());
        return new FamilyData(new Family(id,run,asset,hit.rxn,count.toString(),positions),prefixes);
    }
    public ResultService.Page<Family> families(String run,int offset,int limit)throws IOException {
        require(offset>=0&&offset<=1000&&limit>=1&&limit<=200,"Invalid family page bounds");
        Path file=results.directory(run).resolve("families.json");List<Family> rows=Files.isRegularFile(file)?WorkspaceStore.JSON.readValue(file.toFile(),new TypeReference<>(){}):List.of();
        return new ResultService.Page<>(rows.subList(Math.min(offset,rows.size()),Math.min(offset+limit,rows.size())),offset,limit,rows.size(),"Retained combinatorial families; counts may overlap. No concrete products have been implicitly expanded.");
    }
    public JobService.Job expand(long expected,String actor,String parentRun,List<String> familyIds,int maxProducts,boolean aligned)throws Exception {
        require(familyIds!=null&&!familyIds.isEmpty()&&familyIds.size()<=1000,"Select 1–1000 families");
        require(maxProducts>=1&&maxProducts<=ResultService.MAX_RETAINED,"Expansion attempts must be 1–"+ResultService.MAX_RETAINED);
        var parent=results.getRun(parentRun);require(parent.mode().equals("substructure"),"Expansion requires a substructure run");require(parent.state().equals("completed"),"Wait for family search to finish");
        List<FamilyData> families=new ArrayList<>();for(String id:familyIds.stream().distinct().toList())families.add(readFamily(parentRun,id));
        var settings=Map.of("parentRunId",parentRun,"families",familyIds,"maxProducts",maxProducts,"alignedDepictions",aligned,"policy","Position order then synthon ID/IDCode lexicographic prefix; global attempted-product cap; exact tuple deduplication.");
        synchronized(store){var run=results.create(expected,actor,"bounded_expansion",parent.query(),parent.assetIds(),settings,
                "Bounded prefix from selected retained families; at most "+maxProducts+" assembly attempts. Not exhaustive over all represented combinations.");
            return jobs.submit(store.snapshot().revision,actor,"bounded_expansion",run.id(),settings,List.of(),ctx->{
                ctx.stage("Assembling bounded products");List<ResultService.Hit> hits=new ArrayList<>();Set<String> tuples=new HashSet<>();int attempted=0,rejected=0,duplicate=0;
                for(FamilyData data:families) {
                    int[] cursor=new int[data.prefix().size()];boolean available=data.prefix().stream().noneMatch(List::isEmpty);
                    while(available&&attempted<maxProducts) {
                        ctx.check();List<ResultService.Synthon> synthons=new ArrayList<>();for(int p=0;p<cursor.length;p++)synthons.add(data.prefix().get(p).get(cursor[p]));
                        String tuple=WorkspaceStore.JSON.valueToTree(List.of(data.summary().assetId(),data.summary().reaction(),synthons)).toString();attempted++;
                        if(tuples.add(tuple))try {
                            List<StereoMolecule> parts=new ArrayList<>();for(var synthon:synthons)parts.add(new IDCodeParser().getCompactMolecule(synthon.idcode()));
                            StereoMolecule molecule=SynthonAssembler.assembleSynthons(parts);
                            var hit=ResultService.makeHit(run.id(),tuple,data.summary().assetId(),data.summary().id(),data.summary().reaction(),synthons,molecule,run.query().document(),aligned,Map.of(),
                                    Map.of("parentRunId",parentRun,"familyId",data.summary().id(),"sourceSnapshotId",assets.get(data.summary().assetId()).sourceSnapshotId(),"expansionPolicy",settings.get("policy")));
                            if(hit.matchedAtoms().isEmpty())rejected++;else hits.add(hit);
                        }catch(WorkbenchException|IllegalArgumentException e){rejected++;}else duplicate++;
                        int p=cursor.length-1;while(p>=0&&++cursor[p]>=data.prefix().get(p).size()){cursor[p]=0;p--;}available=p>=0;
                    }
                    if(attempted>=maxProducts)break;
                }
                var diagnostics=Map.of("attemptedProducts",attempted,"rejectedProducts",rejected,"duplicateTuples",duplicate,"retainedOccurrences",hits.size(),"attemptLimitReached",attempted==maxProducts);
                results.complete(run.id(),hits,diagnostics);return diagnostics;
            });}
    }
    private FamilyData readFamily(String run,String id)throws IOException {
        require(id!=null&&id.matches("family-[a-f0-9]{64}"),"Invalid family ID");Path file=results.directory(run).resolve("families").resolve(id+".json.gz");
        if(!Files.isRegularFile(file))throw new WorkbenchException("family_not_found",404,"Family is not retained in this run");
        try(var in=new GZIPInputStream(Files.newInputStream(file))){return WorkspaceStore.JSON.readValue(in,FamilyData.class);}
    }
}
