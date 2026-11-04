package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.*;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.time.Instant;
import java.util.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Bounded retained results, immutable occurrences and exact canonical OCL graph identity. */
public final class ResultService {
    public static final int MAX_RETAINED=20_000;
    public static final String IDENTITY_POLICY="Canonical OpenChemLib IDCode including charge and defined stereochemistry; no tautomer, salt or protonation normalization. Source occurrences remain separate.";
    public record Synthon(int position,String id,String idcode){}
    public record Hit(String id,String runId,String moleculeId,QueryChemistry.Document document,String smiles,
                      String assetId,String familyId,String reaction,List<Synthon> synthons,
                      Map<String,Double> scores,Map<String,Double> properties,List<Integer> matchedAtoms,JsonNode provenance){}
    public record Run(String id,String mode,String queryId,WorkspaceStore.QueryRevision query,List<String> assetIds,
                      JsonNode settings,String actor,String createdAt,String scope,String state,int retainedCount,
                      int uniqueMoleculeCount,JsonNode diagnostics){}
    public record CollectionRecord(String id,String name,String notes,List<String> hitIds,String actor,String createdAt,String updatedAt){}
    public record Page<T>(List<T> items,int offset,int limit,int total,String scope){}
    public record Filters(Double minWeight,Double maxWeight,Integer minAtoms,Integer maxAtoms){}
    public record Comparison(String leftRunId,String rightRunId,int shared,int onlyLeft,int onlyRight,
                             List<String> sharedIdentities,List<String> leftIdentities,List<String> rightIdentities,String identityPolicy){}
    private final WorkspaceStore store;
    public ResultService(WorkspaceStore store)throws IOException{this.store=store;Files.createDirectories(store.directory().resolve("runs"));}
    public Run create(long expected,String actor,String mode,WorkspaceStore.QueryRevision query,List<String> assets,Object settings,String scope) {
        Run run=new Run("run-"+UUID.randomUUID(),mode,query.id(),query,List.copyOf(assets),WorkspaceStore.JSON.valueToTree(settings),actor,
                Instant.now().toString(),scope,"queued",0,0,WorkspaceStore.JSON.createObjectNode());
        return store.mutate(expected,actor,"run.create",s->{s.runs.put(run.id(),WorkspaceStore.JSON.valueToTree(run));s.view="screening";s.viewEntityId=run.id();return run;});
    }
    public Run getRun(String id){JsonNode n=store.snapshot().runs.get(id);if(n==null)throw new WorkbenchException("run_not_found",404,"Unknown run: "+id);return WorkspaceStore.JSON.convertValue(n,Run.class);}
    public Path directory(String id)throws IOException {getRun(id);require(id.matches("run-[a-f0-9-]{36}"),"Invalid run ID");Path d=store.directory().resolve("runs").resolve(id);Files.createDirectories(d);return d;}
    public void state(String id,String state,Object diagnostics){store.internal("run."+state,s->{Run r=WorkspaceStore.JSON.convertValue(s.runs.get(id),Run.class);
        require(r!=null,"Run is absent");s.runs.put(id,WorkspaceStore.JSON.valueToTree(new Run(r.id(),r.mode(),r.queryId(),r.query(),r.assetIds(),r.settings(),r.actor(),r.createdAt(),r.scope(),state,r.retainedCount(),r.uniqueMoleculeCount(),WorkspaceStore.JSON.valueToTree(diagnostics))));return null;});}
    public void complete(String id,List<Hit> hits,Object diagnostics)throws IOException {
        require(hits.size()<=MAX_RETAINED,"Retained result limit exceeded");
        require(hits.stream().allMatch(h->h.runId().equals(id)),"Result belongs to another run");
        require(hits.stream().map(Hit::id).distinct().count()==hits.size(),"Duplicate hit occurrence IDs");
        Path target=directory(id).resolve("hits.json"),part=target.resolveSibling("hits.json.partial");
        WorkspaceStore.JSON.writeValue(part.toFile(),hits);Files.move(part,target,StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
        store.internal("run.completed",s->{Run r=WorkspaceStore.JSON.convertValue(s.runs.get(id),Run.class);
            s.runs.put(id,WorkspaceStore.JSON.valueToTree(new Run(r.id(),r.mode(),r.queryId(),r.query(),r.assetIds(),r.settings(),r.actor(),r.createdAt(),r.scope(),"completed",hits.size(),
                    (int)hits.stream().map(Hit::moleculeId).distinct().count(),WorkspaceStore.JSON.valueToTree(diagnostics))));return null;});
    }
    public List<Hit> all(String runId)throws IOException {
        Path file=directory(runId).resolve("hits.json");if(!Files.isRegularFile(file))return List.of();
        List<Hit> rows=WorkspaceStore.JSON.readValue(file.toFile(),new TypeReference<>(){});
        require(rows.size()<=MAX_RETAINED,"Stored results exceed configured bound");return rows;
    }
    public Page<Hit> page(String runId,int offset,int limit,String sort,boolean descending,Filters filters)throws IOException {
        bounds(offset,limit);List<Hit> rows=new ArrayList<>(all(runId));
        var scoreFields=new HashSet<String>();rows.forEach(h->scoreFields.addAll(h.scores().keySet()));
        if(filters!=null)rows.removeIf(h->!accept(h,filters));
        if(sort!=null&&!sort.isBlank()) {
            require(Set.of("weight","heavyAtoms","rotors").contains(sort)||getRun(runId).settings().path("scoreTargets").toString().contains("\""+sort+"\"")
                    ||scoreFields.contains(sort),"Unknown sortable field");
            Comparator<Hit> comparator=Comparator.comparingDouble(h->h.properties().getOrDefault(sort,h.scores().getOrDefault(sort,Double.NaN)));
            if(descending)comparator=comparator.reversed();rows.sort(comparator.thenComparing(Hit::id));
        }
        return new Page<>(List.copyOf(rows.subList(Math.min(offset,rows.size()),Math.min(offset+limit,rows.size()))),offset,limit,rows.size(),
                "Retained results only. Filtering or sorting this bounded list does not establish the global top-K satisfying a filter.");
    }
    public Hit hit(String id)throws IOException {
        require(id!=null&&id.matches("run-[a-f0-9-]{36}/hit-[a-f0-9]{64}"),"Invalid hit ID");String run=id.substring(0,id.indexOf('/'));
        return all(run).stream().filter(h->h.id().equals(id)).findFirst().orElseThrow(()->new WorkbenchException("hit_not_found",404,"Unknown hit occurrence"));
    }
    public void validateHitIds(List<String> ids,int maximum)throws IOException {
        require(ids!=null&&ids.size()<=maximum,"Too many selected source occurrences");
        var byRun=new HashMap<String,Set<String>>();
        for(String id:ids){require(id!=null&&id.matches("run-[a-f0-9-]{36}/hit-[a-f0-9]{64}"),"Invalid hit ID");
            String run=id.substring(0,id.indexOf('/'));
            if(!byRun.containsKey(run))byRun.put(run,new HashSet<>(all(run).stream().map(Hit::id).toList()));
            require(byRun.get(run).contains(id),"Unknown selected source occurrence");
        }
    }
    public CollectionRecord saveCollection(long expected,String actor,String id,String name,String notes,List<String> hitIds)throws IOException {
        require(name!=null&&!name.isBlank()&&name.length()<=200,"Collection name required (max 200 characters)");require(notes==null||notes.length()<=20_000,"Notes exceed limit");
        if(hitIds==null){require(id!=null,"New collections require explicit occurrence IDs");hitIds=collection(id).hitIds();}
        require(hitIds.size()<=MAX_RETAINED,"Collection exceeds retained-result limit");
        List<String> unique=hitIds.stream().distinct().toList();Map<String,Set<String>> runs=new HashMap<>();
        for(String hit:unique){require(hit!=null&&hit.contains("/"),"Invalid hit ID");String r=hit.substring(0,hit.indexOf('/'));
            if(!runs.containsKey(r))runs.put(r,new HashSet<>(all(r).stream().map(Hit::id).toList()));require(runs.get(r).contains(hit),"Unknown selected hit");}
        return store.mutate(expected,actor,"collection.save",s->{if(id!=null)require(s.collections.containsKey(id),"Collection is absent");
            String now=Instant.now().toString();String created=id==null?now:s.collections.get(id).path("createdAt").asText(now);
            var c=new CollectionRecord(id==null?"collection-"+UUID.randomUUID():id,name,Objects.toString(notes,""),unique,actor,created,now);
            s.collections.put(c.id(),WorkspaceStore.JSON.valueToTree(c));s.view="saved";s.viewEntityId=c.id();return c;});
    }
    public CollectionRecord collection(String id){var n=store.snapshot().collections.get(id);if(n==null)throw new WorkbenchException("collection_not_found",404,"Unknown collection");return WorkspaceStore.JSON.convertValue(n,CollectionRecord.class);}
    public List<Hit> collectionHits(String id)throws IOException {var selected=collection(id).hitIds();var byRun=new HashMap<String,Map<String,Hit>>();var rows=new ArrayList<Hit>();
        for(String ref:selected){String run=ref.substring(0,ref.indexOf('/'));if(!byRun.containsKey(run)){Map<String,Hit> map=new HashMap<>();for(Hit h:all(run))map.put(h.id(),h);byRun.put(run,map);}Hit h=byRun.get(run).get(ref);require(h!=null,"Collection references a missing occurrence");rows.add(h);}return rows;}
    public Comparison compare(String left,String right,int sampleLimit)throws IOException {
        require(sampleLimit>=0&&sampleLimit<=200,"Comparison identity sample limit is 0–200");
        Set<String>a=new TreeSet<>(all(left).stream().map(Hit::moleculeId).toList()),b=new TreeSet<>(all(right).stream().map(Hit::moleculeId).toList());
        Set<String> shared=new TreeSet<>(a);shared.retainAll(b);Set<String> onlyA=new TreeSet<>(a);onlyA.removeAll(b);Set<String> onlyB=new TreeSet<>(b);onlyB.removeAll(a);
        return new Comparison(left,right,shared.size(),onlyA.size(),onlyB.size(),shared.stream().limit(sampleLimit).toList(),onlyA.stream().limit(sampleLimit).toList(),onlyB.stream().limit(sampleLimit).toList(),IDENTITY_POLICY);
    }
    public List<String> diverse(String runId,int count)throws IOException {
        List<Hit> rows=new ArrayList<>(all(runId).stream().collect(java.util.stream.Collectors.toMap(Hit::moleculeId,h->h,(a,b)->a,LinkedHashMap::new)).values());
        require(count>=1&&count<=500,"Diversity selection is bounded to 1–500 molecules");if(rows.isEmpty())return List.of();
        long[][] vectors=new long[rows.size()][];var descriptor=new DescriptorHandlerLongFFP512();
        for(int i=0;i<rows.size();i++)vectors[i]=descriptor.createDescriptor(QueryChemistry.parse(rows.get(i).document()));
        int[] chosen=new DiversitySelector().select(vectors,Math.min(count,rows.size()));return Arrays.stream(chosen).mapToObj(i->rows.get(i).id()).toList();
    }
    public Map<String,Object> export(String runId,String collectionId)throws Exception {
        require((runId==null)!=(collectionId==null),"Choose exactly one run or collection");
        List<Hit> rows=runId==null?collectionHits(collectionId):all(runId);String id="export-"+UUID.randomUUID();Path dir=store.directory().resolve("exports").resolve(id);Files.createDirectories(dir);
        Map<String,Run> runs=new LinkedHashMap<>();for(Hit h:rows)runs.putIfAbsent(h.runId(),getRun(h.runId()));
        if(runId!=null)runs.putIfAbsent(runId,getRun(runId));
        Map<String,Object> provenance=new LinkedHashMap<>();provenance.put("formatVersion",1);provenance.put("workspaceId",store.snapshot().workspaceId);
        provenance.put("exportedAt",Instant.now().toString());provenance.put("identityPolicy",IDENTITY_POLICY);provenance.put("runs",runs);provenance.put("occurrenceCount",rows.size());
        provenance.put("uniqueMoleculeCount",rows.stream().map(Hit::moleculeId).distinct().count());provenance.put("hits",rows);
        if(collectionId!=null)provenance.put("collection",collection(collectionId));
        WorkspaceStore.JSON.writeValue(dir.resolve("provenance.json").toFile(),provenance);
        try(var tsv=Files.newBufferedWriter(dir.resolve("results.tsv"));var sdf=Files.newBufferedWriter(dir.resolve("results.sdf"))) {
            tsv.write("hit_id\tmolecule_id\tsmiles\tasset_id\treaction\tsynthons_json\tscores_json\n");
            for(Hit h:rows) {
                tsv.write(String.join("\t",List.of(h.id(),h.moleculeId(),h.smiles(),h.assetId(),Objects.toString(h.reaction(),""),compact(h.synthons()),compact(h.scores())).stream().map(ResultService::tsv).toList())+"\n");
                sdf.write(new MolfileCreator(QueryChemistry.parse(h.document())).getMolfile());
                for(var field:Map.of("HIT_ID",h.id(),"MOLECULE_ID",h.moleculeId(),"ASSET_ID",h.assetId(),"SCORES",compact(h.scores()),"PROVENANCE",compact(h.provenance()),"SYNTHONS",compact(h.synthons())).entrySet())
                    sdf.write("\n>  <"+field.getKey()+">\n"+field.getValue()+"\n\n");
                sdf.write("$$$$\n");
            }
        }
        Files.writeString(dir.resolve("summary.md"),"# Hyperspace workbench export\n\n"+rows.size()+" source occurrences; "+rows.stream().map(Hit::moleculeId).distinct().count()+" unique molecular identities.\n\n"+IDENTITY_POLICY+"\n\n"+
                "Scores retain each run's method and objective; scores from different methods are not averaged.\n\n"+runs.values().stream().map(r->"- `"+r.id()+"` ("+r.mode()+"): "+r.scope()+"; query `"+r.queryId()+"`.\n").reduce("",String::concat));
        return Map.of("exportId",id,"directory",dir.toString(),"files",List.of("results.sdf","results.tsv","summary.md","provenance.json"),"occurrenceCount",rows.size());
    }
    public static Hit makeHit(String run,String occurrenceKey,String asset,String family,String reaction,List<Synthon> synthons,
                              StereoMolecule molecule,QueryChemistry.Document query,boolean aligned,Map<String,Double> scores,Object provenance) {
        molecule.setFragment(false);molecule.ensureHelperArrays(Molecule.cHelperCIP);
        if(aligned&&query!=null) {
            var q=QueryChemistry.parse(query);q.setFragment(true);var inventor=new CoordinateInventor();
            if(!QueryChemistry.inspect("depiction",query).bonds().stream().anyMatch(QueryChemistry.Bond::bridge))
                inventor.setCustomTemplateList(List.of(new InventorTemplate(q,new DescriptorHandlerLongFFP512().createDescriptor(q),true)));
            inventor.invent(molecule);
        }
        var document=QueryChemistry.canonical(molecule);var canonical=QueryChemistry.parse(document);List<Integer> match=List.of();
        if(query!=null)try {
            var q=QueryChemistry.parse(query);q.setFragment(true);
            var search=new SSSearcher();search.setMol(q,canonical);
            if(search.findFragmentInMolecule(SSSearcher.cCountModeFirstMatch,SSSearcher.cDefaultMatchMode)>0)
                match=Arrays.stream(search.getMatchList().get(0)).boxed().toList();
        }catch(Exception e){throw new WorkbenchException("highlight_failed",400,"Cannot interpret query for highlighting: "+e.getMessage());}
        return new Hit(run+"/hit-"+hash(occurrenceKey),run,"mol-"+hash(document.idcode()),document,IsomericSmilesCreator.createSmiles(canonical),asset,family,reaction,List.copyOf(synthons),Map.copyOf(scores),
                Map.of("heavyAtoms",(double)canonical.getAtoms(),"weight",(double)canonical.getMolweight(),"rotors",(double)canonical.getRotatableBondCount()),match,WorkspaceStore.JSON.valueToTree(provenance));
    }
    public static String hash(String value){try{return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value.getBytes(StandardCharsets.UTF_8)));}catch(Exception e){throw new IllegalStateException(e);}}
    private static String compact(Object value)throws IOException{return WorkspaceStore.JSON.writer().without(com.fasterxml.jackson.databind.SerializationFeature.INDENT_OUTPUT).writeValueAsString(value);}
    private static String tsv(String value){return value.replace("\t"," ").replace("\r"," ").replace("\n"," ");}
    private static void bounds(int offset,int limit){require(offset>=0&&offset<=MAX_RETAINED&&limit>=1&&limit<=200,"Page offset must be nonnegative and limit 1–200");}
    private static boolean accept(Hit h,Filters f){double w=h.properties().get("weight"),a=h.properties().get("heavyAtoms");return(f.minWeight()==null||w>=f.minWeight())&&(f.maxWeight()==null||w<=f.maxWeight())&&(f.minAtoms()==null||a>=f.minAtoms())&&(f.maxAtoms()==null||a<=f.maxAtoms());}
}
