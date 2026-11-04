package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.*;
import com.idorsia.research.chem.hyperspace.*;
import com.idorsia.research.chem.hyperspace.rawspace.*;
import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.io.*;
import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.zip.GZIPInputStream;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Catalogue and explicit resource lifecycle. Registration never materializes a large index. */
public final class AssetService implements AutoCloseable {
    public enum Kind {RAW_SPACE,DOWNSAMPLED_SPACE,SUBSTRUCTURE_INDEX,SIMILARITY_INDEX,PRODUCT_INDEX,FLAT_INDEX,PYTHON_CACHE,GRAPH_MODEL,COMPACT_MODEL,POSE_MODEL}
    public record Detection(String path,Kind kind,List<String> capabilities,Long recordCount,String coverage,JsonNode metadata,List<String> diagnostics){}
    public record Asset(String id,String libraryId,String sourceSnapshotId,String name,String path,Kind kind,
                        List<String> capabilities,Long recordCount,String coverage,JsonNode metadata,
                        String status,String error,String registeredAt){}
    public record Library(String id,String name){}
    public record SourceSnapshot(String id,String libraryId,String version,String provenance){}
    private static final ObjectMapper JSON=WorkspaceStore.JSON;
    private final WorkspaceStore store;
    private final Map<String,Resource> resources=new HashMap<>();
    private final Map<String,String> bindings=new HashMap<>();
    private static final class Resource {final Object value;int active;Resource(Object value){this.value=value;}}
    public final class Lease implements AutoCloseable {
        private final Resource resource;private boolean closed;
        private Lease(Resource resource){this.resource=resource;resource.active++;}
        public Object value(){if(closed)throw new IllegalStateException("Resource lease is closed");return resource.value;}
        @Override public void close(){synchronized(AssetService.this){if(!closed){closed=true;resource.active--;}}}
    }
    public AssetService(WorkspaceStore store) {
        this.store=store;
        if(store.snapshot().assets.values().stream().anyMatch(a->Set.of("ready","loading").contains(a.path("status").asText())))
            store.internal("resources.owner_restarted",s->{s.assets.replaceAll((id,node)->{
                var a=JSON.convertValue(node,Asset.class);return JSON.valueToTree(withStatus(a,"unloaded",null));
            });return null;});
    }
    public Detection inspect(Path input,Kind requested)throws Exception {
        Path path=input.toRealPath();List<String> diagnostics=new ArrayList<>();JsonNode metadata=JSON.createObjectNode();Kind kind=requested;Long count=null;
        Path manifest=Files.isDirectory(path)?path.resolve("manifest.json"):null;
        if(manifest!=null&&Files.isRegularFile(manifest)) {
            require(Files.size(manifest)<=16*1024*1024,"Manifest exceeds 16 MiB inspection bound");
            metadata=JSON.readTree(manifest.toFile());
            String type=metadata.path("artifactType").asText(metadata.path("artifact_type").asText());
            Kind detected=switch(type) {
                case "hyperspace-molecule-fingerprint-index"->Kind.FLAT_INDEX;
                case "deepspace7_streamed_supplier_fingerprints_v1"->Kind.PYTHON_CACHE;
                case "hyperspace3d-product-fingerprint-index"->Kind.PRODUCT_INDEX;
                case "deepspace7-v1-onnx-bundle"->Kind.GRAPH_MODEL;
                case "deepspace7-skelspheres16-onnx-bundle"->Kind.COMPACT_MODEL;
                case "deepspace7-shared-candidate-pose-onnx-bundle"->Kind.POSE_MODEL;
                default->null;
            };
            require(detected!=null,"Unsupported prepared-asset manifest type: "+type);
            require(kind==null||kind==detected,"Requested asset type differs from its manifest");kind=detected;
            switch(kind) {
                case FLAT_INDEX -> {var m=MoleculeFingerprintIndexReader.loadManifest(path);count=m.recordCount;}
                case PYTHON_CACHE -> {try(var source=MoleculeFingerprintDataSources.open(path)){count=source.recordCount();}}
                case PRODUCT_INDEX -> {var m=JSON.treeToValue(metadata,ProductFingerprintIndexManifest.class);m.validate();count=m.recordCount;}
                case GRAPH_MODEL -> JSON.treeToValue(metadata,DeepSpaceModelManifest.class).validate();
                case COMPACT_MODEL -> JSON.treeToValue(metadata,CompactSkelSpheresManifest.class).validate();
                case POSE_MODEL -> require(metadata.path("architecture").asText().equals("shared_candidate_pose_v1")
                        &&metadata.path("atomFeatureDim").asInt()==50&&metadata.path("embeddingDim").asInt()==128,"Unsupported pose-model contract");
                default -> throw new IllegalArgumentException("Unexpected manifest kind");
            }
            diagnostics.add("Manifest/startup checks only; no full fingerprint payload scan or hash pass.");
        } else {
            require(Files.isRegularFile(path),"Expected a prepared asset file or a directory with manifest.json");
            if(kind==null) {
                String name=path.getFileName().toString();
                if(name.endsWith(".rawspace.gz")||name.endsWith(".rawspace.json"))kind=Kind.RAW_SPACE;
                else if(name.endsWith("_similarity3.data"))kind=Kind.SIMILARITY_INDEX;
                else if(name.endsWith(".data"))kind=Kind.SUBSTRUCTURE_INDEX;
            }
            require(kind!=null&&Set.of(Kind.RAW_SPACE,Kind.DOWNSAMPLED_SPACE,Kind.SUBSTRUCTURE_INDEX,Kind.SIMILARITY_INDEX).contains(kind),"Specify a supported rawspace or index file type");
            if(Set.of(Kind.SUBSTRUCTURE_INDEX,Kind.SIMILARITY_INDEX).contains(kind)) {
                try(var in=Files.newInputStream(path)){require(in.read()==0x1f&&in.read()==0x8b,"Expected gzip-serialized Hyperspace index");}
            }
            diagnostics.add("File exists; chemical contents/counts are validated during explicit loading, not registration.");
        }
        String coverage=switch(kind) {
            case RAW_SPACE -> "Combinatorial combinations defined by the raw synthon sets; not enumerated rows.";
            case DOWNSAMPLED_SPACE -> "Reduced representative synthon sets; does not cover all full-space combinations.";
            case SUBSTRUCTURE_INDEX,SIMILARITY_INDEX -> "Combinatorial coverage of this particular prepared index; counts available after loading.";
            case PRODUCT_INDEX -> "Sampled complete product tuples only; not exhaustive combinatorial coverage.";
            case FLAT_INDEX,PYTHON_CACHE -> "Accepted molecule rows in this cache; rejected/source rows and provenance remain distinct.";
            default -> "Model parameters; no molecular-library coverage.";
        };
        return new Detection(path.toString(),kind,capabilities(kind),count,coverage,metadata,diagnostics);
    }
    public Asset register(long expected,String actor,Detection detection,String name,String libraryId,String libraryName,String snapshotId,String sourceVersion,String provenance) {
        require(name!=null&&!name.isBlank()&&name.length()<=200,"Asset name required (max 200 characters)");
        return store.mutate(expected,actor,"asset.register",state->{
            String library=libraryId;
            if(library==null||library.isBlank()){library="library-"+UUID.randomUUID();state.libraries.put(library,JSON.valueToTree(new Library(library,libraryName==null?name:libraryName)));}
            else require(state.libraries.containsKey(library),"Logical library is absent");
            String snapshot=snapshotId;
            if(snapshot==null||snapshot.isBlank()) {
                snapshot="source-"+UUID.randomUUID();state.sourceSnapshots.put(snapshot,JSON.valueToTree(new SourceSnapshot(snapshot,library,
                        sourceVersion==null?"unspecified":sourceVersion,provenance==null?"User-registered local asset; no row-equivalence asserted.":provenance)));
            } else {
                require(state.sourceSnapshots.containsKey(snapshot),"Source snapshot is absent");
                require(state.sourceSnapshots.get(snapshot).path("libraryId").asText().equals(library),"Source snapshot belongs to a different logical library");
            }
            var asset=new Asset("asset-"+UUID.randomUUID(),library,snapshot,name,detection.path(),detection.kind(),detection.capabilities(),
                    detection.recordCount(),detection.coverage(),detection.metadata().deepCopy(),"registered",null,Instant.now().toString());
            state.assets.put(asset.id(),JSON.valueToTree(asset));return asset;
        });
    }
    public Asset get(String id) {
        JsonNode value=store.snapshot().assets.get(id);if(value==null)throw new WorkbenchException("asset_not_found",404,"Unknown asset: "+id);
        return JSON.convertValue(value,Asset.class);
    }
    /** Call from an application job, never from registration or the browser rendering path. */
    public synchronized Asset load(long expected,String actor,String id,String graphModelId)throws Exception {
        Asset asset=get(id);store.mutate(expected,actor,"asset.load",s->{s.assets.put(id,JSON.valueToTree(withStatus(asset,"loading",null)));return null;});
        String key=asset.kind()+":"+Path.of(asset.path()).toRealPath();
        try {
            if(!resources.containsKey(key)) {
                Path path=Path.of(asset.path());Object resource=switch(asset.kind()) {
                    case RAW_SPACE,DOWNSAMPLED_SPACE -> RawSynthonSpaceIO.read(path);
                    case SUBSTRUCTURE_INDEX,SIMILARITY_INDEX -> loadIndex(path,asset.kind());
                    case FLAT_INDEX,PYTHON_CACHE -> MoleculeFingerprintDataSources.open(path);
                    case PRODUCT_INDEX -> JSON.treeToValue(asset.metadata(),ProductFingerprintIndexManifest.class);
                    case GRAPH_MODEL -> DeepSpaceModelBundle.load(path);
                    case COMPACT_MODEL -> CompactSkelSpheresModelBundle.load(path,graph(graphModelId).manifest());
                    case POSE_MODEL -> SharedCandidatePoseBundle.load(path,graph(graphModelId));
                };
                resources.put(key,new Resource(resource));
            }
            bindings.put(id,key);var ready=withStatus(asset,"ready",null);
            store.internal("asset.ready",s->{s.assets.put(id,JSON.valueToTree(ready));return null;});return ready;
        } catch(Exception e) {
            var failed=withStatus(asset,"failed",e.getMessage());store.internal("asset.failed",s->{s.assets.put(id,JSON.valueToTree(failed));return null;});throw e;
        }
    }
    private DeepSpaceModelBundle graph(String id) {
        require(id!=null,"A loaded canonical graph-model asset is required for this bundle");
        Resource resource=resources.get(bindings.get(id));
        require(resource!=null&&resource.value instanceof DeepSpaceModelBundle,"Load the canonical graph-model asset first");
        return (DeepSpaceModelBundle)resource.value;
    }
    public synchronized Lease acquire(String id) {
        Resource resource=resources.get(bindings.get(id));
        if(resource==null)throw new WorkbenchException("resource_not_ready",409,"Explicitly load this asset before searching");
        return new Lease(resource);
    }
    public synchronized Asset unload(long expected,String actor,String id)throws Exception {
        Asset asset=get(id);String key=bindings.get(id);Resource resource=resources.get(key);
        if(resource!=null&&resource.active>0)throw new WorkbenchException("resource_busy",409,"Active jobs hold this resource; cancel or finish them before unloading");
        var result=withStatus(asset,"unloaded",null);
        store.mutate(expected,actor,"asset.unload",s->{s.assets.put(id,JSON.valueToTree(result));return null;});
        bindings.remove(id);
        if(resource!=null&&!bindings.containsValue(key)){resources.remove(key);if(resource.value instanceof AutoCloseable closeable)closeable.close();}
        return result;
    }
    public synchronized Map<String,Object> resourceStatus(String id){Resource resource=resources.get(bindings.get(id));return Map.of("loaded",resource!=null,"activeJobs",resource==null?0:resource.active);}
    private static Object loadIndex(Path path,Kind kind)throws Exception {
        try(var objects=new ObjectInputStream(new GZIPInputStream(Files.newInputStream(path)))) {
            Object value=objects.readObject();
            require(value instanceof SynthonSpace,"File does not contain a supported SynthonSpace index");
            require(kind!=Kind.SIMILARITY_INDEX||value instanceof SynthonSimilaritySpace3,"Expected SynthonSimilaritySpace3");
            ((SynthonSpace)value).initAfterJavaDeserialization();return value;
        }
    }
    private static List<String> capabilities(Kind kind){return switch(kind){
        case RAW_SPACE->List.of("prepare_substructure","downsample","sample_products");
        case DOWNSAMPLED_SPACE->List.of("sample_products");
        case SUBSTRUCTURE_INDEX->List.of("substructure","bounded_expansion");
        case SIMILARITY_INDEX->List.of("legacy_similarity");
        case PRODUCT_INDEX->List.of("sampled_product_3d");
        case FLAT_INDEX,PYTHON_CACHE->List.of("skelspheres_2d","phesa_3d","skelspheres_then_phesa","pose_phesa_3d");
        case GRAPH_MODEL->List.of("graph_encoder","graph_comparator");
        case COMPACT_MODEL->List.of("skelspheres_projection");
        case POSE_MODEL->List.of("pose_encoder","directed_comparator");};}
    private static Asset withStatus(Asset a,String status,String error){return new Asset(a.id(),a.libraryId(),a.sourceSnapshotId(),a.name(),a.path(),a.kind(),a.capabilities(),a.recordCount(),a.coverage(),a.metadata(),status,error,a.registeredAt());}
    @Override public synchronized void close()throws Exception {
        for(Resource resource:resources.values())if(resource.active>0)throw new IllegalStateException("Active jobs must stop before closing resources");
        for(Resource resource:resources.values())if(resource.value instanceof AutoCloseable closeable)closeable.close();resources.clear();bindings.clear();
    }
}
