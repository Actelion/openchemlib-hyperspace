package com.idorsia.research.chem.hyperspace.workbench;

import com.idorsia.research.chem.hyperspace3d.index.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.*;
import java.nio.file.*;
import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/** Exports the real public test fingerprints into the existing Python supplier-cache wire format. */
final class PythonCacheFixture {
    static Path copy(Path nativeCache,Path output,DeepSpaceModelBundle model,CompactSkelSpheresModelBundle compact)throws Exception {
        Files.createDirectories(output);Path shard=output.resolve("shard_000000");Files.createDirectories(shard);
        List<float[]> base=new ArrayList<>(),projection=new ArrayList<>();List<MoleculeFingerprintMetadata> rows=new ArrayList<>();
        try(var source=MoleculeFingerprintDataSources.open(nativeCache)){
            for(int s=0;s<source.shardCount();s++)try(var a=source.openShard(s,MoleculeFingerprintColumn.BASE_128);var b=source.openShard(s,MoleculeFingerprintColumn.COMPACT_16)){
                for(long row=0;row<a.recordCount();row++){base.add(a.readVector(row));projection.add(b.readVector(row));var ref=new MoleculeVectorReference(s,row);rows.add(source.resolveMetadata(List.of(ref)).get(ref));}
            }
        }
        Path a=shard.resolve("fingerprints_128d_fp16.bin"),b=shard.resolve("fingerprints_16d_fp16.bin"),table=shard.resolve("records.tsv.gz");write(a,base,128);write(b,projection,16);
        try(var out=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(Files.newOutputStream(table)),java.nio.charset.StandardCharsets.UTF_8))){
            out.write("local_index\tsource_row\tmolecule_id\tsmiles\tcanonical_smiles\tnum_atoms\n");for(int n=0;n<rows.size();n++){var r=rows.get(n);out.write(n+"\t"+r.sourceRow()+"\t"+r.moleculeId()+"\t"+r.smiles()+"\t"+r.canonicalSmiles()+"\t"+r.heavyAtomCount()+"\n");}
        }
        var provenance=Map.of("integrated_checkpoint_sha256",model.manifest().foundationCheckpointSha256,"predictor_checkpoint_sha256",model.manifest().predictorCheckpointSha256,"projection_checkpoint_sha256",compact.manifest().projectionCheckpointSha256);
        int count=rows.size();Map<String,Object> sm=new LinkedHashMap<>();sm.put("artifact_type","deepspace7_streamed_supplier_fingerprint_shard_v1");sm.put("shard_index",0);sm.put("source_row_start",0);sm.put("source_row_stop",count);sm.put("source_rows",count);sm.put("valid_molecules",count);sm.put("rejected_molecules",0);sm.put("provenance",provenance);
        sm.put("fingerprints",Map.of("128d",vector(a,count,128),"16d",vector(b,count,16)));sm.put("records",Map.of("path",table.getFileName().toString(),"sha256",SharedCandidatePoseBundle.sha256(table)));
        WorkspaceStore.JSON.writeValue(shard.resolve("manifest.json").toFile(),sm);Files.writeString(shard.resolve(".complete"),"complete\n");
        var root=new LinkedHashMap<String,Object>();root.put("artifact_type","deepspace7_streamed_supplier_fingerprints_v1");root.put("complete_source",true);root.put("fingerprints",Map.of("128d",Map.of("bytes_per_molecule",256,"dtype","float16"),"16d",Map.of("bytes_per_molecule",32,"dtype","float16")));root.put("provenance",provenance);root.put("source_rows",count);root.put("valid_molecules",count);root.put("rejected_molecules",0);root.put("shards",1);
        WorkspaceStore.JSON.writeValue(output.resolve("manifest.json").toFile(),root);return output;
    }
    private static void write(Path path,List<float[]> rows,int dim)throws Exception{var bytes=ByteBuffer.allocate(rows.size()*dim*2).order(ByteOrder.LITTLE_ENDIAN);for(float[] row:rows)for(float v:row)bytes.putShort(Float.floatToFloat16(v));Files.write(path,bytes.array());}
    private static Map<String,Object> vector(Path path,int count,int dim)throws Exception{return Map.of("path",path.getFileName().toString(),"dtype","float16","shape",List.of(count,dim),"sha256",SharedCandidatePoseBundle.sha256(path));}
}
