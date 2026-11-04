package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.*;
import com.idorsia.research.chem.hyperspace3d.feature.DeepSpaceTensorBatch;
import java.nio.FloatBuffer;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;

/** Reusable directed batch scorer for indexed and future assembled candidate embeddings. */
public final class SharedCandidatePoseRuntime implements EmbeddingComparator, AutoCloseable {
    private final DeepSpaceOnnxEnvironment runtime;
    private final OrtSession query, comparator;
    private final float ppWeight;

    public SharedCandidatePoseRuntime(DeepSpaceOnnxEnvironment runtime, SharedCandidatePoseBundle bundle, double ppWeight) {
        if (!Double.isFinite(ppWeight) || ppWeight < 0 || ppWeight > 1)
            throw new IllegalArgumentException("pose ppWeight must be in [0,1]");
        this.runtime = runtime; this.ppWeight = (float)ppWeight;
        OrtSession q = null, c = null;
        try {
            q=runtime.open(bundle.directory().resolve("query.onnx"));
            c=runtime.open(bundle.directory().resolve("comparator.onnx"));
            if (!q.getInputNames().equals(Set.of("atom_x","pair_x","atom_mask","distances"))
                    || !q.getOutputNames().equals(Set.of("embedding"))
                    || !c.getInputNames().equals(Set.of("query_embedding","candidate_embedding","pp_weight"))
                    || !c.getOutputNames().equals(Set.of("scores")))
                throw new IllegalArgumentException("pose ONNX names mismatch");
            OnnxContract.requireInput(q,"atom_x",OnnxJavaType.FLOAT,32,50);
            OnnxContract.requireInput(q,"pair_x",OnnxJavaType.FLOAT,32,32,36);
            OnnxContract.requireInput(q,"atom_mask",OnnxJavaType.BOOL,32);
            OnnxContract.requireInput(q,"distances",OnnxJavaType.FLOAT,32,32);
            OnnxContract.requireOutput(q,"embedding",OnnxJavaType.FLOAT,128);
            OnnxContract.requireInput(c,"query_embedding",OnnxJavaType.FLOAT,128);
            OnnxContract.requireInput(c,"candidate_embedding",OnnxJavaType.FLOAT,128);
            OnnxContract.requireInput(c,"pp_weight",OnnxJavaType.FLOAT,1);
            OnnxContract.requireOutput(c,"scores",OnnxJavaType.FLOAT,3);
        } catch (RuntimeException e) {
            try { if (q!=null) q.close(); if(c!=null)c.close(); } catch(OrtException suppressed) { e.addSuppressed(suppressed); }
            throw e;
        }
        query=q; comparator=c;
    }
    public synchronized float[][] encode(DeepSpaceTensorBatch batch, float[] distances) {
        int b=batch.batchSize();
        if(distances.length!=b*32*32) throw new IllegalArgumentException("pose distances shape mismatch");
        for(float value:distances) if(!Float.isFinite(value) || value<0)
            throw new IllegalArgumentException("pose distances must be finite and nonnegative");
        boolean[][] mask=new boolean[b][32];
        for(int i=0;i<b;i++)System.arraycopy(batch.atomMask(),i*32,mask[i],0,32);
        try(var atoms=tensor(batch.atomFeatures(),b,32,50);
            var pairs=tensor(batch.pairFeatures(),b,32,32,36);
            var masks=OnnxTensor.createTensor(runtime.environment(),mask);
            var d=tensor(distances,b,32,32);
            var result=query.run(Map.of("atom_x",atoms,"pair_x",pairs,"atom_mask",masks,"distances",d))) {
            return checked((float[][])result.get(0).getValue(),b,128);
        } catch(OrtException e){throw new DeepSpaceInferenceException("pose query inference failed",e);}
    }
    @Override public float[][] compare(float[] q,float[][] candidates) {
        float[] flat=new float[candidates.length*128];
        for(int i=0;i<candidates.length;i++) {
            if(candidates[i].length!=128)throw new IllegalArgumentException("candidate must be 128D");
            System.arraycopy(candidates[i],0,flat,i*128,128);
        }
        return compareFlat(q,flat,candidates.length);
    }
    @Override public synchronized float[][] compareFlat(float[] q,float[] candidates,int count) {
        if(q.length!=128 || count<0 || candidates.length!=Math.multiplyExact(count,128))
            throw new IllegalArgumentException("pose comparison shape mismatch");
        if(count==0)return new float[0][3];
        float[] repeated=new float[count*128],weights=new float[count]; Arrays.fill(weights,ppWeight);
        for(int i=0;i<count;i++)System.arraycopy(q,0,repeated,i*128,128);
        try(var qt=tensor(repeated,count,128);var ct=tensor(candidates,count,128);var w=tensor(weights,count,1);
            var result=comparator.run(Map.of("query_embedding",qt,"candidate_embedding",ct,"pp_weight",w))) {
            return checked((float[][])result.get(0).getValue(),count,3);
        }catch(OrtException e){throw new DeepSpaceInferenceException("pose comparator inference failed",e);}
    }
    private OnnxTensor tensor(float[] values,long... shape)throws OrtException {
        return OnnxTensor.createTensor(runtime.environment(),FloatBuffer.wrap(values),shape);
    }
    private static float[][] checked(float[][] values,int rows,int columns) {
        if(values.length!=rows)throw new DeepSpaceInferenceException("pose output batch mismatch");
        for(float[] row:values){
            if(row.length!=columns)throw new DeepSpaceInferenceException("pose output width mismatch");
            for(float value:row)if(!Float.isFinite(value))throw new DeepSpaceInferenceException("nonfinite pose output");
        }
        return values;
    }
    @Override public synchronized void close() {
        try { query.close(); } catch(OrtException e){throw new DeepSpaceInferenceException("closing pose query",e);}
        finally {try{comparator.close();}catch(OrtException e){throw new DeepSpaceInferenceException("closing pose comparator",e);}}
    }
}
