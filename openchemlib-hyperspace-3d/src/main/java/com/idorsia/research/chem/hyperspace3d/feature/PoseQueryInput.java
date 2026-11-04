package com.idorsia.research.chem.hyperspace3d.feature;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace3d.model.SharedCandidatePoseBundle;
import java.nio.file.*;
import java.util.*;

/** Supplied coordinates in angstrom, with an explicit source-to-heavy-atom mapping. */
public record PoseQueryInput(StereoMolecule molecule, DeepSpaceTensorBatch batch,
        float[] distances, Map<String,Object> provenance) {
    /** record is zero based; null is allowed only for a single-record input. */
    public static PoseQueryInput read(Path path, Integer record) throws Exception {
        if(record!=null && record<0)throw new IllegalArgumentException("query record must be zero based and nonnegative");
        String selected=null;
        int count=0;
        try(var reader=Files.newBufferedReader(path)) {
            StringBuilder current=new StringBuilder();
            for(String line;(line=reader.readLine())!=null;) {
                if(line.equals("$$$$")) {
                    if(count==(record==null?0:record))selected=current.toString();
                    count++; current.setLength(0);
                } else {
                    current.append(line).append('\n');
                    if(current.length()>4*1024*1024)throw new IllegalArgumentException("query record exceeds 4 MiB limit");
                }
            }
            if(!current.toString().isBlank()) {
                if(count==(record==null?0:record))selected=current.toString();
                count++;
            }
        }
        if(record==null && count!=1)throw new IllegalArgumentException("multiple query records require explicit record selection");
        if(selected==null)throw new IllegalArgumentException("selected query record is absent");
        StereoMolecule input=new StereoMolecule();
        if(!new MolfileParser().parse(input,selected))throw new IllegalArgumentException("invalid MOL/SDF query record");
        // OCL's MOL reader rotates input axes to (x,-y,-z). Restore the supplied
        // frame for coordinate provenance; this proper rotation leaves stereo/distances intact.
        for(int atom=0;atom<input.getAllAtoms();atom++) {
            input.setAtomY(atom,-input.getAtomY(atom));input.setAtomZ(atom,-input.getAtomZ(atom));
        }
        return prepare(input,Map.of("input",path.toAbsolutePath().toString(),
                "inputSha256",SharedCandidatePoseBundle.sha256(path),"record",record==null?0:record,
                "recordCount",count,"coordinateSource","supplied_mol_sdf","distanceUnits","angstrom",
                "coordinateFrame","input_mol_axes"));
    }
    public static PoseQueryInput prepare(StereoMolecule input, Map<String,Object> source) {
        StereoMolecule heavy=new StereoMolecule(input);
        int n=heavy.getAllAtoms();
        double[][] original=new double[n][3];
        boolean[] hydrogen=new boolean[n];
        for(int a=0;a<n;a++) {
            original[a]=new double[]{heavy.getAtomX(a),heavy.getAtomY(a),heavy.getAtomZ(a)};
            for(double v:original[a])if(!Double.isFinite(v))throw new IllegalArgumentException("nonfinite query coordinates at atom "+a);
            hydrogen[a]=heavy.getAtomicNo(a)==1;
            if(hydrogen[a])heavy.markAtomForDeletion(a);
            // Atom maps are input annotations, not inference features. Never insert temporary map numbers.
            else heavy.setAtomMapNo(a,0,false);
        }
        int[] mapping=heavy.deleteMarkedAtomsAndBonds();
        if(mapping==null){mapping=new int[n];for(int a=0;a<n;a++)mapping[a]=a;}
        heavy.ensureHelperArrays(Molecule.cHelperCIP);
        int count=heavy.getAtoms();
        for(int a=0;a<n;a++) {
            if(hydrogen[a])continue;
            int b=mapping[a];
            if(b<0 || b>=count || heavy.getAtomicNo(b)!=input.getAtomicNo(a)
                    || heavy.getAtomX(b)!=original[a][0] || heavy.getAtomY(b)!=original[a][1]
                    || heavy.getAtomZ(b)!=original[a][2])
                throw new IllegalArgumentException("inconsistent heavy-atom coordinate mapping at atom "+a);
        }
        var batch=new DeepSpaceTensorBatchBuilder(new OCLDeepSpaceFeaturizer()).build(List.of(heavy));
        float[] distances=new float[32*32];
        boolean nonzero=false;
        for(int i=0;i<count;i++)for(int j=0;j<count;j++) {
            double dx=heavy.getAtomX(i)-heavy.getAtomX(j),dy=heavy.getAtomY(i)-heavy.getAtomY(j),dz=heavy.getAtomZ(i)-heavy.getAtomZ(j);
            float d=(float)Math.sqrt(dx*dx+dy*dy+dz*dz);
            if(!Float.isFinite(d))throw new IllegalArgumentException("query coordinate distances overflow float32");
            distances[i*32+j]=d; nonzero|=d>0;
        }
        if(!nonzero)throw new IllegalArgumentException("missing or collapsed query coordinates");
        Map<String,Object> provenance=new LinkedHashMap<>(source);
        provenance.put("heavyAtomCount",count);
        provenance.put("sourceToHeavyAtom",mapping);
        provenance.put("explicitHydrogensRemoved",n-count);
        provenance.put("atomMapsCleared",true);
        provenance.put("queryIdcode",new Canonizer(heavy).getIDCode());
        provenance.put("coordinates",java.util.stream.IntStream.range(0,count)
                .mapToObj(a->List.of(heavy.getAtomX(a),heavy.getAtomY(a),heavy.getAtomZ(a))).toList());
        return new PoseQueryInput(heavy,batch,distances,provenance);
    }
}
