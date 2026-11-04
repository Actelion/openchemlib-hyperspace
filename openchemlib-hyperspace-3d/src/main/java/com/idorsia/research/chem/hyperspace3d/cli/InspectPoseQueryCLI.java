package com.idorsia.research.chem.hyperspace3d.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.feature.PoseQueryInput;
import java.nio.file.Path;
import java.util.LinkedHashMap;

public final class InspectPoseQueryCLI {
    public static void main(String[] args) throws Exception {
        if(args.length!=1 && args.length!=2)throw new IllegalArgumentException("Usage: InspectPoseQueryCLI query.sdf [zero-based-record]");
        var pose=PoseQueryInput.read(Path.of(args[0]),args.length==2?Integer.valueOf(args[1]):null);
        var summary=new LinkedHashMap<>(pose.provenance());
        summary.remove("coordinates");
        summary.put("valid",true);
        new ObjectMapper().writeValue(System.out,summary);
    }
}
