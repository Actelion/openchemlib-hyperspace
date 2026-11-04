package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.*;
import com.fasterxml.jackson.databind.*;
import com.idorsia.research.chem.hyperspace3d.feature.*;
import com.idorsia.research.chem.hyperspace3d.model.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class PoseDeploymentTest {
    @TempDir Path temp;
    @Test void coordinatesAndHydrogenMapping()throws Exception {
        StereoMolecule m=new StereoMolecule();new SmilesParser().parse(m,"CCCCCC");
        for(int a=0;a<m.getAllAtoms();a++){m.setAtomX(a,a);m.setAtomY(a,a%2);m.setAtomMapNo(a,a+1,false);}
        var planar=PoseQueryInput.prepare(m,Map.of());
        assertEquals(0,m.getAtomZ(0)); // legitimate planar coordinates are accepted
        int h=m.addAtom(1);m.addBond(0,h,Molecule.cBondTypeSingle);m.setAtomY(h,-1);
        var explicit=PoseQueryInput.prepare(m,Map.of());
        assertArrayEquals(planar.distances(),explicit.distances());
        assertArrayEquals(planar.batch().atomFeatures(),explicit.batch().atomFeatures());
        for(int a=0;a<m.getAllAtoms();a++){double x=m.getAtomX(a),y=m.getAtomY(a);m.setAtomX(a,10-y);m.setAtomY(a,x-3);m.setAtomZ(a,5);}
        assertArrayEquals(explicit.distances(),PoseQueryInput.prepare(m,Map.of()).distances(),1e-6f);
        m.setAtomX(0,Double.NaN);
        assertThrows(IllegalArgumentException.class,()->PoseQueryInput.prepare(m,Map.of()));
        StereoMolecule absent=new StereoMolecule();new SmilesParser().parse(absent,"CCCCCC");
        for(int a=0;a<absent.getAllAtoms();a++){absent.setAtomX(a,0);absent.setAtomY(a,0);}
        assertThrows(IllegalArgumentException.class,()->PoseQueryInput.prepare(absent,Map.of()));
    }
    @Test void publicFixtureGraphAndDistanceParityWithoutModels()throws Exception {
        Path golden=Path.of("src/test/resources/pose-golden");
        try(var in=new GZIPInputStream(Files.newInputStream(golden.resolve("golden.json.gz")))) {
            for(JsonNode r:new ObjectMapper().readTree(in).path("records")) {
                var input=PoseQueryInput.read(golden.resolve(r.path("file").asText()),null);
                assertArrayEquals(floats(r.path("atomX")),input.batch().atomFeatures(),0f);
                assertArrayEquals(floats(r.path("pairX")),input.batch().pairFeatures(),0f);
                assertArrayEquals(floats(r.path("distances")),input.distances(),2e-6f);
            }
        }
    }
    @Test void realModelGoldenParity()throws Exception {
        String root=System.getProperty("pose.deployment");
        Assumptions.assumeTrue(root!=null,"set -Dpose.deployment to the generated deployment directory");
        Path dir=Path.of(root),golden=dir.resolve("golden");
        var graph=DeepSpaceModelBundle.load(Path.of("model-bundles/deepspace7-v1"));
        var bundle=SharedCandidatePoseBundle.load(dir.resolve("shared-candidate-seed17-verified"),graph);
        JsonNode fixture;
        try(var in=new GZIPInputStream(Files.newInputStream(golden.resolve("golden.json.gz")))) {
            fixture=new ObjectMapper().readTree(in);
        }
        String originalHash=graph.manifest().predictorCheckpointSha256;
        graph.manifest().predictorCheckpointSha256="0".repeat(64);
        assertThrows(DeepSpaceInferenceException.class,()->SharedCandidatePoseBundle.load(bundle.directory(),graph));
        graph.manifest().predictorCheckpointSha256=originalHash;
        var device=DeepSpaceOnnxEnvironment.Device.valueOf(System.getProperty("pose.device","CPU"));
        try(var runtime=new SharedCandidatePoseRuntime(new DeepSpaceOnnxEnvironment(device,0),bundle,.5)) {
            List<float[]> embeddings=new ArrayList<>();
            for(JsonNode record:fixture.path("records")) {
                var input=PoseQueryInput.read(golden.resolve(record.path("file").asText()),null);
                assertArrayEquals(floats(record.path("atomX")),input.batch().atomFeatures(),0f);
                assertArrayEquals(floats(record.path("pairX")),input.batch().pairFeatures(),0f);
                assertArrayEquals(floats(record.path("distances")),input.distances(),2e-6f);
                var embedding=runtime.encode(input.batch(),input.distances())[0];embeddings.add(embedding);
                assertArrayEquals(floats(record.path("embedding")),embedding,2e-5f);
                var reordered=new StereoMolecule(input.molecule());reordered.swapAtoms(0,reordered.getAllAtoms()-1);
                var reorderedInput=PoseQueryInput.prepare(reordered,Map.of());
                assertArrayEquals(embedding,runtime.encode(reorderedInput.batch(),reorderedInput.distances())[0],2e-5f);
                var candidate=floats(record.path("candidate"));
                var scores=runtime.compare(embedding,new float[][]{candidate})[0];
                assertArrayEquals(floats(record.path("scores")),scores,2e-5f);
                assertEquals(.5*(scores[1]+scores[2]),scores[0],1e-6);
                assertFalse(Arrays.equals(scores,runtime.compare(candidate,new float[][]{embedding})[0]),"directed comparator must preserve query/candidate roles");
                float[][] many=new float[7][];Arrays.fill(many,candidate);
                for(float[] row:runtime.compare(embedding,many))assertArrayEquals(scores,row,2e-5f);
            }
            float difference=0;
            for(int i=0;i<128;i++)difference=Math.max(difference,Math.abs(embeddings.get(0)[i]-embeddings.get(1)[i]));
            assertTrue(difference>1e-3,"distinct hexane conformations must affect encoding beyond numerical noise");
        }
    }
    @Test void requiresRecordSelectionAndRejectsMalformedInput()throws Exception {
        Path fixture=Path.of("src/test/resources/pose-golden/query-0.mol");
        String mol=Files.readString(fixture);
        var single=PoseQueryInput.read(fixture,null);
        assertEquals(-.0046,single.molecule().getAtomY(0),1e-7);
        Path multiple=temp.resolve("poses.sdf");Files.writeString(multiple,mol+"\n$$$$\n"+mol+"\n$$$$\n");
        assertThrows(IllegalArgumentException.class,()->PoseQueryInput.read(multiple,null));
        assertEquals(0,PoseQueryInput.read(multiple,0).provenance().get("record"));
        assertEquals(1,PoseQueryInput.read(multiple,1).provenance().get("record"));
        assertThrows(IllegalArgumentException.class,()->PoseQueryInput.read(multiple,2));
        Path bad=temp.resolve("missing.mol");Files.writeString(bad,"no coordinates");
        assertThrows(IllegalArgumentException.class,()->PoseQueryInput.read(bad,null));
    }
    private static float[] floats(JsonNode array) {
        float[] values=new float[array.size()];for(int i=0;i<values.length;i++)values[i]=(float)array.get(i).asDouble();return values;
    }
}
