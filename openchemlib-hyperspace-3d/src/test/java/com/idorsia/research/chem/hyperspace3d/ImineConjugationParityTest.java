package com.idorsia.research.chem.hyperspace3d;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.idorsia.research.chem.hyperspace3d.feature.OCLDeepSpaceFeaturizer;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/** Independent Python/RDKit expectations for conjugated N/O donors and saturated controls. */
class ImineConjugationParityTest {
    @Test void graphTensorsMatchRdkitForImineDonorsAndSaturatedControls()throws Exception {
        var input=getClass().getResourceAsStream("/com/idorsia/research/chem/hyperspace3d/imine-conjugation-v2-features.json.gz");assertNotNull(input);
        JsonNode fixture;try(var gzip=new GZIPInputStream(input)){fixture=new ObjectMapper().readTree(gzip);}
        assertEquals(1,fixture.path("formatVersion").asInt());assertEquals("deepspace7-v2",fixture.path("graphFeatureContract").asText());assertEquals(50,fixture.path("atomFeatureDim").asInt());
        var featurizer=new OCLDeepSpaceFeaturizer();
        for(var row:fixture.path("records")) {
            var mol=new StereoMolecule();new SmilesParser().parse(mol,row.path("smiles").asText());var actual=featurizer.featurize(mol);assertTrue(actual.accepted(),row.path("id").asText());
            String name=row.path("smiles").asText();assertValues(row.path("atomX"),actual.atomFeatures(),name+" atom");assertValues(row.path("pairX"),actual.pairFeatures(),name+" pair");
            for(int i=0;i<32;i++)assertEquals(row.path("atomMask").get(i).asBoolean(),actual.atomMask()[i],name+" mask "+i);
        }
    }
    private void assertValues(JsonNode expected,float[] actual,String message){assertEquals(expected.size(),actual.length);for(int i=0;i<actual.length;i++)assertEquals(expected.get(i).floatValue(),actual[i],0f,message+" feature "+i);}
}
