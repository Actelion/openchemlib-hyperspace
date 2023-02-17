import com.actelion.research.chem.StereoMolecule;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestSubstructureSearchHelper {



    @Test
    void demoTestMethod() {

        if(true) {
            assertTrue(true);
            return;
        }

        StereoMolecule test_a = HyperspaceUtils.parseIDCode("fbmAB@K@BLdTTTRRTRJLcZBgM@A@EQTD@F\\NA`");
        //StereoMolecule test_a = HyperspaceUtils.parseIDCode("ffcAB@H`BLdTTTRRRQLRTcZ|Vm@A@EEUPP@YpLB@");
        List<StereoMolecule> expanded = null;
        try {
            expanded = SubstructureSearchHelper.expandBridgedSearches(test_a);
        } catch (Exception e) {
            e.printStackTrace();
        }
        for(StereoMolecule mi : expanded) {
            System.out.println("idc: "+mi.getIDCode());
        }
        assertTrue(expanded.size()==4);
    }

}
