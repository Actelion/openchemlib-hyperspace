package com.idorsia.research.chem.hyperspace.descriptor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerStandard2DFactory;

public class MyDescriptorResolver {

    public static DescriptorHandler<long[], StereoMolecule> resolveDescriptorHandlerFromName(String shortName) {
        if(shortName.equals("FFP1024_plus_ffp")) {
            return new DescriptorHandlerLongFFP1024_plus("ffp");
        }
        else if(shortName.equals("FFP1024_plus_pfp")) {
            return new DescriptorHandlerLongFFP1024_plus("pfp");
        }
        else if(shortName.equals("PPC2048")) {
            return new DescriptorHandlerPPCore();
        }
        else {
            return DescriptorHandlerStandard2DFactory.getFactory().create(shortName);
        }
    }

}
