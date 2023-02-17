package com.idorsia.research.chem.hyperspace;

import com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI;
import com.idorsia.research.chem.hyperspace.gui.HyperspaceServer;
import com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser2;

public class HyperspaceLauncher {

    /**
     * Entrypoint for all Hyperspace functionality (SpaceCreation, GUI, Server, etc..)
     *
     * First parameter can be:
     * CREATESPACE (then it uses the SynthonSpaceParser2 )
     * GUI (then it uses the HyperspaceSearchGUI )
     * SERVER (then it uses the HyperspaceServer)
     *
     */
    public static void main(String args[]) {
        String args2[] = new String[args.length-1];
        for(int zi=1;zi<args.length;zi++) {args2[zi-1]=args[zi];}

        if(args[0].equalsIgnoreCase("CREATESPACE")) {
            SynthonSpaceParser2.main(args2);
        }
        else if(args[0].equalsIgnoreCase("GUI")) {
            HyperspaceSearchGUI.main(args2);
        }
        else if(args[0].equalsIgnoreCase("SERVER")) {
            HyperspaceServer.main(args2);
        }
    }

}
