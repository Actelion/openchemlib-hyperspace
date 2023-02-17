package com.idorsia.research.chem.hyperspace.server;

public class HyperspaceServer {

    public static void main(String args[]) {
        String args_start_server[] = new String[]{"-t" , "32" , "--init" , "/opt/chemspaceserver/data/hyperspace_config_server_cloudp02.json" };
        com.idorsia.research.chem.hyperspace.gui.HyperspaceServer.main(args_start_server);
    }

}
