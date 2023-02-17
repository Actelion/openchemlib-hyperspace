package com.idorsia.research.chem.hyperspace.gui;

import com.idorsia.research.chem.hyperspace.HyperspaceIOUtils;
import com.idorsia.research.chem.hyperspace.HyperspaceUtils;
import com.idorsia.research.chem.hyperspace.gui.search.AbstractSearchProvider;
import org.json.JSONObject;

import java.util.HashMap;
import java.util.Map;

public class HyperspaceServerEngine {

    Map<String, AbstractSearchProvider> searchProvider = new HashMap<>();


    public HyperspaceServerEngine() {

    }

    public void loadServices(String config_file) {
        String si = HyperspaceUtils.readFileIntoString(config_file);

        //JSONObject jo = new JSONObject(si);
        //HyperspaceInit.loadSearchProviderInitFile(null,si);



    }




}
