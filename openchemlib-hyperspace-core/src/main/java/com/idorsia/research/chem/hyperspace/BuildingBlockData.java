package com.idorsia.research.chem.hyperspace;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

public class BuildingBlockData implements Serializable {

    public static final String BB_STRUCTURE          = "Structure";
    public static final String BB_STRUCTURE_SYNTHON  = "SynthonStructure";
    public static final String BB_BARCODE    = "Barcode";
    public static final String BB_LOCATION   = "Location";
    public static final String BB_SYNTHON_TYPE = "SynthonType";



    /**
     * Map from BB-Id
     */
    Map<String, Map<String,String>> mData = new HashMap<>();

    public BuildingBlockData(){
    }

    public synchronized String getBuildingBlockData(String frag_id, String field_id) {
        if(!mData.containsKey(frag_id)) {
            return null;
        }
        if(!mData.get(frag_id).containsKey(field_id)) {
            return null;
        }
        return mData.get(frag_id).get(field_id);
    }

    public synchronized void putBuildingBlockData(String frag_id, String field_id, String value) {
        if(!mData.containsKey(frag_id)) {
            mData.put(frag_id,new HashMap<>());
        }
        mData.get(frag_id).put(field_id,value);
    }

}
