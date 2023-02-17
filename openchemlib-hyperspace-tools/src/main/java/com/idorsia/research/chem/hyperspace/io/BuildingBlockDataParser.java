package com.idorsia.research.chem.hyperspace.io;

import com.actelion.research.chem.io.DWARFileParser;
import com.idorsia.research.chem.hyperspace.BuildingBlockData;

import java.io.*;

public class BuildingBlockDataParser {

    /**
     * Arguments: 1. path to input folder 2. output filename
     *
     * @param args
     */
    public static void main(String args[]) {
        BuildingBlockData data = new BuildingBlockData();
        String folder = args[0];
        String output = args[1];
        File ff       = new File(folder);
        for(File fi : ff.listFiles() ) {
            DWARFileParser p = new DWARFileParser(fi);

            int sidx_bb_structure         = p.getSpecialFieldIndex("Building Block");
            int sidx_bb_synthon_structure = p.getSpecialFieldIndex("Synthon");
            int idx_loc                  = p.getFieldIndex("Location");
            int idx_barcode              = p.getFieldIndex("Barcode");
            int idx_type                 = p.getFieldIndex("Type");

            while(p.next()) {
                String structure_bb = p.getSpecialFieldData( sidx_bb_structure );
                String synthon_structure = p.getSpecialFieldData(sidx_bb_synthon_structure);
                String loc           = p.getFieldData(idx_loc);
                String barcode       = p.getFieldData(idx_barcode);
                String type          = "";
                if(idx_type>=0){type = p.getFieldData(idx_type);}
                data.putBuildingBlockData(structure_bb,BuildingBlockData.BB_STRUCTURE,structure_bb);
                data.putBuildingBlockData(structure_bb,BuildingBlockData.BB_STRUCTURE_SYNTHON,synthon_structure);
                data.putBuildingBlockData(structure_bb,BuildingBlockData.BB_LOCATION,loc);
                data.putBuildingBlockData(structure_bb,BuildingBlockData.BB_BARCODE,barcode);
                data.putBuildingBlockData(structure_bb,BuildingBlockData.BB_SYNTHON_TYPE,type);
            }
        }
        // write output file:
        FileOutputStream f_out = null;
        try {
            f_out = new FileOutputStream(new File(output));
            ObjectOutputStream out = new ObjectOutputStream(f_out);
            out.writeObject(data);
            out.flush();
            out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static BuildingBlockData loadBuildingBlockData(File fi) {
        BuildingBlockData data = new BuildingBlockData();
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream(fi));
            data = (BuildingBlockData) in.readObject();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        return data;
    }

}
