package com.idorsia.research.chem.hyperspace.services;

import com.actelion.research.chem.StereoMolecule;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.IOException;
import java.io.StringWriter;
import java.util.List;

public abstract class HSService {
    public static abstract class Config{
//        public String serializeToJSON() throws IOException {
//            return new ObjectMapper().writeValueAsString(this);
//        }
//        public Config deserializeFromJSON(String json) throws IOException {
//            return (new ObjectMapper()).reader().readValue(json,this.getClass());
//        }
    };

    public abstract void runService(Config conf);


}
