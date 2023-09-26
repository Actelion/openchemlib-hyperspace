package com.idorsia.research.chem.hyperspace.importer;

import com.actelion.research.chem.io.DWARFileParser;

import java.util.*;
import java.util.stream.Collectors;

public class ImporterTool {

    /**
     * virtual means that the assembled synthons will still have non-assembled "out"-connectors.
     * virtual= false means that a concrete library
     *
     * @param paths
     *
     * @return
     */
    public static ImporterSynthonReaction importSynthonReaction( List<String> paths , String idfield_name) throws Exception {

        List<ImporterSynthonSet> ssets = new ArrayList<>();
        for(String pi : paths) {
            ImporterSynthonSet iss = importSynthonSet(pi,idfield_name);
            ssets.add(iss);
        }

        ImporterSynthonReaction isr = new ImporterSynthonReaction(ssets);
        isr.checkValidity();

        System.out.println("Import successful!");
        return isr;
    }

    public static ImporterSynthonSet importSynthonSet(String path, String idfield_name) throws Exception {
        DWARFileParser in = new DWARFileParser(path);
        int si_synthon = in.getSpecialFieldIndex("Synthon");
        int si_bb      = in.getSpecialFieldIndex("Building Block");
        int i_id       = in.getFieldIndex(idfield_name);

        if(si_synthon<0 || si_bb < 0 || i_id < 0) {
            System.out.println("[INFO] parsing error -> skip file: "+path);
            return null;
        }

        List<ImporterSynthon> synthons_unfiltered = new ArrayList<>();
        while(in.next()) {
            String idc_synthon = in.getSpecialFieldData(si_synthon);
            String idc_bb      = in.getSpecialFieldData(si_bb);
            String synthon_id  = in.getFieldData(i_id);
            synthons_unfiltered.add(new ImporterSynthon(synthon_id,idc_synthon,idc_bb));
        }

        // now determine consensus synthon count..
        Map<List<Long>,List<ImporterSynthon>> sorted_by_cc = new HashMap<>();

        for(ImporterSynthon is : synthons_unfiltered) {
            List<Long> ci = Arrays.stream(is.getConnectorCounts()).asLongStream().boxed().collect(Collectors.toList());
            if(!sorted_by_cc.containsKey(ci)) {sorted_by_cc.put(ci,new ArrayList<>());}
            sorted_by_cc.get(ci).add(is);
        }

        int maxSynthons = 0;
        List<ImporterSynthon> synthons_filtered = null;
        for(List<Long> ki : sorted_by_cc.keySet()) {
            if(sorted_by_cc.get(ki).size() > maxSynthons) {
                synthons_filtered = sorted_by_cc.get(ki);
                maxSynthons = synthons_filtered.size();
            }
        }

        if(synthons_filtered == null) {
            throw new Exception("Did not find synthons..");
        }

        if(Arrays.stream(synthons_filtered.iterator().next().getConnectorCounts()).asLongStream().anyMatch(xi -> xi>1 ) ) {
            throw new Exception("More than one connector per type..");
        }

        return new ImporterSynthonSet(path, new ArrayList<>(synthons_filtered));
    }

    public static MultiStepSynthonReaction combineTwoSynthonReactions(MultiStepSynthonReaction ra, MultiStepSynthonReaction rb) throws Exception {
        return CombinedMultiStepSynthonReaction.combineTwoMultiStepReactions(ra,rb);
    }



}
