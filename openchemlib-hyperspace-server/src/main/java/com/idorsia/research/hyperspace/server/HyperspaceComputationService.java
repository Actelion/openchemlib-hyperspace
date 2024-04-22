package com.idorsia.research.hyperspace.server;


import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.hyperspace.SimpleCombinatorialHit;
import com.actelion.research.chem.hyperspace.SimpleSynthon;
import com.idorsia.research.chem.hyperspace.CachedDescriptorProvider;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.service.SubstructureSearchTask;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.messaging.simp.SimpMessagingTemplate;
import org.springframework.scheduling.annotation.Async;
import org.springframework.scheduling.annotation.EnableAsync;
import org.springframework.stereotype.Service;

import java.util.List;
import java.util.Map;

@EnableAsync
@Service
public class HyperspaceComputationService {

    private int numThreadsPerQuery = 8;

    private final SimpMessagingTemplate messagingTemplate;
    private final ConfigurationService configurationService;

    @Autowired
    public HyperspaceComputationService(SimpMessagingTemplate messagingTemplate, ConfigurationService configurationService) {
        this.messagingTemplate = messagingTemplate;
        this.configurationService = configurationService;
    }

    @Async("computationExecutor")
    public void executeComputation(SubstructureSearchTask task, boolean returnAssembled) {
        // Your computation logic here
        doComputation(task, returnAssembled);
    }

    private void doComputation(SubstructureSearchTask task, boolean returnAssembled) {
        // TODO: implement this correctly..
        SynthonSpace space = configurationService.getLoadedSpaces().entrySet().iterator().next().getValue();

        CachedDescriptorProvider cdp = new CachedDescriptorProvider("FragFp");
        StereoMolecule qi = new StereoMolecule();
        IDCodeParser icp = new IDCodeParser();
        icp.parse(qi,task.getQueryIDCode());
        qi.ensureHelperArrays(Molecule.cHelperCIP);

        SynthonSpace.CombinatorialHitReceiver receiver = new SynthonSpace.CombinatorialHitReceiver() {
            @Override
            public void addCombinatorialHits(List<SynthonSpace.CombinatorialHit> hits, List<SynthonSpace.SplitPatternWithConnectorProximityPruningStatistics> stats) {
                if(!returnAssembled) {
                    // using the messagingTemplate
                    for(SynthonSpace.CombinatorialHit hi : hits) {
                        SimpleCombinatorialHit shit = convertToSimpleCombinatorialHit(hi);
                        messagingTemplate.convertAndSend("/sss/results_combinatorial",shit);
                    }
                }
                else {
                    // TODO: assemble hits..
                    //messagingTemplate.convertAndSend("/sss/results_expanded",);
                }
            }
        };
        try {
            SubstructureSearchHelper.run_substructure_search_streaming_01_withBridgedBondsExpansion(space,
                    cdp, qi,
                    numThreadsPerQuery,
                    true, true,
                    receiver);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static SimpleCombinatorialHit convertToSimpleCombinatorialHit(SynthonSpace.CombinatorialHit hit) {
        SimpleSynthon[][] synthons = new SimpleSynthon[hit.hit_fragments.size()][];
        int scnt = 0;
        for(Map.Entry<SynthonSpace.FragType, List<SynthonSpace.FragId>> set_i : hit.hit_fragments.entrySet()) {
            SimpleSynthon[] synthons_i = new SimpleSynthon[set_i.getValue().size()];
            String synthon_set = "SynthonSet:"+set_i.getKey().frag;
            for(int zi=0;zi<set_i.getValue().size();zi++) {
                SynthonSpace.FragId fid = set_i.getValue().get(zi);
                synthons_i[zi] = new SimpleSynthon(fid.idcode,fid.fragment_id,fid.rxn_id,synthon_set);
            }
            scnt++;
            synthons[scnt] = synthons_i;
        }
        SimpleCombinatorialHit shit = new SimpleCombinatorialHit(hit.rxn,synthons);
        return shit;
    }

}
