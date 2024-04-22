package com.idorsia.research.hyperspace.server;

import com.idorsia.research.chem.hyperspace.service.SubstructureSearchTask;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.messaging.handler.annotation.MessageMapping;
import org.springframework.messaging.handler.annotation.SendTo;
import org.springframework.messaging.simp.SimpMessagingTemplate;

public class Hyperspace3Controller {

    private final HyperspaceComputationService computationService;
    private final SimpMessagingTemplate messagingTemplate;

    @Autowired
    public Hyperspace3Controller(HyperspaceComputationService computationService, SimpMessagingTemplate messagingTemplate) {
        this.computationService = computationService;
        this.messagingTemplate = messagingTemplate;
    }

    @MessageMapping("/runSubstructureSearch")
    public void processComputationTask(SubstructureSearchTask task) {
        // Perform computation...
        computationService.executeComputation(task,false);
        //String result = doComputation(task);
        //messagingTemplate.convertAndSend("/sss/results", result);
        //doComputation(task)
        // Return value will be sent to the subscriber in "/topic/results" as a message
        //return result;
    }

    @MessageMapping("/runSubstructureSearch_ExpandedResults")
    public void processComputationTask_Expanded(SubstructureSearchTask task) {
        // Perform computation...
        //String result = doComputation(task);
        computationService.executeComputation(task,true);
        //messagingTemplate.convertAndSend("/sss/results_expanded", result);
        // Return value will be sent to the subscriber in "/topic/results" as a message
        //return result;
    }

}
