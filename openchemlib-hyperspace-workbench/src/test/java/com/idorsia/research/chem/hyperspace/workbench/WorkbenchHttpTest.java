package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.JsonNode;
import java.net.*;
import java.net.http.*;
import java.nio.file.*;
import java.time.Duration;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class WorkbenchHttpTest {
    @TempDir Path temp;
    private final HttpClient client=HttpClient.newHttpClient();
    private HttpResponse<String> call(WorkbenchServer server,String operation,Object args,long revision)throws Exception {
        return client.send(HttpRequest.newBuilder(URI.create(server.url()+"/api/call")).header("X-Workbench-Token",server.token()).header("Content-Type","application/json")
            .POST(HttpRequest.BodyPublishers.ofString(WorkspaceStore.JSON.writeValueAsString(Map.of("operation",operation,"arguments",args,"actor","agent:test","expectedRevision",revision)))).build(),HttpResponse.BodyHandlers.ofString());
    }
    private JsonNode result(HttpResponse<String> response)throws Exception {assertEquals(200,response.statusCode(),response.body());return WorkspaceStore.JSON.readTree(response.body()).path("result");}
    @Test void authenticatedOwnerSharesAgentMutationsAndRunsAfterDisconnect()throws Exception {
        Path workspace=temp.resolve("workspace");String query;
        try(var server=new WorkbenchServer(workspace,0,temp.resolve("frontend"))) {
            assertEquals(401,client.send(HttpRequest.newBuilder(URI.create(server.url()+"/api/state")).GET().build(),HttpResponse.BodyHandlers.ofString()).statusCode());
            assertEquals(400,client.send(HttpRequest.newBuilder(URI.create(server.url()+"/api/state")).header("X-Workbench-Token",server.token()).header("Origin","https://foreign.example").GET().build(),HttpResponse.BodyHandlers.ofString()).statusCode());
            var connection=WorkspaceStore.JSON.readTree(workspace.resolve("connection.json").toFile());assertEquals(server.url(),connection.path("url").asText());
            var published=result(call(server,"query.create",Map.of("structure","CCO","format","smiles","fragment",true,"name","Shared query"),0));query=published.path("query").path("id").asText();
            var state=result(call(server,"workspace.state",Map.of(),-1));assertEquals(query,state.path("currentQuery").path("id").asText());assertEquals("agent:test",state.path("activity").get(0).path("actor").asText());
            assertEquals(409,call(server,"query.create",Map.of("structure","CO"),0).statusCode());
            assertEquals(400,call(server,"query.create",Map.of("structure","CO","unknown",true),1).statusCode());
            assertEquals(400,call(server,"query.create",Map.of("structure","CO","fragment","yes"),1).statusCode());
            String atom=published.path("inspection").path("atoms").get(0).path("ref").asText();
            result(call(server,"selection.change",Map.of("queryRevisionId",query,"atoms",List.of(atom),"bonds",List.of(),"hits",List.of()),1));
            var changed=result(call(server,"query.edit_selection",Map.of("operation","ATOM_ELEMENTS","elements",List.of(7,8)),2));
            assertTrue(changed.path("oldReferencesInvalidated").asBoolean());
            assertEquals(409,call(server,"query.edit",Map.of("queryId",changed.path("query").path("id").asText(),"edits",List.of(Map.of("operation","ATOM_ELEMENTS","refs",List.of(atom),"elements",List.of(6)))),3).statusCode());
            var job=result(call(server,"toy.prepare",Map.of(),3));String jobId=job.path("id").asText();
            // No browser/session context is retained. A separate request later discovers the completed owner job.
            long deadline=System.nanoTime()+Duration.ofSeconds(30).toNanos();JsonNode status=null;
            while(System.nanoTime()<deadline){status=result(call(server,"job.status",Map.of("jobId",jobId),-1));if(Set.of("completed","failed").contains(status.path("state").asText()))break;Thread.sleep(10);}
            assertNotNull(status);assertEquals("completed",status.path("state").asText(),status.toString());
            assertEquals(2,result(call(server,"workspace.state",Map.of(),-1)).path("assets").size());
        }
        assertFalse(Files.exists(workspace.resolve("connection.json")));
        try(var server=new WorkbenchServer(workspace,0,temp.resolve("frontend"))) {
            var state=result(call(server,"workspace.state",Map.of(),-1));assertEquals(2,state.path("assets").size());assertEquals(1,state.path("jobs").size());assertEquals("completed",state.path("jobs").get(0).path("state").asText());
            assertNotEquals(query,state.path("currentQuery").path("id").asText());
        }
    }
}
