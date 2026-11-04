package com.idorsia.research.chem.hyperspace.workbench;

import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class PoseImportTest {
    @TempDir Path temp;
    @Test void uploadAndExplicitPoseRecordKeepImmutableCoordinatesAndActorHistory()throws Exception {
        String mol=Files.readString(Path.of("../openchemlib-hyperspace-3d/src/test/resources/pose-golden/query-0.mol"));
        try(var app=new WorkbenchApplication(temp)) {
            var upload=WorkspaceStore.JSON.valueToTree(app.invoke("query.upload",WorkspaceStore.JSON.valueToTree(Map.of("content",mol+"\n$$$$\n"+mol+"\n$$$$\n","format","sdf")),0,"scientist"));
            String path=upload.path("path").asText();assertTrue(Files.isRegularFile(Path.of(path)));assertEquals(1,app.store.snapshot().revision);
            assertEquals("query.upload",app.store.snapshot().activity.get(0).operation());
            var arguments=WorkspaceStore.JSON.valueToTree(Map.of("path",path,"name","Pose fixture"));
            assertThrows(IllegalArgumentException.class,()->app.invoke("query.import_pose",arguments,1,"scientist"));
            ((com.fasterxml.jackson.databind.node.ObjectNode)arguments).put("record",1);
            var result=(QueryService.PublishResult)app.invoke("query.import_pose",arguments,1,"scientist");
            var pose=result.query().pose();assertEquals(1,pose.path("record").asInt());assertEquals(2,pose.path("recordCount").asInt());assertTrue(pose.path("coordinates").size()>0);
            String snapshot=pose.path("input").asText();assertNotEquals(path,snapshot);assertEquals(Files.readString(Path.of(path)),Files.readString(Path.of(snapshot)));
            assertFalse(result.query().document().fragment());
            Files.writeString(Path.of(path),"external upload changed after import");assertTrue(Files.readString(Path.of(snapshot)).contains("$$$$"));
            var original=app.queries.revision(result.query().id());var inspected=app.queries.inspect(original.id());
            var changed=app.queries.edit(2,"agent:mcp",original.id(),List.of(new QueryChemistry.Edit(QueryChemistry.Operation.ATOM_ELEMENTS,List.of(inspected.atoms().get(0).ref()),List.of(6,7),false,null,null,null,null,null)),"Graph branch");
            assertTrue(changed.poseDetached());assertNull(changed.query().pose());assertEquals(pose,app.queries.revision(original.id()).pose());
        }
    }
}
