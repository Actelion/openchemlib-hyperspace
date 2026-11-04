package com.idorsia.research.chem.hyperspace.workbench;

import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class CollaborationHistoryTest {
    @TempDir Path temp;
    @Test void runBranchesUndoAndNamedWorkspacePreserveOriginalQueriesAndSelections()throws Exception {
        String queryId,runId;try(var app=new WorkbenchApplication(temp)){
            app.invoke("workspace.rename",WorkspaceStore.JSON.valueToTree(Map.of("name","Scientist and agent")),0,"scientist");
            var first=app.queries.create(1,"scientist","smiles","CCO",true,null,"First query",null);queryId=first.query().id();
            var run=app.results.create(2,"scientist","substructure",first.query(),List.of(),Map.of("limit",20),"Saved test scope");runId=run.id();
            app.queries.create(3,"agent:mcp","smiles","c1ccccc1",true,null,"Newer query",null);
            var branch=(QueryService.PublishResult)app.invoke("run.branch_query",WorkspaceStore.JSON.valueToTree(Map.of("runId",runId)),4,"scientist");
            assertEquals(queryId,branch.query().parentId());assertEquals(first.query().document(),branch.query().document());assertFalse(branch.poseDetached());
            assertEquals(queryId,app.results.getRun(runId).queryId());assertEquals(3,app.queries.history(0,100).total());assertEquals(1,app.queries.history(1,1).items().size());
            assertThrows(WorkbenchException.class,()->app.queries.branch(4,"agent:mcp",queryId,"Stale"));
            app.queries.undo(5,"scientist");assertEquals("Newer query",app.queries.revision(app.store.snapshot().currentQueryId).name());assertEquals(queryId,app.results.getRun(runId).queryId());
        }
        try(var app=new WorkbenchApplication(temp)){assertEquals("Scientist and agent",app.store.snapshot().name);assertEquals(queryId,app.results.getRun(runId).queryId());}
    }
    @Test void oldSavedRunsRemainNavigableBeyondTheRecentSummaryWindow()throws Exception {
        try(var app=new WorkbenchApplication(temp)){
            var q=app.queries.create(0,"scientist","smiles","CCO",true,null,"History query",null).query();String first=null;
            for(int i=0;i<105;i++){var r=app.results.create(app.store.snapshot().revision,"scientist","substructure",q,List.of(),Map.of(),"History boundary");if(i==0)first=r.id();}
            assertEquals(100,((List<?>)app.state().get("runs")).size());assertEquals(105,app.state().get("runCount"));
            var page=app.history("runs",100,20);assertEquals(5,page.items().size());assertEquals(first,page.items().get(4).path("id").asText());
            app.invoke("workspace.navigate",WorkspaceStore.JSON.valueToTree(Map.of("view","screening","entityId",first)),app.store.snapshot().revision,"agent:mcp");
            String selected=first;assertTrue(((List<?>)app.state().get("runs")).stream().anyMatch(n->((com.fasterxml.jackson.databind.JsonNode)n).path("id").asText().equals(selected)));
            assertThrows(WorkbenchException.class,()->app.history("runs",0,101));assertThrows(WorkbenchException.class,()->app.history("unknown",0,20));
        }
    }
    @Test void disappearedAssetsFailCleanlyInsteadOfRemainingLoading()throws Exception {
        try(var app=new WorkbenchApplication(temp.resolve("workspace"))){Path raw=temp.resolve("public.rawspace");com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpaceIO.write(ToySetup.raw(),raw);
            var detection=app.assets.inspect(raw,null);assertEquals(AssetService.Kind.RAW_SPACE,detection.kind());var asset=app.assets.register(0,"scientist",detection,"Raw",null,null,null,"v1","Toy");
            Files.delete(raw);assertThrows(Exception.class,()->app.assets.load(1,"scientist",asset.id(),null));assertEquals("failed",app.assets.get(asset.id()).status());
        }
    }
    @Test void editingNotesAndMembersRetainsCollectionIdAndCreationTime()throws Exception {
        try(var app=new WorkbenchApplication(temp)){
            var q=app.queries.create(0,"scientist","smiles","CCCCCCO",false,null,"Concrete",null).query();
            var run=app.results.create(1,"scientist","test",q,List.of(),Map.of(),"Retained example");var m=QueryChemistry.parse(q.document());
            var a=ResultService.makeHit(run.id(),"occurrence-a","asset",null,null,List.of(),m,null,false,Map.of("phesa_total",.7),Map.of());var b=ResultService.makeHit(run.id(),"occurrence-b","asset",null,null,List.of(),m,null,false,Map.of("phesa_total",.6),Map.of());
            app.results.complete(run.id(),List.of(a,b),Map.of());var original=app.results.saveCollection(2,"scientist",null,"Candidates","Initial",List.of(a.id(),b.id()));
            var edited=app.results.saveCollection(3,"agent:mcp",original.id(),"Reviewed candidates","Updated notes",List.of(b.id()));assertEquals(original.id(),edited.id());assertEquals(original.createdAt(),edited.createdAt());assertEquals(List.of(b.id()),app.results.collectionHits(edited.id()).stream().map(ResultService.Hit::id).toList());
            assertEquals(0,app.results.page(run.id(),0,20,"phesa_total",true,new ResultService.Filters(10000.,null,null,null)).total());
            assertThrows(WorkbenchException.class,()->app.results.saveCollection(3,"scientist",edited.id(),"Stale","Do not overwrite",List.of(a.id())));
        }
    }
}
