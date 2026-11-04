package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.node.JsonNodeFactory;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import static org.junit.jupiter.api.Assertions.*;

class QueryWorkspaceTest {
    @TempDir Path temporary;
    private QueryChemistry.Edit atom(String ref,List<Integer> elements){return new QueryChemistry.Edit(QueryChemistry.Operation.ATOM_ELEMENTS,List.of(ref),elements,false,null,null,null,null,null);}
    @Test void actualAtomAlternativesAndBridgeSemanticsRoundTrip()throws Exception {
        var graph=QueryChemistry.importStructure("smiles","CC",true,null);
        var changed=QueryChemistry.edit("r1",graph,List.of(atom("r1:a0",List.of(7,8))));
        var matches=QueryChemistry.test(changed,List.of("CN","CO","CC"));
        assertEquals(List.of(true,true,false),matches.stream().map(QueryChemistry.ExampleResult::matches).toList());
        var linker=QueryChemistry.importStructure("smiles","CO",true,null);
        var bridged=QueryChemistry.edit("r2",linker,List.of(new QueryChemistry.Edit(QueryChemistry.Operation.BOND_BRIDGE,
                List.of("r2:b0"),null,null,null,null,null,1,2)));
        var info=QueryChemistry.inspect("r3",bridged);assertTrue(info.bonds().get(0).bridge());
        assertEquals(1,info.bonds().get(0).minAtoms());assertEquals(2,info.bonds().get(0).maxAtoms());
        assertEquals(List.of(false,true,true),QueryChemistry.test(bridged,List.of("CO","CCO","CCCO")).stream().map(QueryChemistry.ExampleResult::matches).toList());
        assertThrows(WorkbenchException.class,()->QueryChemistry.edit("r2",linker,List.of(new QueryChemistry.Edit(
                QueryChemistry.Operation.BOND_BRIDGE,List.of("r2:b0"),null,null,null,null,null,0,0))));
        var copied=QueryChemistry.importStructure("idcode",bridged.idcode(),null,bridged.coordinates());
        assertEquals(bridged.idcode(),copied.idcode());assertTrue(copied.fragment());
        // Consumed by the real OCL JS round-trip verification, not a handwritten encoded query.
        Files.createDirectories(Path.of("target"));
        WorkspaceStore.JSON.writeValue(Path.of("target/query-roundtrip.json").toFile(),List.of(graph,changed,bridged));
    }
    @Test void staleFeatureReferencesAndFailedMutationsCannotChangeState()throws Exception {
        try(var store=new WorkspaceStore(temporary)) {
            var service=new QueryService(store);
            var first=service.create(0,"scientist","smiles","CC",true,null,"first",null);
            var second=service.edit(1,"agent:local",first.query().id(),List.of(atom(first.inspection().atoms().get(0).ref(),List.of(7,8))),"edited");
            assertEquals(2,second.workspaceRevision());assertTrue(second.oldReferencesInvalidated());
            assertThrows(WorkbenchException.class,()->service.edit(2,"agent:local",second.query().id(),
                    List.of(atom(first.inspection().atoms().get(0).ref(),List.of(6))),"stale"));
            assertEquals(2,store.snapshot().revision);assertEquals(second.query().id(),store.snapshot().currentQueryId);
            assertThrows(WorkbenchException.class,()->service.create(0,"scientist","smiles","CO",true,null,"conflict",null));
            assertEquals(2,store.snapshot().queries.size());
        }
    }
    @Test void poseDetachesOnGraphEditsAndUndoRestoresSelectionWithoutRemovingRevisions()throws Exception {
        try(var store=new WorkspaceStore(temporary)) {
            var service=new QueryService(store);
            var concrete=QueryChemistry.importStructure("smiles","CCO",false,null);
            var pose=JsonNodeFactory.instance.objectNode().put("input","query.sdf").put("record",0);
            var original=service.attachPose(0,"scientist",concrete,pose,"pose");
            var selected=new WorkspaceStore.Selection(original.query().id(),List.of(original.inspection().atoms().get(0).ref()),List.of(),List.of());
            service.select(1,"scientist",selected);
            var edited=service.editSelection(2,"agent:local",QueryChemistry.Operation.ATOM_ELEMENTS,List.of(6,7),null,null,null,null);
            assertTrue(edited.poseDetached());assertNull(edited.query().pose());assertEquals(0,store.snapshot().selection.atoms().size());
            var undone=service.undo(3,"scientist");assertEquals(original.query().id(),undone.currentQueryId);
            assertEquals(selected,undone.selection);assertEquals(2,undone.queries.size());
            assertEquals("workspace.undo",undone.activity.get(3).operation());
        }
        try(var reopened=new WorkspaceStore(temporary)) {
            assertEquals(4,reopened.snapshot().revision);assertEquals(2,reopened.snapshot().queries.size());
            assertFalse(reopened.snapshot().selection.atoms().isEmpty());
        }
    }
    @Test void onlyOneOwnerAndSnapshotsDoNotMutateLiveState()throws Exception {
        String id;
        try(var store=new WorkspaceStore(temporary)) {
            id=store.snapshot().workspaceId;
            assertThrows(java.io.IOException.class,()->new WorkspaceStore(temporary));
            var snapshot=store.snapshot();snapshot.name="uncommitted change";assertEquals("Hyperspace workbench",store.snapshot().name);
        }
        try(var reopened=new WorkspaceStore(temporary)){assertEquals(id,reopened.snapshot().workspaceId);}
    }
}
