package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.JsonNode;
import java.time.Instant;
import java.util.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Shared browser/agent operations with immutable graph revisions and optimistic writes. */
public final class QueryService {
    public record PublishResult(long workspaceRevision,WorkspaceStore.QueryRevision query,
                                QueryChemistry.Inspection inspection,boolean oldReferencesInvalidated,boolean poseDetached){}
    private final WorkspaceStore store;
    public QueryService(WorkspaceStore store){this.store=store;}
    public PublishResult create(long expected,String actor,String format,String structure,Boolean fragment,
                                String coordinates,String name,String branchFrom)throws Exception {
        QueryChemistry.Document document=QueryChemistry.importStructure(format,structure,fragment,coordinates);
        return publish(expected,actor,document,null,name,branchFrom,"query.create");
    }
    public QueryChemistry.Inspection inspect(String queryId){return QueryChemistry.inspect(queryId,revision(queryId).document());}
    public WorkspaceStore.QueryRevision revision(String id) {
        var query=store.snapshot().queries.get(id);
        if(query==null)throw new WorkbenchException("query_not_found",404,"Unknown query revision: "+id);
        return query;
    }
    public PublishResult edit(long expected,String actor,String baseId,List<QueryChemistry.Edit> edits,String name) {
        var original=revision(baseId);
        QueryChemistry.Document changed=QueryChemistry.edit(baseId,original.document(),edits);
        return publish(expected,actor,changed,null,name,baseId,"query.edit");
    }
    public PublishResult branch(long expected,String actor,String queryId,String name) {
        var source=revision(queryId);var pose=source.pose()==null||source.pose().isNull()?null:source.pose().deepCopy();
        return publish(expected,actor,source.document(),pose,name==null?"Branch · "+source.name():name,source.id(),"query.branch");
    }
    public ResultService.Page<WorkspaceStore.QueryRevision> history(int offset,int limit) {
        require(offset>=0&&limit>=1&&limit<=100,"History pages contain 1–100 revisions");var rows=new ArrayList<>(store.snapshot().queries.values());Collections.reverse(rows);
        int start=Math.min(offset,rows.size());return new ResultService.Page<>(List.copyOf(rows.subList(start,Math.min(start+limit,rows.size()))),offset,limit,rows.size(),"Immutable query revisions; branching preserves the original run and pose provenance.");
    }
    public PublishResult attachPose(long expected,String actor,QueryChemistry.Document document,JsonNode pose,String name) {
        require(pose!=null&&!pose.isNull(),"Pose provenance is required");
        require(!document.fragment(),"A supplied pose must be a concrete molecule");
        return publish(expected,actor,document,pose.deepCopy(),name,store.snapshot().currentQueryId,"query.attach_pose");
    }
    private PublishResult publish(long expected,String actor,QueryChemistry.Document document,JsonNode pose,String name,String parent,String operation) {
        return store.mutate(expected,actor,operation,state->{
            if(parent!=null)require(state.queries.containsKey(parent),"Parent query revision is absent");
            // Typed edits always edit the visible revision. Branching from history is a separate create operation.
            if(operation.equals("query.edit")&&!Objects.equals(parent,state.currentQueryId))
                throw new WorkbenchException("stale_query",409,"The scientist has selected a newer query; inspect it before applying edits");
            WorkspaceStore.pushUndo(state);
            var previous=state.currentQueryId==null?null:state.queries.get(state.currentQueryId);
            var query=new WorkspaceStore.QueryRevision("query-"+UUID.randomUUID(),parent,
                    name==null||name.isBlank()?"Untitled query":name,document,pose,actor,Instant.now().toString());
            require(query.name().length()<=200,"Query name is too long");
            state.queries.put(query.id(),query);state.currentQueryId=query.id();state.selection=WorkspaceStore.Selection.empty(query.id());
            return new PublishResult(state.revision+1,query,QueryChemistry.inspect(query.id(),document),true,
                    previous!=null&&previous.pose()!=null&&!previous.pose().isNull()&&pose==null);
        });
    }
    public WorkspaceStore.Selection select(long expected,String actor,WorkspaceStore.Selection selection) {
        require(selection!=null,"Selection required");
        require(selection.hits().size()<=10000,"At most 10,000 hit occurrences may be selected");
        return store.mutate(expected,actor,"selection.change",state->{
            require(Objects.equals(selection.queryRevisionId(),state.currentQueryId),"Selection must refer to the current query revision");
            if(selection.queryRevisionId()!=null) {
                var query=state.queries.get(selection.queryRevisionId());
                var inspection=QueryChemistry.inspect(query.id(),query.document());
                Set<String> atoms=new HashSet<>(inspection.atoms().stream().map(QueryChemistry.Atom::ref).toList());
                Set<String> bonds=new HashSet<>(inspection.bonds().stream().map(QueryChemistry.Bond::ref).toList());
                if(!atoms.containsAll(selection.atoms())||!bonds.containsAll(selection.bonds()))
                    throw new WorkbenchException("stale_reference",409,"Selected atom/bond references do not belong to the visible revision");
            } else require(selection.atoms().isEmpty()&&selection.bonds().isEmpty(),"No query is selected");
            WorkspaceStore.pushUndo(state);state.selection=selection;return selection;
        });
    }
    public PublishResult editSelection(long expected,String actor,QueryChemistry.Operation operation,List<Integer> elements,
                                      String stateValue,Integer min,Integer max,Boolean enabled) {
        return editSelection(expected,actor,operation,elements,stateValue,min,max,enabled,null,false);
    }
    public PublishResult editSelection(long expected,String actor,QueryChemistry.Operation operation,List<Integer> elements,
                                      String stateValue,Integer min,Integer max,Boolean enabled,List<String> types,Boolean exclude) {
        var selection=store.snapshot().selection;
        require(selection.queryRevisionId()!=null,"Select a query first");
        boolean atom=operation.name().startsWith("ATOM_");
        var edit=new QueryChemistry.Edit(operation,atom?selection.atoms():selection.bonds(),elements,exclude,stateValue,types,enabled,min,max);
        return edit(expected,actor,selection.queryRevisionId(),List.of(edit),"Edited selection");
    }
    public WorkspaceStore.State undo(long expected,String actor) {
        store.mutate(expected,actor,"workspace.undo",state->{
            require(!state.undo.isEmpty(),"No query/selection edit can be undone");
            var frame=state.undo.remove(state.undo.size()-1);state.currentQueryId=frame.queryId();state.selection=frame.selection();return null;
        });
        return store.snapshot();
    }
}
