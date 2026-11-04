package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.*;
import java.io.IOException;
import java.nio.channels.*;
import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.function.Function;

/** One process owns the metadata and loaded-resource lifecycle for a workspace. */
public final class WorkspaceStore implements AutoCloseable {
    public static final ObjectMapper JSON=new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT);
    public record QueryRevision(String id,String parentId,String name,QueryChemistry.Document document,
                                JsonNode pose,String actor,String createdAt){}
    public record Selection(String queryRevisionId,List<String> atoms,List<String> bonds,List<String> hits){
        public Selection{atoms=List.copyOf(atoms);bonds=List.copyOf(bonds);hits=List.copyOf(hits);}
        public static Selection empty(String query){return new Selection(query,List.of(),List.of(),List.of());}
    }
    public record Activity(long revision,String actor,String operation,String at){}
    public record UndoFrame(String queryId,Selection selection){}
    public static final class State {
        public int schemaVersion=1;
        public String workspaceId="workspace-"+UUID.randomUUID(),name="Hyperspace workbench",currentQueryId;
        public long revision;
        public Map<String,QueryRevision> queries=new LinkedHashMap<>();
        public Map<String,JsonNode> libraries=new LinkedHashMap<>(),sourceSnapshots=new LinkedHashMap<>(),assets=new LinkedHashMap<>(),runs=new LinkedHashMap<>(),collections=new LinkedHashMap<>();
        public Selection selection=Selection.empty(null);
        public List<Activity> activity=new ArrayList<>();
        public List<UndoFrame> undo=new ArrayList<>();
        public String view="screening",viewEntityId;
    }
    private final Path directory;
    private final FileChannel lockChannel;
    private final FileLock lock;
    private State state;
    private boolean closed;
    public WorkspaceStore(Path directory)throws IOException {
        this.directory=directory.toAbsolutePath().normalize();Files.createDirectories(this.directory);
        lockChannel=FileChannel.open(this.directory.resolve("owner.lock"),StandardOpenOption.CREATE,StandardOpenOption.WRITE);
        FileLock acquired=null;
        try {
            acquired=lockChannel.tryLock();
            if(acquired==null)throw new IOException("Workspace already has a live owner; attach to its HTTP service");
            Path file=this.directory.resolve("workspace.json");
            state=Files.isRegularFile(file)?JSON.readValue(file.toFile(),State.class):new State();
            if(state.schemaVersion!=1)throw new IOException("Unsupported workspace schema version");
            if(!Files.exists(file))persist(state);
        } catch(Exception e) {
            if(acquired!=null)acquired.close();lockChannel.close();
            if(e instanceof IOException io)throw io;
            throw new IOException("Cannot acquire/read workspace owner state",e);
        }
        lock=acquired;
    }
    public Path directory(){return directory;}
    public synchronized State snapshot() {
        ensureOpen();return copy(state);
    }
    public synchronized <T> T mutate(long expectedRevision,String actor,String operation,Function<State,T> operationBody) {
        ensureOpen();
        if(expectedRevision!=state.revision)throw new WorkbenchException("revision_conflict",409,"Workspace changed; expected "+expectedRevision+", current "+state.revision);
        WorkbenchException.require(actor!=null&&!actor.isBlank()&&actor.length()<=80,"An actor of at most 80 characters is required");
        State next=copy(state);T result=operationBody.apply(next);
        next.revision=state.revision+1;
        next.activity.add(new Activity(next.revision,actor,operation,Instant.now().toString()));
        if(next.activity.size()>1000)next.activity=new ArrayList<>(next.activity.subList(next.activity.size()-1000,next.activity.size()));
        try{persist(next);}catch(IOException e){throw new WorkbenchException("workspace_write_failed",500,"Workspace state was not changed: "+e.getMessage());}
        state=next;return result;
    }
    public synchronized <T> T internal(String operation,Function<State,T> body){return mutate(state.revision,"system",operation,body);}
    public static void pushUndo(State state) {
        state.undo.add(new UndoFrame(state.currentQueryId,state.selection));
        if(state.undo.size()>100)state.undo.remove(0);
    }
    private void persist(State value)throws IOException {
        Path temporary=directory.resolve("workspace.json.partial");JSON.writeValue(temporary.toFile(),value);
        Files.move(temporary,directory.resolve("workspace.json"),StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);
    }
    private static State copy(State s){return JSON.convertValue(s,State.class);}
    private void ensureOpen(){if(closed)throw new IllegalStateException("Workspace is closed");}
    @Override public synchronized void close()throws IOException{if(!closed){closed=true;try{lock.close();}finally{lockChannel.close();}}}
}
