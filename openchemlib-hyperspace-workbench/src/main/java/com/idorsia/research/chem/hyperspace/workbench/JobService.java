package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.JsonNode;
import java.io.IOException;
import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Workspace-owned computations survive browser disconnects. Restart never silently replays a job. */
public final class JobService implements AutoCloseable {
    public record Job(String id,String kind,String runId,String state,String stage,String actor,
                      String createdAt,String finishedAt,JsonNode request,JsonNode output,String error,List<String> log){}
    @FunctionalInterface public interface Work {Object run(Context context)throws Exception;}
    private static final Set<String> ACTIVE=Set.of("queued","running","cancelling");
    private final WorkspaceStore store;
    private final ThreadPoolExecutor executor;
    private final Map<String,Control> controls=new ConcurrentHashMap<>();
    private final AtomicBoolean closing=new AtomicBoolean();
    private static final class Control {final AtomicBoolean cancelled=new AtomicBoolean();volatile Thread thread;}
    public final class Context {
        private final String id;private final Control control;
        private Context(String id,Control control){this.id=id;this.control=control;}
        public String id(){return id;}
        public boolean cancelled(){return control.cancelled.get()||Thread.currentThread().isInterrupted();}
        public void check(){if(cancelled())throw new CancellationException("Cancellation requested");}
        public void stage(String stage){check();update(id,j->copy(j,j.state(),stage,null,null,null));}
        public void detail(String key,Object value){
            update(id,j->{var request=j.request().deepCopy();((com.fasterxml.jackson.databind.node.ObjectNode)request).set(key,WorkspaceStore.JSON.valueToTree(value));
                return new Job(j.id(),j.kind(),j.runId(),j.state(),j.stage(),j.actor(),j.createdAt(),j.finishedAt(),request,j.output(),j.error(),j.log());});
        }
        public void log(String message){
            require(message!=null&&message.length()<=4000,"Log line exceeds bound");
            update(id,j->{var lines=new ArrayList<>(j.log());lines.add(Instant.now()+" "+message);if(lines.size()>100)lines.remove(0);
                return new Job(j.id(),j.kind(),j.runId(),j.state(),j.stage(),j.actor(),j.createdAt(),j.finishedAt(),j.request(),j.output(),j.error(),lines);});
        }
    }
    public JobService(WorkspaceStore store,int workers,int queueCapacity)throws IOException {
        require(workers>=1&&workers<=8&&queueCapacity>=1&&queueCapacity<=64,"Invalid application job limits");
        this.store=store;Files.createDirectories(store.directory().resolve("jobs"));
        executor=new ThreadPoolExecutor(workers,workers,0,TimeUnit.SECONDS,new ArrayBlockingQueue<>(queueCapacity),
                r->{var t=new Thread(r,"workbench-job");t.setDaemon(false);return t;},new ThreadPoolExecutor.AbortPolicy());
        if(store.snapshot().jobs.values().stream().anyMatch(n->ACTIVE.contains(n.path("state").asText())))
            store.internal("jobs.owner_restarted",s->{s.jobs.replaceAll((id,node)->{
                Job j=WorkspaceStore.JSON.convertValue(node,Job.class);
                return ACTIVE.contains(j.state())?WorkspaceStore.JSON.valueToTree(copy(j,"interrupted","owner restarted",Instant.now().toString(),null,
                        "Application owner stopped. No computation was restarted; submit a new job or resume its documented preparation artifacts.")):node;
            });return null;});
    }
    /** Leases are acquired by the application before submit, held through queued/running work, and always released. */
    public Job submit(long expected,String actor,String kind,String runId,Object request,List<? extends AutoCloseable> leases,Work work) {
        require(!closing.get(),"Application is shutting down");
        require(kind!=null&&!kind.isBlank()&&work!=null,"Job kind and work required");
        String id="job-"+UUID.randomUUID();Control control=new Control();
        Job job=new Job(id,kind,runId,"queued","queued",actor,Instant.now().toString(),null,WorkspaceStore.JSON.valueToTree(request),null,null,List.of());
        try {
            store.mutate(expected,actor,"job.submit",s->{s.jobs.put(id,WorkspaceStore.JSON.valueToTree(job));return null;});
            controls.put(id,control);
            executor.execute(()->execute(job,control,leases,work));
        } catch(RuntimeException e) {
            controls.remove(id);closeLeases(leases);
            if(e instanceof RejectedExecutionException){update(id,j->copy(j,"failed","queue full",Instant.now().toString(),null,"Job queue is full; retry when a job finishes"));
                throw new WorkbenchException("queue_full",429,"Job queue is full");}
            throw e;
        }
        return job;
    }
    private void execute(Job initial,Control control,List<? extends AutoCloseable> leases,Work work) {
        control.thread=Thread.currentThread();String terminal="failed",error=null;JsonNode output=null;
        try {
            new Context(initial.id(),control).check();update(initial.id(),j->copy(j,"running","starting",null,null,null));
            Object result=work.run(new Context(initial.id(),control));
            new Context(initial.id(),control).check();output=WorkspaceStore.JSON.valueToTree(result);terminal="completed";
        } catch(CancellationException|InterruptedException e){terminal="cancelled";error="Cancellation requested; completed artifacts remain available where documented.";
        } catch(Exception e){if(control.cancelled.get())terminal="cancelled";error=e.getClass().getSimpleName()+": "+Objects.toString(e.getMessage(),"No detail");
        } finally {
            try{closeLeases(leases);}catch(RuntimeException e){terminal="failed";error="Resource cleanup failed: "+e.getMessage();}
            final String proposedState=terminal,detail=error;final JsonNode result=output;
            update(initial.id(),j->{String state=control.cancelled.get()&&!proposedState.equals("failed")?"cancelled":proposedState;
                return copy(j,state,state,Instant.now().toString(),result,state.equals("cancelled")&&detail==null?"Cancellation requested; completed artifacts remain available.":detail);});
            final Job finished=get(initial.id());final String state=finished.state();
            if(initial.runId()!=null&&!state.equals("completed"))store.internal("run."+state,s->{
                var node=s.runs.get(initial.runId());if(node instanceof com.fasterxml.jackson.databind.node.ObjectNode run){
                    run.put("state",state);run.set("diagnostics",WorkspaceStore.JSON.valueToTree(Map.of("jobId",initial.id(),"error",Objects.toString(finished.error(),""))));
                }return null;
            });
            controls.remove(initial.id());
            // An interrupted worker must not contaminate the next computation on this pooled thread.
            Thread.interrupted();
        }
    }
    public Job get(String id){JsonNode n=store.snapshot().jobs.get(id);if(n==null)throw new WorkbenchException("job_not_found",404,"Unknown job: "+id);return WorkspaceStore.JSON.convertValue(n,Job.class);}
    public List<Job> list(){return store.snapshot().jobs.values().stream().map(n->WorkspaceStore.JSON.convertValue(n,Job.class)).toList();}
    public Job cancel(long expected,String actor,String id) {
        get(id);Control control=controls.get(id);require(control!=null,"No live job handle; refresh status");
        store.mutate(expected,actor,"job.cancel",s->{
            Job current=WorkspaceStore.JSON.convertValue(s.jobs.get(id),Job.class);require(ACTIVE.contains(current.state()),"Job is already terminal");
            // Set the flag while holding the same lock as terminal publication. A late cancel cannot resurrect a completed job or drop a worker handle.
            control.cancelled.set(true);s.jobs.put(id,WorkspaceStore.JSON.valueToTree(copy(current,"cancelling","cancellation requested",null,null,null)));return null;
        });
        // Cooperative cancellation keeps the lease until the underlying engine's workers have actually returned.
        return get(id);
    }
    private void update(String id,java.util.function.UnaryOperator<Job> operation) {
        store.internal("job.progress",s->{JsonNode node=s.jobs.get(id);require(node!=null,"Job metadata missing");
            s.jobs.put(id,WorkspaceStore.JSON.valueToTree(operation.apply(WorkspaceStore.JSON.convertValue(node,Job.class))));return null;});
    }
    private static Job copy(Job j,String state,String stage,String finished,JsonNode output,String error){return new Job(j.id(),j.kind(),j.runId(),state,stage,j.actor(),j.createdAt(),finished,j.request(),output,error,j.log());}
    private static void closeLeases(List<? extends AutoCloseable> leases){RuntimeException failed=null;
        for(int i=leases.size()-1;i>=0;i--)try{leases.get(i).close();}catch(Exception e){if(failed==null)failed=new IllegalStateException(e);else failed.addSuppressed(e);}
        if(failed!=null)throw failed;
    }
    @Override public void close() {
        closing.set(true);controls.values().forEach(c->c.cancelled.set(true));executor.shutdown();
        boolean interrupted=false;
        while(!executor.isTerminated())try{executor.awaitTermination(1,TimeUnit.SECONDS);}catch(InterruptedException e){interrupted=true;}
        if(interrupted)Thread.currentThread().interrupt();
    }
}
