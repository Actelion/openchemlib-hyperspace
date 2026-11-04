package com.idorsia.research.chem.hyperspace.workbench;

import com.fasterxml.jackson.databind.JsonNode;
import com.sun.net.httpserver.*;
import java.io.*;
import java.net.*;
import java.nio.file.*;
import java.nio.file.attribute.PosixFilePermissions;
import java.security.*;
import java.util.*;
import java.util.concurrent.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Loopback-only HTTP transport. A token is required for all API operations, including reads. */
public final class WorkbenchServer implements AutoCloseable {
    public final WorkbenchApplication application;
    private final HttpServer server;private final ExecutorService transport;
    private final String token;private final Path frontend,connectionFile;private final String authority;
    public WorkbenchServer(Path workspace,int port,Path frontend)throws Exception {
        this.application=new WorkbenchApplication(workspace);this.frontend=frontend.toAbsolutePath().normalize();
        try{server=HttpServer.create(new InetSocketAddress(InetAddress.getByName("127.0.0.1"),port),32);
        }catch(Exception e){application.close();throw e;}
        byte[] secret=new byte[32];new SecureRandom().nextBytes(secret);token=HexFormat.of().formatHex(secret);
        authority="127.0.0.1:"+server.getAddress().getPort();
        transport=Executors.newFixedThreadPool(8,r->{var t=new Thread(r,"workbench-http");t.setDaemon(true);return t;});server.setExecutor(transport);
        server.createContext("/",this::handle);server.start();
        connectionFile=application.store.directory().resolve("connection.json");
        try {
            if(!Files.exists(connectionFile))Files.createFile(connectionFile,PosixFilePermissions.asFileAttribute(PosixFilePermissions.fromString("rw-------")));
            else Files.setPosixFilePermissions(connectionFile,PosixFilePermissions.fromString("rw-------"));
            WorkspaceStore.JSON.writeValue(connectionFile.toFile(),Map.of("formatVersion",1,"workspaceId",application.store.snapshot().workspaceId,
                    "name",application.store.snapshot().name,"url",url(),"token",token,"pid",ProcessHandle.current().pid()));
        }catch(Exception e){server.stop(0);transport.shutdownNow();application.close();throw e;}
    }
    public String url(){return "http://"+authority;}
    public String token(){return token;}
    public String browserUrl(){return url()+"/#token="+token;}
    private void handle(HttpExchange exchange)throws IOException {
        try {
            String host=exchange.getRequestHeaders().getFirst("Host");require(authority.equals(host),"Unexpected Host header; use the configured loopback URL");
            String origin=exchange.getRequestHeaders().getFirst("Origin");require(origin==null||origin.equals(url()),"Foreign origins are not allowed");
            String path=exchange.getRequestURI().getPath();
            if(path.startsWith("/api/")) {
                String supplied=exchange.getRequestHeaders().getFirst("X-Workbench-Token");
                if(supplied==null||!MessageDigest.isEqual(token.getBytes(java.nio.charset.StandardCharsets.UTF_8),supplied.getBytes(java.nio.charset.StandardCharsets.UTF_8)))
                    throw new WorkbenchException("unauthorized",401,"Attach using this workspace's connection.json or launch URL");
                if(path.equals("/api/state")&&exchange.getRequestMethod().equals("GET")){json(exchange,200,Map.of("result",application.state(),"revision",application.store.snapshot().revision));return;}
                if(path.equals("/api/operations")&&exchange.getRequestMethod().equals("GET")){json(exchange,200,Map.of("operations",WorkbenchApplication.operations()));return;}
                if(path.startsWith("/api/export/")&&exchange.getRequestMethod().equals("GET")){download(exchange,path);return;}
                if(!path.equals("/api/call")||!exchange.getRequestMethod().equals("POST"))throw new WorkbenchException("endpoint_not_found",404,"Unknown API endpoint");
                byte[] data=exchange.getRequestBody().readNBytes(2_000_001);require(data.length<=2_000_000,"Request exceeds 2 MB");
                JsonNode request=WorkspaceStore.JSON.readTree(data);require(request!=null&&request.isObject(),"Expected a JSON object");
                String actor=request.path("actor").asText();require(!actor.isBlank()&&actor.length()<=80,"Provide an actor (max 80 characters)");
                String name=request.path("operation").asText();var operation=WorkbenchApplication.operations().stream().filter(o->o.name().equals(name)).findFirst()
                        .orElseThrow(()->new WorkbenchException("operation_not_found",404,"Unknown operation"));
                JsonNode args=request.path("arguments");require(args.isObject(),"Arguments must be an object");validate(operation,args);
                long revision=request.path("expectedRevision").asLong(-1);if(operation.writes())require(revision>=0,"Writes require expectedRevision from current workspace state");
                Object result=application.invoke(name,args,revision,actor);var response=new LinkedHashMap<String,Object>();response.put("result",result);response.put("revision",application.store.snapshot().revision);json(exchange,200,response);return;
            }
            require(exchange.getRequestMethod().equals("GET"),"Static assets require GET");
            Path file=frontend.resolve(path.equals("/")?"index.html":path.substring(1)).normalize();require(file.startsWith(frontend),"Invalid asset path");
            if(!Files.isRegularFile(file)){if(!path.contains("."))file=frontend.resolve("index.html");else throw new WorkbenchException("not_found",404,"Static asset not found");}
            if(!Files.isRegularFile(file))throw new WorkbenchException("frontend_missing",503,"Build frontend: cd openchemlib-hyperspace-workbench/frontend && npm ci && npm run build");
            String name=file.getFileName().toString();String type=name.endsWith(".html")?"text/html; charset=utf-8":name.endsWith(".js")?"text/javascript; charset=utf-8":name.endsWith(".css")?"text/css; charset=utf-8":name.endsWith(".svg")?"image/svg+xml":"application/octet-stream";
            exchange.getResponseHeaders().set("Content-Type",type);exchange.getResponseHeaders().set("Content-Security-Policy","default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'");
            exchange.getResponseHeaders().set("X-Content-Type-Options","nosniff");exchange.sendResponseHeaders(200,Files.size(file));try(var out=exchange.getResponseBody()){Files.copy(file,out);}
        }catch(WorkbenchException e){json(exchange,e.status(),Map.of("error",e.code(),"message",e.getMessage(),"revision",application.store.snapshot().revision));
        }catch(IllegalArgumentException e){json(exchange,400,Map.of("error","invalid_argument","message",Objects.toString(e.getMessage(),"Invalid argument"),"revision",application.store.snapshot().revision));
        }catch(Exception e){json(exchange,500,Map.of("error",e.getClass().getSimpleName(),"message",Objects.toString(e.getMessage(),"Operation failed"),"revision",application.store.snapshot().revision));
        }finally{exchange.close();}
    }
    @SuppressWarnings("unchecked") private static void validate(WorkbenchApplication.Operation operation,JsonNode args) {
        Map<String,Object> properties=(Map<String,Object>)operation.arguments().get("properties");
        args.fieldNames().forEachRemaining(key->require(properties.containsKey(key),"Unknown argument: "+key));
        for(String key:(List<String>)operation.arguments().get("required"))require(args.hasNonNull(key),"Missing argument: "+key);
        for(var entry:properties.entrySet())if(args.has(entry.getKey())) {
            JsonNode value=args.get(entry.getKey());Object type=((Map<String,Object>)entry.getValue()).get("type");
            List<?> types=type instanceof List<?> list?list:List.of(type);boolean valid=false;
            for(Object t:types)valid|=switch(t.toString()){case "null"->value.isNull();case "string"->value.isTextual();case "boolean"->value.isBoolean();case "integer"->value.isIntegralNumber();case "array"->value.isArray();case "object"->value.isObject();default->false;};
            require(valid,"Incorrect type for argument: "+entry.getKey());
        }
    }
    private void download(HttpExchange exchange,String path)throws IOException {
        String[] parts=path.split("/");require(parts.length==5&&parts[3].matches("export-[a-f0-9-]{36}")&&Set.of("results.sdf","results.tsv","summary.md","provenance.json").contains(parts[4]),"Unknown export file");
        Path file=application.store.directory().resolve("exports").resolve(parts[3]).resolve(parts[4]);if(!Files.isRegularFile(file))throw new WorkbenchException("not_found",404,"Export file is absent");
        exchange.getResponseHeaders().set("Content-Type","application/octet-stream");exchange.getResponseHeaders().set("Content-Disposition","attachment; filename=\""+parts[4]+"\"");
        exchange.sendResponseHeaders(200,Files.size(file));try(var out=exchange.getResponseBody()){Files.copy(file,out);}
    }
    private static void json(HttpExchange exchange,int status,Object result)throws IOException {
        byte[] bytes=WorkspaceStore.JSON.writeValueAsBytes(result);exchange.getResponseHeaders().set("Content-Type","application/json; charset=utf-8");exchange.getResponseHeaders().set("Cache-Control","no-store");exchange.getResponseHeaders().set("X-Content-Type-Options","nosniff");
        exchange.sendResponseHeaders(status,bytes.length);try(var out=exchange.getResponseBody()){out.write(bytes);}
    }
    @Override public void close()throws Exception{server.stop(1);transport.shutdownNow();try{application.close();}finally{Files.deleteIfExists(connectionFile);}}
    public static void main(String[] args)throws Exception {
        Path workspace=Path.of(".hyperspace-workbench"),frontend=Path.of("openchemlib-hyperspace-workbench/frontend/dist");int port=8765;
        for(int i=0;i<args.length;i+=2){require(i+1<args.length,"Options require values");switch(args[i]){case "--workspace"->workspace=Path.of(args[i+1]);case "--port"->port=Integer.parseInt(args[i+1]);case "--frontend"->frontend=Path.of(args[i+1]);default->throw new IllegalArgumentException("Unknown option: "+args[i]);}}
        var app=new WorkbenchServer(workspace,port,frontend);Runtime.getRuntime().addShutdownHook(new Thread(()->{try{app.close();}catch(Exception e){System.err.println(e.getMessage());}}));
        System.out.println("Hyperspace workspace: "+workspace.toAbsolutePath());System.out.println("Open browser: "+app.browserUrl());System.out.println("MCP connection: "+app.application.store.directory().resolve("connection.json"));
    }
}
