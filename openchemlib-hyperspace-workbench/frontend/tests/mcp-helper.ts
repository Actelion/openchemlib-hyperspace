import {spawn} from 'node:child_process';
import {createInterface} from 'node:readline';
import {fileURLToPath} from 'node:url';
export async function connectMcp(connectionFile:string){
 const jar=fileURLToPath(new URL('../../../openchemlib-hyperspace-mcp/target/openchemlib-hyperspace-mcp.jar',import.meta.url));
 const child=spawn('java',['-jar',jar,'--workbench',connectionFile],{stdio:['pipe','pipe','pipe']});
 let sequence=0,stderr='';const pending=new Map<number,{resolve:(value:any)=>void;reject:(error:Error)=>void;timer:ReturnType<typeof setTimeout>}>();
 child.stderr.on('data',data=>{stderr=(stderr+data.toString()).slice(-8192);});
 createInterface({input:child.stdout}).on('line',line=>{let value;try{value=JSON.parse(line);}catch{return;}const item=pending.get(value.id);if(!item)return;pending.delete(value.id);clearTimeout(item.timer);value.error?item.reject(new Error(JSON.stringify(value.error))):item.resolve(value.result);});
 child.on('exit',()=>{for(const item of pending.values()){clearTimeout(item.timer);item.reject(new Error('MCP exited: '+stderr));}pending.clear();});
 const request=(method:string,params:Record<string,unknown>={})=>new Promise<any>((resolve,reject)=>{const id=++sequence;const timer=setTimeout(()=>{pending.delete(id);reject(new Error(`MCP timeout: ${method} ${stderr}`));},30000);pending.set(id,{resolve,reject,timer});child.stdin.write(JSON.stringify({jsonrpc:'2.0',id,method,params})+'\n');});
 await request('initialize',{protocolVersion:'2025-11-25',capabilities:{},clientInfo:{name:'workbench-browser-verifier',version:'1'}});
 child.stdin.write(JSON.stringify({jsonrpc:'2.0',method:'notifications/initialized'})+'\n');
 return {request,async call(name:string,args:Record<string,unknown>={}){const response=await request('tools/call',{name,arguments:args});if(response.isError)throw new Error(JSON.stringify(response));return response.structuredContent??JSON.parse(response.content[0].text);},close(){child.stdin.end();setTimeout(()=>{if(child.exitCode===null)child.kill('SIGTERM');},1000).unref();}};
}
