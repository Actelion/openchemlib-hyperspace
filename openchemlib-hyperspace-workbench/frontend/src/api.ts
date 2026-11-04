export type Document={idcode:string;coordinates:string;fragment:boolean};
export type Query={id:string;parentId:string|null;name:string;document:Document;pose:Record<string,unknown>|null;actor:string;createdAt:string};
export type Atom={ref:string;index:number;element:string;charge:number;elements:number[];featuresHex:string;aromatic:boolean;notAromatic:boolean;ring:boolean;matchStereo:boolean;excludeElements:boolean;ringState:string};
export type Bond={ref:string;index:number;atomA:string;atomB:string;type:number;features:number;bridge:boolean;minAtoms:number|null;maxAtoms:number|null;allowedTypes:string[];aromaticState:string;ringState:string;matchStereo:boolean};
export type Inspection={revisionId:string;document:Document;atoms:Atom[];bonds:Bond[];warnings:string[]};
export type Selection={queryRevisionId:string|null;atoms:string[];bonds:string[];hits:string[]};
export type Asset={id:string;name:string;path:string;kind:string;libraryId:string;sourceSnapshotId:string;status:string;error:string|null;capabilities:string[];recordCount:number|null;coverage:string;metadata:Record<string,unknown>};
export type Run={id:string;mode:string;queryId:string;query:Query;assetIds:string[];settings:Record<string,unknown>;scope:string;state:string;retainedCount:number;uniqueMoleculeCount:number;diagnostics:Record<string,unknown>;createdAt:string};
export type Job={id:string;kind:string;runId:string|null;state:string;stage:string;actor:string;createdAt:string;finishedAt:string|null;request:Record<string,unknown>;output:Record<string,unknown>|null;error:string|null;log:string[]};
export type Collection={id:string;name:string;notes:string;hitCount:number;actor:string;createdAt:string;updatedAt:string};
export type Hit={id:string;runId:string;moleculeId:string;document:Document;smiles:string;assetId:string;familyId:string|null;reaction:string|null;synthons:{position:number;id:string;idcode:string}[];scores:Record<string,number>;properties:Record<string,number>;matchedAtoms:number[];provenance:Record<string,unknown>};
export type Family={id:string;runId:string;assetId:string;reaction:string;representedCombinations:string;positions:{position:number;count:number;sample:{position:number;id:string;idcode:string}[]}[]};
export type Page<T>={items:T[];offset:number;limit:number;total:number;scope:string};
export type State={runCount:number;collectionCount:number;jobCount:number;workspaceId:string;name:string;revision:number;sequence:number;currentQuery:Query|null;inspection:Inspection|null;selection:Selection;assets:Record<string,Asset>;libraries:Record<string,{id:string;name:string}>;sourceSnapshots:Record<string,{version:string;provenance:string}>;resourceStatus:Record<string,{loaded:boolean;activeJobs:number}>;runs:Run[];collections:Collection[];jobs:Job[];activity:{sequence:number;revision:number;actor:string;operation:string;at:string}[];queryHistory:{id:string;name:string;createdAt:string}[];view:'libraries'|'screening'|'saved';viewEntityId:string|null;undoAvailable:boolean;identityPolicy:string};
const fragment=new URLSearchParams(location.hash.slice(1));
if(fragment.has('token')){sessionStorage.setItem('hyperspace-token',fragment.get('token')!);history.replaceState(null,'',location.pathname);}
export const token=()=>sessionStorage.getItem('hyperspace-token')??'';
export class ApiError extends Error{constructor(public code:string,message:string,public revision:number){super(message);}}
export async function state():Promise<State>{const r=await fetch('/api/state',{headers:{'X-Workbench-Token':token()}});const body=await r.json();if(!r.ok)throw new ApiError(body.error,body.message,body.revision);return body.result;}
export async function call<T>(operation:string,args:Record<string,unknown>,expectedRevision:number,onRevision?:(revision:number)=>void):Promise<T>{
  const response=await fetch('/api/call',{method:'POST',headers:{'Content-Type':'application/json','X-Workbench-Token':token()},body:JSON.stringify({operation,arguments:args,expectedRevision,actor:'scientist:browser'})});
  const body=await response.json();if(!response.ok)throw new ApiError(body.error,body.message,body.revision);onRevision?.(body.revision);return body.result;
}
export async function download(exportId:string,file:string){const r=await fetch(`/api/export/${exportId}/${file}`,{headers:{'X-Workbench-Token':token()}});if(!r.ok)throw new Error('Export download failed');const url=URL.createObjectURL(await r.blob());const a=document.createElement('a');a.href=url;a.download=file;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
