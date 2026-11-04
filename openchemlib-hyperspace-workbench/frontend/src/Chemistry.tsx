import {useEffect,useRef,useState} from 'react';
import {CanvasEditor,Molecule} from 'openchemlib';
import type {Document,Query,Selection} from './api';

export function Structure({document:doc,atoms=[],height=150}:{document:Document;atoms?:number[];height?:number}){
  try{const m=Molecule.fromIDCode(doc.idcode,doc.coordinates);m.ensureHelperArrays(Molecule.cHelperCIP);for(const a of atoms)if(a<m.getAllAtoms())m.setAtomColor(a,Molecule.cAtomColorMagenta);
    const svg=m.toSVG(320,height,undefined,{strokeWidth:1.5});return <img className="structure" height={height} src={`data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`} alt="Molecular structure with query matches in magenta"/>;
  }catch{return <div className="structure-error">Structure could not be depicted</div>;}
}
export function QueryCanvas({query,selection,onPublish,onSelection}:{query:Query|null;selection:Selection;onPublish:(doc:Document,baseId:string|null)=>Promise<unknown>;onSelection:(atoms:string[],bonds:string[])=>void}){
  const host=useRef<HTMLDivElement>(null),editor=useRef<CanvasEditor|null>(null),latest=useRef({query,selection,onSelection});latest.current={query,selection,onSelection};
  const base=useRef<Query|null>(null),dirtyRef=useRef(false),publishing=useRef(false);
  const [dirty,setDirty]=useState(false),[error,setError]=useState(''),[saving,setSaving]=useState(false),[conflict,setConflict]=useState(false);
  function load(current:Query|null){const e=editor.current;if(!e)return;const m=current?Molecule.fromIDCode(current.document.idcode,current.document.coordinates):new Molecule(0,0);if(!current)m.setFragment(true);e.setMolecule(m);base.current=current;dirtyRef.current=false;setDirty(false);setConflict(false);setError('');}
  useEffect(()=>{if(!host.current)return;const e=new CanvasEditor(host.current,{initialFragment:true});editor.current=e;
    e.setOnChangeListener(event=>{if(!event.isUserEvent)return;
      if(event.type==='molecule'){dirtyRef.current=true;setDirty(true);}
      if(event.type==='selection'&&!dirtyRef.current&&latest.current.query?.id===base.current?.id&&latest.current.query){const m=e.getMolecule(),id=latest.current.query.id;
        latest.current.onSelection(Array.from({length:m.getAllAtoms()},(_,i)=>i).filter(i=>m.isSelectedAtom(i)).map(i=>`${id}:a${i}`),Array.from({length:m.getAllBonds()},(_,i)=>i).filter(i=>m.isSelectedBond(i)).map(i=>`${id}:b${i}`));}
    });return ()=>{e.destroy();editor.current=null;};},[]);
  useEffect(()=>{if(dirtyRef.current&&!publishing.current&&query?.id!==base.current?.id){setConflict(true);return;}load(query);},[query?.id]);
  useEffect(()=>{if(dirtyRef.current||!editor.current||!query||query.id!==base.current?.id)return;const m=editor.current.getMolecule();for(let i=0;i<m.getAllAtoms();i++)m.setAtomSelection(i,selection.atoms.includes(`${query.id}:a${i}`));editor.current.moleculeChanged();},[selection.atoms.join('|'),query?.id]);
  async function publish(){const m=editor.current?.getMolecule();if(!m)return;setSaving(true);setError('');publishing.current=true;try{if(m.getAllAtoms()===0)throw new Error('Draw a structure first');await onPublish({idcode:m.getIDCode(),coordinates:m.getIDCoordinates(),fragment:m.isFragment()},base.current?.id??null);dirtyRef.current=false;setDirty(false);setConflict(false);load(latest.current.query);}catch(e){setError(String(e instanceof Error?e.message:e));}finally{publishing.current=false;setSaving(false);}}
  return <>{conflict&&<div className="notice draft-conflict" role="status"><b>Shared query changed; your drawing is preserved.</b><p>This unpublished draft is based on {base.current?.name??'an empty editor'}. Searches use the current published shared query. Publish this draft as a branch, or load the shared query.</p><button className="small" onClick={()=>load(query)}>Discard draft and load shared query</button></div>}<div className="canvas-wrap" ref={host} data-testid="query-canvas"/><div className="editor-foot"><span className={dirty?'draft':'muted'}>{dirty?'Unpublished drawing':query?.document.fragment?'OCL fragment query':'Concrete molecule'}</span><button className="small" onClick={publish} disabled={!dirty||saving}>{saving?'Publishing…':conflict?'Publish draft as branch':'Publish drawing'}</button></div>{error&&<p className="inline-error">{error}</p>}</>;
}
