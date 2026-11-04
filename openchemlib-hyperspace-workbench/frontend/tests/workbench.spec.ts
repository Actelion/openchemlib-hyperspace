import {test,expect} from '@playwright/test';
import {readFileSync} from 'node:fs';
const connection=JSON.parse(readFileSync(process.env.WORKBENCH_CONNECTION??'/tmp/hyperspace-workbench-browser/connection.json','utf8'));
test('scientist toy search, expansion, collection, export and reconnect',async({page})=>{
 const errors:string[]=[];page.on('pageerror',error=>errors.push(error.message));
 await page.goto(connection.url+'/#token='+connection.token);
 await expect(page.getByRole('heading',{name:'Query studio'})).toBeVisible();
 await page.getByRole('button',{name:'Libraries',exact:true}).click();
 await page.getByRole('button',{name:'Prepare toy library',exact:true}).click();
 await expect(page.getByRole('heading',{name:'Toy · substructure',exact:true})).toBeVisible({timeout:30000});
 const card=page.locator('.asset-card').filter({has:page.getByRole('heading',{name:'Toy · substructure',exact:true})});
 await card.getByRole('button',{name:'Load resource',exact:true}).click();
 await expect(card.getByText('ready',{exact:true})).toBeVisible();
 await card.getByRole('button',{name:'Use in screening'}).click();
 // Hide the job drawer to verify the full layout and canvas.
 await page.locator('.drawer-bar>button').click();
 await page.getByText('Paste or import a structure',{exact:true}).click();
 await page.getByRole('textbox',{name:'Query structure',exact:true}).fill('CCO');
 await page.getByRole('button',{name:'Publish query',exact:true}).click();
 await expect(page.getByText('OCL fragment query',{exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Search combinatorial space'}).click();
 await expect(page.locator('.family-card').first()).toBeVisible({timeout:30000});
 const familyBoxes=page.locator('.family-card input[type=checkbox]');
 for(let i=0;i<await familyBoxes.count();i++)await familyBoxes.nth(i).check();
 await page.getByRole('spinbutton',{name:'Expansion limit'}).fill('4');
 await page.getByRole('button',{name:'Expand selected'}).click();
 await expect(page.locator('.hit-card')).toHaveCount(4,{timeout:30000});
 await page.locator('.drawer-bar>button').click();
 await page.getByRole('button',{name:'Select page',exact:true}).click();
 await expect(page.getByRole('button',{name:'Save 4 to collection',exact:true})).toBeEnabled();
 await page.locator('.hit-card').first().click();
 await expect(page.getByRole('heading',{name:'Selected occurrence',exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Select page',exact:true}).click();
 await page.locator('.query-panel').evaluate(element=>{element.scrollTop=0;});
 await page.screenshot({path:'test-results/workbench-structures-light.png',fullPage:true});
 await page.getByRole('button',{name:'Toggle appearance'}).click();
 await expect(page.getByRole('button',{name:'Export ↗',exact:true})).toHaveCSS('background-color','rgb(28, 32, 42)');
 await page.screenshot({path:'test-results/workbench-structures-dark.png',fullPage:true});
 await page.getByRole('button',{name:'Toggle appearance'}).click();
 await page.getByRole('button',{name:'☷ Table',exact:true}).click();
 await expect(page.locator('tbody tr')).toHaveCount(4);
 await page.getByRole('button',{name:'⋮ Property map',exact:true}).click();
 await expect(page.getByRole('img',{name:'Linked molecular weight versus rotatable bonds plot'})).toBeVisible();
 await page.getByRole('button',{name:'Save 4 to collection',exact:true}).click();
 const dialog=page.getByRole('dialog',{name:'Save a collection'});
 await dialog.getByLabel('Name',{exact:true}).fill('Toy alcohol collection');
 await dialog.getByLabel('Scientific notes',{exact:true}).fill('Bounded example; source tuples preserved.');
 await dialog.getByRole('button',{name:'Save collection',exact:true}).click();
 await expect(page.getByRole('heading',{name:'Saved runs & collections',exact:true})).toBeVisible();
 await page.locator('.collection-card').filter({hasText:'Toy alcohol collection'}).getByRole('button',{name:'Open'}).click();
 await expect(page.getByRole('heading',{name:'Toy alcohol collection',exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Export ↗',exact:true}).click();
 await expect(page.getByRole('dialog',{name:'Export ready'})).toBeVisible();
 const saved=page.waitForEvent('download');await page.getByRole('button',{name:'↓ results.tsv',exact:true}).click();
 expect((await saved).suggestedFilename()).toBe('results.tsv');
 await page.getByRole('button',{name:'Close dialog'}).click();await page.reload();
 await expect(page.getByRole('heading',{name:'Toy alcohol collection',exact:true})).toBeVisible();
 await expect(page.getByText('Bounded example; source tuples preserved.',{exact:true})).toBeVisible();
 expect(errors).toEqual([]);
});

test('MCP agent edits scientist selection, searches and updates the visible collection',async({page})=>{
 const {connectMcp}=await import('./mcp-helper');
 const mcp=await connectMcp(process.env.WORKBENCH_CONNECTION??'/tmp/hyperspace-workbench-browser/connection.json');
 try {
  const tools=await mcp.request('tools/list');expect(tools.tools.some((t:{name:string})=>t.name==='wb_query_edit_selection')).toBe(true);
  let state=(await mcp.call('wb_workspace_state')).result;
  await mcp.call('wb_workspace_navigate',{view:'screening',expectedRevision:state.revision});
  state=(await mcp.call('wb_workspace_state')).result;
  const published=await mcp.call('wb_query_create',{format:'smiles',structure:'CCO',fragment:true,name:'Agent alcohol query',expectedRevision:state.revision});
  const queryId=published.result.query.id;
  await page.goto(connection.url+'/#token='+connection.token);
  await expect(page.locator('.query-name')).toContainText('Agent alcohol query');
  await page.getByText('Query features & shared selection',{exact:false}).click();
  const inspection=(await mcp.call('wb_query_inspect',{queryId})).result;
  const oxygen=inspection.atoms.find((a:{element:string})=>a.element==='O');
  await page.locator('.atom-ref-list button').filter({hasText:new RegExp(`^${oxygen.index+1} O`)}).click();
  await expect(page.locator('.atom-ref-list button.selected')).toHaveCount(1);
  state=(await mcp.call('wb_workspace_state')).result;
  expect(state.selection.atoms).toContain(oxygen.ref);
  const changed=await mcp.call('wb_query_edit_selection',{operation:'ATOM_ELEMENTS',elements:[7,8],expectedRevision:state.revision});
  expect(changed.result.oldReferencesInvalidated).toBe(true);
  await expect(page.locator('.query-name')).toContainText('Edited selection');
  await expect(page.locator('.atom-ref-list button').filter({hasText:'[7,8]'})).toBeVisible();
  state=(await mcp.call('wb_workspace_state')).result;
  const index=Object.values(state.assets).find((a:any)=>a.kind==='SUBSTRUCTURE_INDEX'&&a.status==='ready') as {id:string};expect(index).toBeTruthy();
  const started=await mcp.call('wb_search_substructure',{queryId:changed.result.query.id,assetIds:[index.id],threads:2,familyLimit:100,expectedRevision:state.revision});
  const waitJob=async(id:string)=>{const deadline=Date.now()+30000;while(Date.now()<deadline){const job=(await mcp.call('wb_job_status',{jobId:id})).result;if(['completed','failed','cancelled'].includes(job.state)){expect(job.state,job.error??'').toBe('completed');return job;}await new Promise(r=>setTimeout(r,25));}throw new Error('Owner job timeout');};
  await waitJob(started.result.id);
  const families=(await mcp.call('wb_run_families',{runId:started.result.runId,limit:100})).result;
  expect(families.total).toBeGreaterThan(0);
  state=(await mcp.call('wb_workspace_state')).result;
  const expanded=await mcp.call('wb_search_expand',{runId:started.result.runId,familyIds:families.items.map((f:{id:string})=>f.id),maxProducts:4,expectedRevision:state.revision});
  await waitJob(expanded.result.id);
  const hits=(await mcp.call('wb_run_results',{runId:expanded.result.runId,limit:10})).result;
  expect(hits.items.length).toBeGreaterThan(0);
  state=(await mcp.call('wb_workspace_state')).result;
  await mcp.call('wb_selection_change',{queryRevisionId:state.currentQuery.id,atoms:[],bonds:[],hits:[hits.items[0].id],expectedRevision:state.revision});
  await expect(page.getByRole('heading',{name:'Selected occurrence',exact:true})).toBeVisible();
  state=(await mcp.call('wb_workspace_state')).result;
  const collection=await mcp.call('wb_collection_save',{name:'Agent curated collection',notes:'Shared MCP/browser verification',hitIds:hits.items.map((h:{id:string})=>h.id),expectedRevision:state.revision});
  await expect(page.getByRole('heading',{name:'Agent curated collection',exact:true})).toBeVisible();
  state=(await mcp.call('wb_workspace_state')).result;
  await mcp.call('wb_workspace_navigate',{view:'screening',entityId:collection.result.id,expectedRevision:state.revision});
  await expect(page.getByRole('heading',{name:'Agent curated collection',exact:true})).toBeVisible();
  await page.screenshot({path:'test-results/workbench-agent-collection.png',fullPage:true});
 }finally{mcp.close();}
});

test('real learned 2D, ensemble 3D, cascade, bounded ranking and supplied pose in the browser',async({page})=>{
 const {execFileSync}=await import('node:child_process');const {mkdtempSync,writeFileSync,existsSync}=await import('node:fs');const {tmpdir}=await import('node:os');const {resolve}=await import('node:path');
 const root=resolve(process.cwd(),'..','..'),example=mkdtempSync(tmpdir()+'/hyperspace-browser-learned-'),cache=example+'/cache';
 const model=root+'/openchemlib-hyperspace-3d/model-bundles/deepspace7-v1',compact=root+'/openchemlib-hyperspace-3d/model-bundles/deepspace7-skelspheres16';
 writeFileSync(example+'/molecules.tsv','id\tsmiles\nhexane\tCCCCCC\nhexanol\tCCCCCCO\nhexanol-copy\tOCCCCCC\nbenzene\tc1ccccc1\ncaffeine\tCn1c(=O)c2c(ncn2C)n(C)c1=O\n');
 writeFileSync(example+'/build.json',JSON.stringify({formatVersion:1,inputs:{library:example+'/molecules.tsv',modelBundle:model,compactBundle:compact},output:{directory:cache,sourceRowsPerShard:100},runtime:{device:'CPU',cpuWorkers:2,encoderBatchSize:8,queueCapacity:2}}));
 execFileSync('java',['-cp',root+'/openchemlib-hyperspace-3d/target/openchemlib-hyperspace-3d-3.1.0.jar:'+root+'/openchemlib-hyperspace-3d/target/lib/*','com.idorsia.research.chem.hyperspace3d.cli.MoleculeFingerprintIndexBuilderCLI','--config',example+'/build.json'],{timeout:30000,stdio:'pipe'});
 await page.goto(connection.url+'/#token='+connection.token);
 await page.getByRole('button',{name:'Libraries',exact:true}).click();
 async function register(path:string,name:string){await page.getByRole('button',{name:'+ Add data',exact:true}).click();const dialog=page.getByRole('dialog',{name:'Add local data'});await dialog.getByRole('textbox',{name:'Asset path',exact:true}).fill(path);await dialog.getByRole('button',{name:'Inspect prepared asset',exact:true}).click();await dialog.getByRole('textbox',{name:'Asset name',exact:true}).fill(name);await dialog.getByRole('button',{name:'Register without loading',exact:true}).click();const card=page.locator('.asset-card').filter({has:page.getByRole('heading',{name,exact:true})});await expect(card).toBeVisible();await card.getByRole('button',{name:'Load resource',exact:true}).click();await expect(card.getByText('ready',{exact:true})).toBeVisible();}
 await register(model,'Canonical test graph');await register(compact,'Compact test projection');await register(cache,'Public learned test cache');
 const externalPose=process.env.WORKBENCH_POSE_BUNDLE??resolve(root,'..','pose_deployment/shared-candidate-seed17-verified');
 if(existsSync(externalPose+'/manifest.json'))await register(externalPose,'External verified pose');
 await page.getByRole('button',{name:'Screening',exact:true}).click();
 await page.getByText('Paste or import a structure',{exact:true}).click();await page.getByRole('textbox',{name:'Query structure',exact:true}).fill('CCCCCCO');await page.getByRole('checkbox',{name:'Fragment query',exact:true}).uncheck();await page.getByRole('button',{name:'Publish query',exact:true}).click();
 await page.getByText('Learned 2D / 3D screening',{exact:true}).click();await page.getByRole('spinbutton',{name:'Learned retained top K'}).fill('3');
 await page.getByRole('button',{name:'Validate model / cache',exact:true}).click();await expect(page.getByRole('status')).toContainText('Compatible · 5 accepted cache rows');
 async function screen(mode:string,title:string){await page.getByRole('combobox',{name:'Learned search mode',exact:true}).selectOption(mode);if(mode.includes('SKELSPHERES'))await page.getByRole('spinbutton',{name:'Learned shortlist size'}).fill('5');await page.getByRole('button',{name:'Screen fingerprint cache →',exact:true}).click();await expect(page.locator('.results-header .eyebrow')).toHaveText(title);await expect(page.locator('.results-header .tag')).toHaveText('completed');await expect(page.locator('.hit-card')).toHaveCount(3);}
 await screen('PHESA_3D','PHESA 3D');await expect(page.locator('.scope-bar>p')).toContainText('Exhaustive');await page.locator('.hit-card').first().click();await expect(page.getByRole('heading',{name:'Raw score components',exact:true})).toBeVisible();
 await screen('SKELSPHERES_2D','SKELSPHERES 2D');await expect(page.locator('.scope-bar>p')).toContainText('exact OCL SkelSpheres');
 await screen('SKELSPHERES_THEN_PHESA','SKELSPHERES THEN PHESA');await expect(page.locator('.scope-bar>p')).toContainText('16D top-5 shortlist');
 await page.getByRole('button',{name:'Saved & collections',exact:true}).click();await page.locator('.saved-run').filter({hasText:'Imported query'}).filter({hasText:'bounded expansion'}).click();await expect(page.locator('.results-header .eyebrow')).toHaveText('BOUNDED EXPANSION');await page.getByText('Learned 2D / 3D screening',{exact:true}).click();await page.getByRole('spinbutton',{name:'Learned retained top K'}).fill('3');
 await page.getByRole('combobox',{name:'Learned search mode',exact:true}).selectOption('PHESA_3D');await page.getByRole('spinbutton',{name:'Bounded ranking candidate limit'}).fill('3');await page.getByRole('button',{name:'Rank bounded candidates →',exact:true}).click();await expect(page.locator('.results-header .eyebrow')).toHaveText('BOUNDED LEARNED RANKING');await expect(page.locator('.results-header .tag')).toHaveText('completed');await expect(page.locator('.hit-card')).toHaveCount(3);await expect(page.locator('.scope-bar>p')).toContainText('Not exhaustive');await page.locator('.hit-card').first().click();await expect(page.locator('.inspector-synthon')).toHaveCount(2);
 if(existsSync(externalPose+'/manifest.json')){await page.getByRole('button',{name:'Supplied pose',exact:true}).click();await page.getByRole('textbox',{name:'Pose file path'}).fill(root+'/openchemlib-hyperspace-3d/src/test/resources/pose-golden/query-0.mol');await page.getByRole('button',{name:'Validate & attach pose',exact:true}).click();await expect(page.getByText('Pose attached',{exact:true})).toBeVisible();await page.getByRole('combobox',{name:'Learned search mode',exact:true}).selectOption('POSE_PHESA_3D');await page.getByRole('spinbutton',{name:'Pharmacophore weight',exact:true}).fill('0.25');await screen('POSE_PHESA_3D','POSE PHESA 3D');await expect(page.locator('.scope-bar>p')).toContainText('directed');await expect(page.locator('.scope-bar>p')).toContainText('No candidate conformers or computed 3D alignment');}
 await page.locator('.query-panel').evaluate(element=>{element.scrollTop=element.scrollHeight;});await page.screenshot({path:'test-results/workbench-learned-screening.png',fullPage:true});
 await expect(page.getByRole('alert')).toHaveCount(0);
});

test('guided supplier and fingerprint preparation use real detached workers and auto-register outputs',async({page})=>{
 const {mkdtempSync,writeFileSync,existsSync}=await import('node:fs');const {tmpdir}=await import('node:os');const {resolve}=await import('node:path');const root=resolve(process.cwd(),'..','..'),source=mkdtempSync(tmpdir()+'/hyperspace-browser-source-');
 writeFileSync(source+'/synthons.tsv','SMILES\tfragmentId\tsynthonSetIndex\treactionId\n[U]CCC\ta\t0\tpublic\n[U]CCN\tb\t0\tpublic\n[U]CCO\tx\t1\tpublic\n[U]CCC\ty\t1\tpublic\n');
 writeFileSync(source+'/browser-molecules.tsv','id\tsmiles\nhexanol\tCCCCCCO\nbenzene\tc1ccccc1\nsmall\tCCO\n');
 await page.goto(connection.url+'/#token='+connection.token);await page.getByRole('button',{name:'Libraries',exact:true}).click();
 async function prepareDialog(){await page.getByRole('button',{name:'+ Add data',exact:true}).click();const dialog=page.getByRole('dialog',{name:'Add local data'});await dialog.getByRole('button',{name:'Prepare source data',exact:true}).click();await expect(dialog.getByRole('spinbutton',{name:'Preparation CPU workers'})).toBeVisible();return dialog;}
 let dialog=await prepareDialog();await dialog.getByRole('textbox',{name:'Preparation source path'}).fill(source+'/synthons.tsv');await dialog.getByRole('textbox',{name:'Preparation library name'}).fill('Prepared browser source');await dialog.getByRole('textbox',{name:'Preparation child heap'}).fill('512M');
 await dialog.getByRole('button',{name:'Inspect bounded source sample'}).click();await expect(dialog.getByText('Sample and format diagnostics',{exact:true})).toBeVisible();
 await dialog.getByRole('button',{name:'Preview preparation recipe'}).click();await expect(dialog.getByText('Restart / resume',{exact:true})).toBeVisible();await expect(dialog.locator('.recipe-preview')).toContainText('2 CPU workers · 512M child heap');
 await dialog.getByRole('button',{name:'Start detached preparation →',exact:true}).click();
 const imported=page.locator('.asset-card').filter({has:page.getByRole('heading',{name:'Prepared browser source · substructure',exact:true})});await expect(imported).toBeVisible({timeout:30000});await expect(imported.getByText('registered',{exact:true})).toBeVisible();
 await imported.getByRole('button',{name:'Load resource',exact:true}).click();await expect(imported.getByText('ready',{exact:true})).toBeVisible();
 dialog=await prepareDialog();await dialog.getByRole('combobox',{name:'Preparation source type'}).selectOption('fingerprints');await dialog.getByRole('textbox',{name:'Preparation source path'}).fill(source+'/browser-molecules.tsv');await dialog.getByRole('textbox',{name:'Fingerprint output directory'}).fill(source+'/fingerprints');await dialog.getByRole('textbox',{name:'Preparation child heap'}).fill('512M');await dialog.getByRole('spinbutton',{name:'Preparation encoder batch'}).fill('8');
 await dialog.getByRole('button',{name:'Inspect bounded source sample'}).click();await dialog.getByRole('button',{name:'Preview preparation recipe'}).click();await expect(dialog.locator('.recipe-preview')).toContainText('committed source-order fingerprint shards');await page.screenshot({path:'test-results/workbench-preparation-recipe.png',fullPage:true});await dialog.getByRole('button',{name:'Start detached preparation →',exact:true}).click();
 const cache=page.locator('.asset-card').filter({has:page.getByRole('heading',{name:'browser-molecules.tsv · moleculeIndex',exact:true})});await expect(cache).toBeVisible({timeout:30000});await expect(cache.getByText('registered',{exact:true})).toBeVisible();await expect(cache.locator('.metric-row')).toContainText('2');expect(existsSync(source+'/fingerprints/manifest.json')).toBe(true);
 dialog=await prepareDialog();await dialog.getByText('Export a Slurm bundle · no submission',{exact:true}).click();await dialog.getByRole('textbox',{name:'Slurm profile path'}).fill(root+'/fingerprint-slurm-profile.example.json');await dialog.getByRole('textbox',{name:'Slurm bundle output directory'}).fill(source+'/slurm');await dialog.getByRole('button',{name:'Export bundle only',exact:true}).click();await expect(dialog.locator('.slurm-export pre')).toContainText('"submitted": false');expect(existsSync(source+'/slurm/compute.sbatch')).toBe(true);await dialog.getByRole('button',{name:'Close dialog'}).click();
 await page.screenshot({path:'test-results/workbench-prepared-library.png',fullPage:true});await expect(page.getByRole('alert')).toHaveCount(0);
});

test('named workspace, saved-run branching and collection edits preserve shared provenance',async({page})=>{
 const {connectMcp}=await import('./mcp-helper');const mcp=await connectMcp(process.env.WORKBENCH_CONNECTION??'/tmp/hyperspace-workbench-browser/connection.json');
 try{
  await page.goto(connection.url+'/#token='+connection.token);await page.getByRole('button',{name:'Rename workspace',exact:true}).click();const naming=page.getByRole('dialog',{name:'Name this workspace'});await naming.getByLabel('Workspace name',{exact:true}).fill('Shared discovery workspace');await naming.getByRole('button',{name:'Save workspace name'}).click();await expect(page.getByRole('button',{name:'Rename workspace',exact:true})).toContainText('Shared discovery workspace');
  let state=(await mcp.call('wb_workspace_state')).result;expect(state.name).toBe('Shared discovery workspace');const original=state.runs.find((r:any)=>r.mode==='bounded_expansion');expect(original).toBeTruthy();
  await mcp.call('wb_workspace_navigate',{view:'screening',entityId:original.id,expectedRevision:state.revision});await expect(page.getByRole('button',{name:'Branch run query',exact:true})).toBeVisible();await page.getByRole('button',{name:'Branch run query',exact:true}).click();await expect(page.locator('.query-name')).toContainText('Branch ·');
  state=(await mcp.call('wb_workspace_state')).result;expect(state.currentQuery.parentId).toBe(original.queryId);expect(state.runs.find((r:any)=>r.id===original.id).queryId).toBe(original.queryId);
  await page.getByText('Revision history & branches',{exact:true}).click();await page.getByRole('button',{name:'Load query history',exact:true}).click();await expect(page.locator('.history-list>div').first()).toContainText('Branch ·');
  const collection=state.collections.find((c:any)=>c.name==='Toy alcohol collection');expect(collection.hitCount).toBe(4);await page.getByRole('button',{name:'Saved & collections',exact:true}).click();await page.locator('.collection-card').filter({hasText:'Toy alcohol collection'}).getByRole('button',{name:'Open'}).click();await page.getByRole('button',{name:'Edit collection',exact:true}).click();
  let edit=page.getByRole('dialog',{name:'Edit collection',exact:true});await edit.getByLabel('Name',{exact:true}).fill('Reviewed alcohols');await edit.getByLabel('Scientific notes',{exact:true}).fill('Reviewed together; original members retained.');await edit.getByRole('button',{name:'Save collection',exact:true}).click();
  state=(await mcp.call('wb_workspace_state')).result;const updated=state.collections.find((c:any)=>c.id===collection.id);expect(updated.hitCount).toBe(4);expect(updated.createdAt).toBe(collection.createdAt);expect(updated.notes).toContain('Reviewed together');
  await page.locator('.collection-card').filter({hasText:'Reviewed alcohols'}).getByRole('button',{name:'Open'}).click();await expect(page.locator('.hit-card')).toHaveCount(4);await page.locator('.hit-card').first().click();await expect(page.getByRole('button',{name:'Save 1 to collection',exact:true})).toBeEnabled();await page.getByRole('button',{name:'Edit collection',exact:true}).click();edit=page.getByRole('dialog',{name:'Edit collection',exact:true});await edit.getByRole('checkbox',{name:/Replace members with current shared selection/}).check();await edit.getByRole('button',{name:'Save collection',exact:true}).click();
  state=(await mcp.call('wb_workspace_state')).result;expect(state.collections.find((c:any)=>c.id===collection.id).hitCount).toBe(1);await page.locator('.collection-card').filter({hasText:'Reviewed alcohols'}).getByRole('button',{name:'Open'}).click();await expect(page.locator('.hit-card')).toHaveCount(1);await page.reload();await expect(page.locator('.collection-note')).toContainText('Reviewed together');await expect(page.locator('.hit-card')).toHaveCount(1);
 }finally{mcp.close();}
});

test('an unpublished OCL drawing survives an agent edit and publishes as an explicit branch',async({page})=>{
 const {connectMcp}=await import('./mcp-helper');const mcp=await connectMcp(process.env.WORKBENCH_CONNECTION??'/tmp/hyperspace-workbench-browser/connection.json');
 try{
  let state=(await mcp.call('wb_workspace_state')).result;await mcp.call('wb_workspace_navigate',{view:'screening',expectedRevision:state.revision});state=(await mcp.call('wb_workspace_state')).result;
  const base=(await mcp.call('wb_query_create',{format:'smiles',structure:'CCO',fragment:true,name:'Draft base',expectedRevision:state.revision})).result.query;
  await page.goto(connection.url+'/#token='+connection.token);await expect(page.locator('.query-name')).toContainText('Draft base');const canvas=page.locator('[data-testid=query-canvas] canvas').first();const bounds=(await canvas.boundingBox())!;
  await page.mouse.click(bounds.x+26,bounds.y+121);await page.mouse.move(bounds.x+130,bounds.y+58);await page.mouse.down();await page.mouse.move(bounds.x+190,bounds.y+90,{steps:8});await page.mouse.up();await expect(page.getByRole('button',{name:'Publish drawing',exact:true})).toBeEnabled();
  state=(await mcp.call('wb_workspace_state')).result;await mcp.call('wb_query_create',{format:'smiles',structure:'c1ccccc1',fragment:true,name:'Agent replacement',expectedRevision:state.revision});await expect(page.locator('.draft-conflict')).toContainText('your drawing is preserved');await expect(page.getByRole('button',{name:'Publish draft as branch',exact:true})).toBeEnabled();await page.getByRole('button',{name:'Publish draft as branch',exact:true}).click();await expect(page.locator('.query-name')).toContainText('Drawing branch');
  state=(await mcp.call('wb_workspace_state')).result;expect(state.currentQuery.parentId).toBe(base.id);expect(state.currentQuery.document.idcode).not.toBe(base.document.idcode);const history=(await mcp.call('wb_query_history',{limit:20})).result;expect(history.items.some((q:any)=>q.name==='Agent replacement')).toBe(true);await expect(page.locator('.draft-conflict')).toHaveCount(0);
 }finally{mcp.close();}
});
