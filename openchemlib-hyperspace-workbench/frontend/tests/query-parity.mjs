import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {Molecule} from 'openchemlib';

const queries=JSON.parse(readFileSync(new URL('../../target/query-roundtrip.json',import.meta.url),'utf8'));
for(const query of queries) {
  const molecule=Molecule.fromIDCode(query.idcode,query.coordinates);
  assert.equal(molecule.isFragment(),query.fragment);
  assert.equal(molecule.getIDCode(),query.idcode,'OCL JS must preserve the Java query IDCode including all encoded features');
  assert.ok(molecule.getAllAtoms()>0);
}
assert.ok(queries.some(query=>{
  const molecule=Molecule.fromIDCode(query.idcode,query.coordinates);
  return Array.from({length:molecule.getAllBonds()},(_,b)=>molecule.isBondBridge(b)).some(Boolean);
}));
console.log(`Java → OCL JS IDCode/query-feature parity: ${queries.length} real OCL query documents passed.`);

const fixtures=JSON.parse(readFileSync(new URL('../../target/query-feature-contract.json',import.meta.url),'utf8'));
for(const fixture of fixtures){
  const m=Molecule.fromIDCode(fixture.document.idcode,fixture.document.coordinates);m.ensureHelperArrays(Molecule.cHelperCIP);
  assert.equal(m.getIDCode(),fixture.document.idcode,fixture.name+' exact encoded feature contract');
  for(const atom of fixture.inspection.atoms){const features=m.getAtomQueryFeaturesObject(atom.index);assert.deepEqual(Array.from(m.getAtomList(atom.index)??[]),atom.elements,fixture.name+' atom alternatives');assert.equal(features.aromatic,atom.aromatic);assert.equal(features.notAromatic,atom.notAromatic);assert.equal(features.notChain,atom.ring);assert.equal(features.matchStereo,atom.matchStereo);}
  for(const bond of fixture.inspection.bonds){assert.equal(m.getBondQueryFeatures(bond.index),bond.features,fixture.name+' bond feature bits');assert.equal(m.isBondBridge(bond.index),bond.bridge);}
}
console.log(`Java → OCL JS extended parity: ${fixtures.length} aromatic/ring/stereo/alternative/bridge queries passed.`);
