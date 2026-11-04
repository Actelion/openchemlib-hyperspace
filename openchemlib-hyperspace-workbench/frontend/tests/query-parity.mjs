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
