package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class QueryFeatureContractTest {
    record Fixture(String name,QueryChemistry.Document document,QueryChemistry.Inspection inspection,List<String> examples,List<Boolean> expected){}
    private QueryChemistry.Edit edit(QueryChemistry.Operation op,String ref,String state,Boolean enabled,List<String> types,Boolean exclude,List<Integer> elements){return new QueryChemistry.Edit(op,List.of(ref),elements,exclude,state,types,enabled,null,null);}
    private Fixture fixture(String name,String smiles,List<QueryChemistry.Edit> edits,List<String> examples,List<Boolean> expected)throws Exception {
        var document=QueryChemistry.edit("source",QueryChemistry.importStructure("smiles",smiles,true,null),edits);
        assertEquals(expected,QueryChemistry.test(document,examples).stream().map(QueryChemistry.ExampleResult::matches).toList(),name);
        return new Fixture(name,document,QueryChemistry.inspect(name,document),examples,expected);
    }
    @Test void typedAlternativesRingAromaticStereoAndBridgeFeaturesMatchJavaAndGenerateJsFixtures()throws Exception {
        List<Fixture> fixtures=new ArrayList<>();
        fixtures.add(fixture("excluded elements","C",List.of(edit(QueryChemistry.Operation.ATOM_ELEMENTS,"source:a0",null,null,null,true,List.of(7,8))),List.of("C","N","O","F"),List.of(true,false,false,true)));
        fixtures.add(fixture("aromatic atom","C",List.of(edit(QueryChemistry.Operation.ATOM_AROMATIC,"source:a0","aromatic",null,null,null,null)),List.of("c1ccccc1","CC"),List.of(true,false)));
        fixtures.add(fixture("ring atom","C",List.of(edit(QueryChemistry.Operation.ATOM_RING,"source:a0","ring",null,null,null,null)),List.of("C1CCCCC1","CCC"),List.of(true,false)));
        fixtures.add(fixture("chain atom","C",List.of(edit(QueryChemistry.Operation.ATOM_RING,"source:a0","chain",null,null,null,null)),List.of("C1CCCCC1","CCC"),List.of(false,true)));
        fixtures.add(fixture("bond alternatives","CC",List.of(edit(QueryChemistry.Operation.BOND_TYPES,"source:b0",null,null,List.of("single","double"),null,null)),List.of("CC","C=C","C#C"),List.of(true,true,false)));
        fixtures.add(fixture("ring bond","CC",List.of(edit(QueryChemistry.Operation.BOND_RING,"source:b0","ring",null,null,null,null)),List.of("C1CCCCC1","CCC"),List.of(true,false)));
        fixtures.add(fixture("chain bond","CC",List.of(edit(QueryChemistry.Operation.BOND_RING,"source:b0","chain",null,null,null,null)),List.of("C1CCCCC1","CCC"),List.of(false,true)));
        fixtures.add(fixture("aromatic bond","c1ccccc1",List.of(edit(QueryChemistry.Operation.BOND_AROMATIC,"source:b0","aromatic",null,null,null,null)),List.of("c1ccccc1","C1CCCCC1"),List.of(true,false)));
        var atomStereo=QueryChemistry.importStructure("smiles","F[C@](Cl)(Br)I",true,null);var am=QueryChemistry.parse(atomStereo);int center=-1;for(int a=0;a<am.getAtoms();a++)if(am.getAtomicNo(a)==6)center=a;
        fixtures.add(fixture("atom stereochemistry","F[C@](Cl)(Br)I",List.of(edit(QueryChemistry.Operation.ATOM_STEREO,"source:a"+center,null,true,null,null,null)),List.of("F[C@](Cl)(Br)I","F[C@@](Cl)(Br)I"),List.of(true,false)));
        var bondStereo=QueryChemistry.importStructure("smiles","C/C=C/C",true,null);var bm=QueryChemistry.parse(bondStereo);int doubleBond=-1;for(int b=0;b<bm.getBonds();b++)if(bm.getBondOrder(b)==2)doubleBond=b;
        fixtures.add(fixture("bond stereochemistry","C/C=C/C",List.of(edit(QueryChemistry.Operation.BOND_STEREO,"source:b"+doubleBond,null,true,null,null,null)),List.of("C/C=C/C","C/C=C\\C"),List.of(true,false)));
        var bridge=QueryChemistry.edit("source",QueryChemistry.importStructure("smiles","CO",true,null),List.of(new QueryChemistry.Edit(QueryChemistry.Operation.BOND_BRIDGE,List.of("source:b0"),null,null,null,null,null,0,2)));
        fixtures.add(new Fixture("zero-minimum bridge",bridge,QueryChemistry.inspect("bridge",bridge),List.of("CO","CCO","O"),List.of(true,true,false)));
        Files.createDirectories(Path.of("target"));WorkspaceStore.JSON.writeValue(Path.of("target/query-feature-contract.json").toFile(),fixtures);
    }
}
