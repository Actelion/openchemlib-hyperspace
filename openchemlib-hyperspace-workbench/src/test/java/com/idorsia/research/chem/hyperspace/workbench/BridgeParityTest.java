package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class BridgeParityTest {
    private QueryChemistry.Edit bridge(int index,int min,int max){return new QueryChemistry.Edit(QueryChemistry.Operation.BOND_BRIDGE,List.of("base:b"+index),null,null,null,null,null,min,max);}
    private List<Boolean> nativeMatches(QueryChemistry.Document document,List<String> examples)throws Exception {
        List<Boolean> out=new ArrayList<>();for(String example:examples){var m=new StereoMolecule();new SmilesParser().parse(m,example);var s=new SSSearcher();s.setMol(QueryChemistry.parse(document),m);out.add(s.isFragmentInMolecule());}return out;
    }
    private List<Boolean> expandedMatches(QueryChemistry.Document document,List<String> examples)throws Exception {
        var expanded=com.idorsia.research.chem.hyperspace.SubstructureSearchHelper.expandBridgedSearches(QueryChemistry.parse(document));
        List<Boolean> out=new ArrayList<>();for(String example:examples){var m=new StereoMolecule();new SmilesParser().parse(m,example);boolean found=false;
            for(var q:expanded){var s=new SSSearcher();s.setMol(q,m);if(s.isFragmentInMolecule()){found=true;break;}}out.add(found);}return out;
    }
    @Test void zeroIntermediateAtomsKeepsBothBridgeEndpoints()throws Exception {
        var base=QueryChemistry.importStructure("smiles","CO",true,null);var query=QueryChemistry.edit("base",base,List.of(bridge(0,0,2)));
        var examples=List.of("O","CO","C=O","CCO","CCCO","CCCCO","CC");var nativeResult=nativeMatches(query,examples);
        assertEquals(List.of(false,true,true,true,true,true,false),nativeResult);
        var candidateResult=expandedMatches(query,examples);
        assertFalse(candidateResult.get(0),"Zero-length expansion must not remove the carbon endpoint");
        var expanded=com.idorsia.research.chem.hyperspace.SubstructureSearchHelper.expandBridgedSearches(QueryChemistry.parse(query));
        assertEquals(2,expanded.get(0).getAtoms());
        for(int i=0;i<nativeResult.size();i++)if(nativeResult.get(i))assertTrue(candidateResult.get(i),"Expanded paths must cover native positives; extra candidates are checked after assembly");
        assertEquals(nativeResult,QueryChemistry.test(query,examples).stream().map(QueryChemistry.ExampleResult::matches).toList());
    }
    @Test void twoBridgesReferToOriginalBondsAfterExpansion()throws Exception {
        var base=QueryChemistry.importStructure("smiles","CON",true,null);var query=QueryChemistry.edit("base",base,List.of(bridge(0,1,2),bridge(1,0,2)));
        var examples=List.of("CON","CCON","COCN","CCOCN","CCOCCN","CCCOCCN","CCCOCCCN","CCCCOCCCN","O","CN");
        assertEquals(nativeMatches(query,examples),expandedMatches(query,examples));
        assertEquals(nativeMatches(query,examples),QueryChemistry.test(query,examples).stream().map(QueryChemistry.ExampleResult::matches).toList());
    }
    @Test void concreteProductsUseNativeShortestPathBridgeValidation()throws Exception {
        var base=QueryChemistry.importStructure("smiles","CO",true,null);var query=QueryChemistry.edit("base",base,List.of(bridge(0,2,3)));
        var examples=List.of("C1COC1","C1CCOC1","CCCO");
        assertEquals(nativeMatches(query,examples),QueryChemistry.test(query,examples).stream().map(QueryChemistry.ExampleResult::matches).toList());
        for(String example:examples){var molecule=new StereoMolecule();new SmilesParser().parse(molecule,example);
            var hit=ResultService.makeHit("test",example,"asset",null,null,List.of(),molecule,query,false,Map.of(),Map.of());
            assertEquals(nativeMatches(query,List.of(example)).get(0),!hit.matchedAtoms().isEmpty());}
    }
}
