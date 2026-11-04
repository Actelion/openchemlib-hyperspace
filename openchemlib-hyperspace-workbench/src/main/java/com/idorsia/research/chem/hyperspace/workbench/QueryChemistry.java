package com.idorsia.research.chem.hyperspace.workbench;

import com.actelion.research.chem.*;
import com.idorsia.research.chem.hyperspace.SubstructureSearchHelper;
import java.util.*;
import static com.idorsia.research.chem.hyperspace.workbench.WorkbenchException.require;

/** Actual OCL query semantics; all references belong to one immutable revision. */
public final class QueryChemistry {
    public record Document(String idcode,String coordinates,boolean fragment){}
    public record Atom(String ref,int index,String element,int charge,String featuresHex,List<Integer> elements,
                       boolean aromatic,boolean notAromatic,boolean ring,boolean matchStereo){}
    public record Bond(String ref,int index,String atomA,String atomB,int type,int features,
                       boolean bridge,Integer minAtoms,Integer maxAtoms){}
    public record Inspection(String revisionId,Document document,List<Atom> atoms,List<Bond> bonds,List<String> warnings){}
    public enum Operation {ATOM_ELEMENTS,ATOM_AROMATIC,ATOM_RING,ATOM_STEREO,BOND_TYPES,BOND_AROMATIC,BOND_RING,BOND_STEREO,BOND_BRIDGE,CLEAR_BRIDGE}
    public record Edit(Operation operation,List<String> refs,List<Integer> elements,Boolean exclude,
                       String state,List<String> types,Boolean enabled,Integer minAtoms,Integer maxAtoms){}
    public record ExampleResult(String structure,boolean matches){}
    private QueryChemistry(){}

    public static Document importStructure(String format,String value,Boolean fragment,String coordinates)throws Exception {
        require(value!=null&&!value.isBlank()&&value.length()<=1_000_000,"A bounded nonempty structure is required");
        StereoMolecule molecule=new StereoMolecule();
        switch(format==null?"smiles":format) {
            case "smiles" -> new SmilesParser().parse(molecule,value);
            case "idcode" -> new IDCodeParser().parse(molecule,value,coordinates);
            case "mol" -> require(new MolfileParser().parse(molecule,value),"Invalid molfile");
            case "smarts" -> new SmilesParser(SmilesParser.SMARTS_MODE_IS_SMARTS).parse(molecule,value);
            default -> throw new WorkbenchException("invalid_format",400,"Expected smiles, idcode, mol or smarts");
        }
        if(fragment!=null)molecule.setFragment(fragment);
        return canonical(molecule);
    }
    public static StereoMolecule parse(Document document) {
        StereoMolecule molecule=new StereoMolecule();new IDCodeParser().parse(molecule,document.idcode(),document.coordinates());
        molecule.ensureHelperArrays(Molecule.cHelperCIP);return molecule;
    }
    public static Document canonical(StereoMolecule molecule) {
        molecule.ensureHelperArrays(Molecule.cHelperCIP);
        require(molecule.getAtoms()>0&&molecule.getAtoms()<=256,"Query must contain 1–256 heavy atoms");
        var canonizer=new Canonizer(molecule);
        return new Document(canonizer.getIDCode(),canonizer.getEncodedCoordinates(),molecule.isFragment());
    }
    public static Inspection inspect(String revision,Document document) {
        StereoMolecule m=parse(document);List<Atom> atoms=new ArrayList<>();List<Bond> bonds=new ArrayList<>();List<String> warnings=new ArrayList<>();
        for(int a=0;a<m.getAtoms();a++) {
            long features=m.getAtomQueryFeatures(a);int[] alternatives=m.getAtomList(a);
            atoms.add(new Atom(atomRef(revision,a),a,m.getAtomLabel(a),m.getAtomCharge(a),Long.toHexString(features),
                    alternatives==null?List.of():Arrays.stream(alternatives).boxed().toList(),
                    (features&Molecule.cAtomQFAromatic)!=0,(features&Molecule.cAtomQFNotAromatic)!=0,
                    (features&Molecule.cAtomQFNotChain)!=0,(features&Molecule.cAtomQFMatchStereo)!=0));
        }
        int bridges=0;
        for(int b=0;b<m.getBonds();b++) {
            boolean bridge=m.isBondBridge(b);if(bridge)bridges++;
            bonds.add(new Bond(bondRef(revision,b),b,atomRef(revision,m.getBondAtom(0,b)),atomRef(revision,m.getBondAtom(1,b)),
                    m.getBondType(b),m.getBondQueryFeatures(b),bridge,bridge?m.getBondBridgeMinSize(b):null,bridge?m.getBondBridgeMaxSize(b):null));
        }
        if(!m.isFragment())warnings.add("Concrete molecule: substructure execution will use a fragment copy. Learned scoring requires a concrete supported molecule.");
        if(bridges>2)warnings.add("Current Hyperspace bridge expansion supports at most two bridge bonds.");
        return new Inspection(revision,document,atoms,bonds,warnings);
    }
    public static Document edit(String revision,Document document,List<Edit> edits) {
        require(edits!=null&&!edits.isEmpty()&&edits.size()<=128,"Supply 1–128 typed feature edits");
        StereoMolecule m=parse(document);m.setFragment(true);
        for(Edit edit:edits) {
            require(edit.operation()!=null&&edit.refs()!=null&&!edit.refs().isEmpty()&&edit.refs().size()<=256,"Each edit requires an operation and references");
            boolean atom=edit.operation().name().startsWith("ATOM_");
            for(String ref:edit.refs()) {
                int index=reference(revision,ref,atom?"a":"b",atom?m.getAtoms():m.getBonds());
                switch(edit.operation()) {
                    case ATOM_ELEMENTS -> {
                        require(edit.elements()!=null&&!edit.elements().isEmpty()&&edit.elements().size()<=32,"Supply 1–32 atomic numbers");
                        require(edit.elements().stream().allMatch(z->z!=null&&z>0&&z<=118),"Invalid element atomic number");
                        m.setAtomList(index,edit.elements().stream().mapToInt(Integer::intValue).distinct().sorted().toArray(),Boolean.TRUE.equals(edit.exclude()));
                    }
                    case ATOM_AROMATIC -> {
                        long flag=switch(state(edit)){case "any"->0;case "aromatic"->Molecule.cAtomQFAromatic;case "aliphatic"->Molecule.cAtomQFNotAromatic;default->throw invalidState();};
                        m.setAtomQueryFeature(index,Molecule.cAtomQFAromState,false);if(flag!=0)m.setAtomQueryFeature(index,flag,true);
                    }
                    case ATOM_RING -> {
                        long flag=switch(state(edit)){case "any"->0;case "ring"->Molecule.cAtomQFNotChain;case "chain"->Molecule.cAtomQFNot2RingBonds|Molecule.cAtomQFNot3RingBonds|Molecule.cAtomQFNot4RingBonds;default->throw invalidState();};
                        m.setAtomQueryFeature(index,Molecule.cAtomQFRingState,false);if(flag!=0)m.setAtomQueryFeature(index,flag,true);
                    }
                    case ATOM_STEREO -> m.setAtomQueryFeature(index,Molecule.cAtomQFMatchStereo,Boolean.TRUE.equals(edit.enabled()));
                    case BOND_TYPES -> {
                        require(edit.types()!=null&&!edit.types().isEmpty(),"Supply permitted bond types");int mask=0;
                        for(String type:edit.types())mask|=switch(type){case "single"->Molecule.cBondTypeSingle;case "double"->Molecule.cBondTypeDouble;case "triple"->Molecule.cBondTypeTriple;case "aromatic"->Molecule.cBondTypeDelocalized;default->throw new WorkbenchException("invalid_bond_type",400,"Unknown bond type");};
                        m.setBondQueryFeature(index,Molecule.cBondQFAllBondTypes,false);m.setBondQueryFeature(index,mask,true);
                    }
                    case BOND_AROMATIC -> {
                        int flag=switch(state(edit)){case "any"->0;case "aromatic"->Molecule.cBondQFAromatic;case "aliphatic"->Molecule.cBondQFNotAromatic;default->throw invalidState();};
                        m.setBondQueryFeature(index,Molecule.cBondQFAromState,false);if(flag!=0)m.setBondQueryFeature(index,flag,true);
                    }
                    case BOND_RING -> {
                        int flag=switch(state(edit)){case "any"->0;case "ring"->Molecule.cBondQFRing;case "chain"->Molecule.cBondQFNotRing;default->throw invalidState();};
                        m.setBondQueryFeature(index,Molecule.cBondQFRingState,false);if(flag!=0)m.setBondQueryFeature(index,flag,true);
                    }
                    case BOND_STEREO -> m.setBondQueryFeature(index,Molecule.cBondQFMatchStereo,Boolean.TRUE.equals(edit.enabled()));
                    case BOND_BRIDGE -> {
                        require(edit.minAtoms()!=null&&edit.maxAtoms()!=null,"Bridge requires minimum and maximum intervening atoms");
                        int min=edit.minAtoms(),max=edit.maxAtoms();
                        require(min>=0&&min<=15&&max>=min&&max-min<=15,"OCL bridge min is 0–15 and span max-min is 0–15");
                        require(min!=0||max!=0,"OCL 0..0 encodes zero bits, meaning no bridge; a 0..0 bridge cannot be represented");
                        m.setBondQueryFeature(index,Molecule.cBondQFBridge,false);
                        m.setBondQueryFeature(index,(min<<Molecule.cBondQFBridgeMinShift)|((max-min)<<Molecule.cBondQFBridgeSpanShift),true);
                    }
                    case CLEAR_BRIDGE -> m.setBondQueryFeature(index,Molecule.cBondQFBridge,false);
                }
            }
        }
        return canonical(m);
    }
    public static List<ExampleResult> test(Document document,List<String> smiles)throws Exception {
        require(smiles!=null&&smiles.size()<=50,"At most 50 example molecules");
        StereoMolecule fragment=parse(document);fragment.setFragment(true);
        List<StereoMolecule> queries=SubstructureSearchHelper.expandBridgedSearches(fragment);
        List<ExampleResult> results=new ArrayList<>();
        for(String text:smiles) {
            StereoMolecule molecule=new StereoMolecule();new SmilesParser().parse(molecule,text);boolean found=false;
            for(StereoMolecule q:queries){q.setFragment(true);var searcher=new SSSearcher();searcher.setMol(q,molecule);if(searcher.isFragmentInMolecule()){found=true;break;}}
            results.add(new ExampleResult(text,found));
        }
        return results;
    }
    public static String atomRef(String revision,int index){return revision+":a"+index;}
    public static String bondRef(String revision,int index){return revision+":b"+index;}
    private static int reference(String revision,String ref,String type,int size) {
        String prefix=revision+":"+type;
        if(ref==null||!ref.startsWith(prefix))throw new WorkbenchException("stale_reference",409,"Atom/bond references must belong to the edited query revision");
        try{int index=Integer.parseInt(ref.substring(prefix.length()));require(index>=0&&index<size,"Atom/bond index outside query");return index;}
        catch(NumberFormatException e){throw new WorkbenchException("invalid_reference",400,"Invalid atom/bond reference");}
    }
    private static String state(Edit edit){return edit.state()==null?"any":edit.state();}
    private static WorkbenchException invalidState(){return new WorkbenchException("invalid_feature_state",400,"Unsupported query-feature state");}
}
