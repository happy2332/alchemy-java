package org.utd.cs.mln.lmap;

import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.ArrayUtils;
import org.utd.cs.mln.alchemy.util.Parser;

import java.io.FileNotFoundException;
import java.util.*;

public class MlnToHyperCube {

    /**
     */
    public String decToTernary(int d){
        if(d == 0){
            return "0";
        }
        else{
            String result = "";
            while(d > 0){
                int rem = d%3;
                d = d/3;
                result = rem+result;
            }
            return result;
        }
    }


    /**
     * Fills in a hashMap of tuples for each state(true, false, unknown) of each predicate by reading evidence
     * @param evidList List of evidence
     * @param predTuplesHashMap HashMap in which we want to put entries
     */
    public void createPredsTupleHashMap(List<FirstEvidence> evidList, HashMap<PredicateSymbol,
            ArrayList<ArrayList<TupleConstants>>> predTuplesHashMap){
        for(FirstEvidence evid : evidList){
            PredicateSymbol predSymbol = evid.symbol;
            ArrayList<Integer> constants = evid.terms;
            int truthVal = evid.truthValue;
            if(!predTuplesHashMap.containsKey(predSymbol)){
                ArrayList<ArrayList<TupleConstants>> tuples = new ArrayList<ArrayList<TupleConstants>>();
                // Create numPossibleVals+1 tuples, +1 for unknown state
                for(int i = 0 ; i < predSymbol.values.values.size()+1 ; i++){
                    tuples.add(new ArrayList<TupleConstants>());
                }
                predTuplesHashMap.put(predSymbol, tuples);
            }
            predTuplesHashMap.get(predSymbol).get(truthVal).add(new TupleConstants(constants));
        }
    }

    // It creates this list of hyperCubes for each possible value of predicate + unknown state, for each predicate, using evidence and MLN
    // To store these lists, we create a hashMap with key : predicateSymbol, value : List of size (numPossibleVals)+1. This list
    // stores list of hypercubes for each of the state of predicate.

    public HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> createPredsHyperCube(List<FirstEvidence> evidList, MLN mln, Set<String> closedworldpreds){

        // First create a similar hashMap, but instead of hyperCube it stores each tuple
        HashMap<PredicateSymbol,ArrayList<ArrayList<TupleConstants>>> predsTuplesHashMap = new HashMap<PredicateSymbol,ArrayList<ArrayList<TupleConstants>>>();
        createPredsTupleHashMap(evidList, predsTuplesHashMap);
		/*//
		for(PredicateSymbol predSymbol : predsTuplesHashMap.keySet()){
			System.out.println(predsTuplesHashMap.get(predSymbol).get(0));
		}*///

        HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubesHashMap = new HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>>();
        CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
        // Now for each predicate, create its all but unknown hypercubes
        for(PredicateSymbol predSymbol : predsTuplesHashMap.keySet()){
            int numPossibleVals = predSymbol.values.values.size();
            predsHyperCubesHashMap.put(predSymbol, new ArrayList<ArrayList<HyperCube>>());
            for(int i = 0 ; i < numPossibleVals ; i++){
                // Create hyperCubes
                ArrayList<HyperCube>hyperCubes = chcb.createHyperCubesBasic(predsTuplesHashMap.get(predSymbol).get(i));
                ArrayList<HyperCube> mergedHyperCubes = chcb.mergeHyperCubes(hyperCubes);
                ///System.out.println(hyperCubes);
                predsHyperCubesHashMap.get(predSymbol).add(mergedHyperCubes);
            }
            ///System.out.println("pred : "+predSymbol.toString());
            ///System.out.println("True False hypercubes done");
            ///System.out.println("creating unknowns");
            //System.out.println(predSymbol);
            //System.out.println(predsHyperCubesHashMap.get(predSymbol).get(0));
            // Now we have to create hypercubes for unknown. We don't need to calculate all tuples of unknown state
            // we can simply subtract true and false hypercubes from domain hyperCube <0,1,...n>

            // Create this predicate's domain hyperCube
            HyperCube predDomainHyperCube = createPredDomainHyperCube(mln, predSymbol);
            ArrayList<HyperCube>unknownHyperCubes = new ArrayList<HyperCube>();
            ArrayList<HyperCube>allButUnknownHyperCubes = new ArrayList<HyperCube>();
            // Union all but unknown hypercubes
            for (int i = 0; i < numPossibleVals; i++) {
                allButUnknownHyperCubes.addAll(predsHyperCubesHashMap.get(predSymbol).get(i));
            }

            // Merge all but unknown hypercubes so that we get efficient unknown hypercube
            allButUnknownHyperCubes = chcb.mergeHyperCubes(allButUnknownHyperCubes);
            ///System.out.println("both true and false merged and number of hypercubes in them is : "+bothTrueFalseHyperCubes.size());
            //System.out.println(bothTrueFalseHyperCubes);
            //System.out.println(predSymbol);
            // Finally subtract allbutunknown hypercubes from domainHyperCube
            unknownHyperCubes = chcb.createComplementaryHyperCubes(allButUnknownHyperCubes, predDomainHyperCube);
            ///System.out.println("Initial unmerged unknown hypercubes done");
            //System.out.println(unknownHyperCubes);
            // Merge again unknown hyperCubes to get minimal number of hypercubes
            unknownHyperCubes = chcb.mergeHyperCubes(unknownHyperCubes);
            ///System.out.println("merging of unknown hypercubes done, now size of unknown hypercubes is "+ unknownHyperCubes.size());
            //System.out.println(unknownHyperCubes);
            ///System.out.println("creating disjoint hypercubes");
            //chcb.createDisjointHyperCubes(unknownHyperCubes);
            ///System.out.println("created disjoint of unknown hypercubes and now size of unknown hypercubes is "+ unknownHyperCubes.size());
            //System.out.println(bothTrueFalseHyperCubes);
            //System.out.println(unknownHyperCubes);
            if(closedworldpreds.contains(predSymbol.symbol)){
                predsHyperCubesHashMap.get(predSymbol).get(0).addAll(unknownHyperCubes);
                predsHyperCubesHashMap.get(predSymbol).add(new ArrayList<HyperCube>());
            }
            else {
                predsHyperCubesHashMap.get(predSymbol).add(unknownHyperCubes);
            }
        }
        // If there is a predicate which didn't appear in evidence file, its unknown hypercube is only its domain
        // hypercube. Its true and false hypercubes list is empty
        for(PredicateSymbol predSymbol : mln.symbols){
            int numPossibleVals = predSymbol.values.values.size();
            if(!predsHyperCubesHashMap.containsKey(predSymbol)){
                predsHyperCubesHashMap.put(predSymbol, new ArrayList<ArrayList<HyperCube>>());
                for (int i = 0; i < numPossibleVals; i++) {
                    predsHyperCubesHashMap.get(predSymbol).add(new ArrayList<HyperCube>());
                }
                ArrayList<HyperCube> unknownHyperCubes = new ArrayList<HyperCube>();
                unknownHyperCubes.add(createPredDomainHyperCube(mln, predSymbol));
                if(closedworldpreds.contains(predSymbol.symbol)){
                    predsHyperCubesHashMap.get(predSymbol).get(0).addAll(unknownHyperCubes);
                    predsHyperCubesHashMap.get(predSymbol).add(new ArrayList<HyperCube>());
                }
                else {
                    predsHyperCubesHashMap.get(predSymbol).add(unknownHyperCubes);
                }
            }
        }
        return predsHyperCubesHashMap;
    }

    // It creates this list of hyperCubes for each state (True,False,unknown) for each predicate, using evidence and MLN
    // To store these lists, we create a hashMap with key : predicateSymbol, value : List of size 3. This list
    // of size three stores list of hypercubes for each of the state of predicate.
    // For now, we assume openworld
    public HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> createPredsHyperCubeNormalForm(ArrayList<FirstEvidence> evidList, MLN mln){

        // First create a similar hashMap, but instead of hyperCube it stores each tuple
        HashMap<PredicateSymbol,ArrayList<ArrayList<TupleConstants>>> predsTuplesHashMap = new HashMap<PredicateSymbol,ArrayList<ArrayList<TupleConstants>>>();
        createPredsTupleHashMap(evidList, predsTuplesHashMap);
			/*//
			for(PredicateSymbol predSymbol : predsTuplesHashMap.keySet()){
				System.out.println(predsTuplesHashMap.get(predSymbol).get(0));
			}*///

        HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubesHashMap = new HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>>();
        CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
        // Now for each predicate, create its True and false state's hypercubes
        for(PredicateSymbol predSymbol : predsTuplesHashMap.keySet()){
            predsHyperCubesHashMap.put(predSymbol, new ArrayList<ArrayList<HyperCube>>());
            for(int i = 0 ; i < 2 ; i++){
                // Create hyperCubes
                ArrayList<HyperCube>hyperCubes = chcb.createHyperCubesBasic(predsTuplesHashMap.get(predSymbol).get(i));
                //ArrayList<HyperCube> mergedHyperCubes = chcb.mergeHyperCubes(hyperCubes);
                ///System.out.println(hyperCubes);
                //predsHyperCubesHashMap.get(predSymbol).add(mergedHyperCubes);
                predsHyperCubesHashMap.get(predSymbol).add(hyperCubes);
            }
            ///System.out.println("pred : "+predSymbol.toString());
            ///System.out.println("True False hypercubes done");
            ///System.out.println("creating unknowns");
            //System.out.println(predSymbol);
            //System.out.println(predsHyperCubesHashMap.get(predSymbol).get(0));
            // Now we have to create hypercubes for unknown. We don't need to calculate all tuples of unknown state
            // we can simply subtract true and false hypercubes from domain hyperCube <0,1,...n>

            // Create this predicate's domain hyperCube
            HyperCube predDomainHyperCube = createPredDomainHyperCube(mln, predSymbol);
            ArrayList<HyperCube>unknownHyperCubes = new ArrayList<HyperCube>();
            ArrayList<HyperCube>bothTrueFalseHyperCubes = new ArrayList<HyperCube>();
            // Union all hypercubes of true and false
            bothTrueFalseHyperCubes.addAll(predsHyperCubesHashMap.get(predSymbol).get(0));
            bothTrueFalseHyperCubes.addAll(predsHyperCubesHashMap.get(predSymbol).get(1));
            // Merge true and false hypercubes so that we get efficient unknown hypercube
            ///bothTrueFalseHyperCubes = chcb.mergeHyperCubes(bothTrueFalseHyperCubes);
            ///System.out.println("both true and false merged and number of hypercubes in them is : "+bothTrueFalseHyperCubes.size());
            //System.out.println(bothTrueFalseHyperCubes);
            //System.out.println(predSymbol);
            // Finally subtract bothTruefalse hypercubes from domainHyperCube
            unknownHyperCubes = chcb.createComplementaryHyperCubes(bothTrueFalseHyperCubes, predDomainHyperCube);
            ///System.out.println("Initial unmerged unknown hypercubes done");
            //System.out.println(unknownHyperCubes);
            // Merge again unknown hyperCubes to get minimal number of hypercubes
            ///unknownHyperCubes = chcb.mergeHyperCubes(unknownHyperCubes);
            ///System.out.println("merging of unknown hypercubes done, now size of unknown hypercubes is "+ unknownHyperCubes.size());
            //System.out.println(unknownHyperCubes);
            ///System.out.println("creating disjoint hypercubes");
            //chcb.createDisjointHyperCubes(unknownHyperCubes);
            ///System.out.println("created disjoint of unknown hypercubes and now size of unknown hypercubes is "+ unknownHyperCubes.size());
            //System.out.println(bothTrueFalseHyperCubes);
            //System.out.println(unknownHyperCubes);
            predsHyperCubesHashMap.get(predSymbol).add(unknownHyperCubes);
        }
        // If there is a predicate which didn't appear in evidence file, its unknown hypercube is only its domain
        // hypercube. Its true and false hypercubes list is empyy
        for(PredicateSymbol predSymbol : mln.symbols){
            if(!predsHyperCubesHashMap.containsKey(predSymbol)){
                predsHyperCubesHashMap.put(predSymbol, new ArrayList<ArrayList<HyperCube>>());
                predsHyperCubesHashMap.get(predSymbol).add(new ArrayList<HyperCube>());
                predsHyperCubesHashMap.get(predSymbol).add(new ArrayList<HyperCube>());
                ArrayList<HyperCube> unknownHyperCubes = new ArrayList<HyperCube>();
                unknownHyperCubes.add(createPredDomainHyperCube(mln, predSymbol));
                predsHyperCubesHashMap.get(predSymbol).add(unknownHyperCubes);
            }
        }
        return predsHyperCubesHashMap;
    }


    public ArrayList<Formula> createFormulaHyperCube(Formula formula, HashMap<PredicateSymbol, ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap, boolean isNormal) {
        // Create set of Prod(numPossibleVals+1) formulas, where prod is over all atoms. In each formula, set each atom anypossibleValue or unknown, and then find
        // resulting formula's hyperCube. Returns resultant formulas
        ArrayList<Formula> resultFormulas = new ArrayList<Formula>();

        // First find how many potential new formulas can be created
        int numFormulas = 1;
        for(WClause clause : formula.clauses){
            for(Atom atom : clause.atoms){
                numFormulas *= atom.symbol.values.values.size()+1;
            }
        }

        // First filter out all those possible formulas in which an atom has an empty hypercube
        List<Integer> formulaIdsToConsider = new ArrayList<Integer>();
        for(int formulaId = 0 ; formulaId < numFormulas ; formulaId++){
            // Heuristic : We will loop over the clauses in order of their numAtoms i.e process shortest clause first
            // so that we can decide removal of formula as quickly as possible.
            ArrayList<Integer> clauseSizes = new ArrayList<Integer>();
            for(WClause clause : formula.clauses)
            {
                clauseSizes.add(clause.atoms.size());
            }
            int clauseIndices[] = new int[clauseSizes.size()];
            for (int i = 0; i < clauseIndices.length; i++) {
                clauseIndices[i] = i;
            }
            ArrayUtils.bubbleSort(clauseSizes,clauseIndices);
            int formulaNum = formulaId;
            boolean formulaExists = true;
            for(int clauseId : clauseIndices)
            {
                WClause clause = formula.clauses.get(clauseId);
                for(Atom atom : clause.atoms)
                {
                    int numPossibleVals = atom.symbol.values.values.size();
                    int rem = formulaNum%(numPossibleVals+1);
                    formulaNum = formulaNum/(numPossibleVals+1);
                    List<HyperCube> hypercubes = predsHyperCubeHashMap.get(atom.symbol).get(rem);
                    if(hypercubes.isEmpty())
                    {
                        formulaExists = false;
                        break;
                    }
                }
                if(formulaExists == false)
                {
                    break;
                }
            }
            if(formulaExists)
                formulaIdsToConsider.add(formulaId);
        }

        // Now loop over all potential formulas, and find which one will remain and which one will be removed.
        for(int formulaId : formulaIdsToConsider){
            // Create new formula which is a copy of the original formula
            Formula newFormula = MLN.create_new_formula(formula);
            newFormula.formulaId = formula.formulaId;

            // Stores which clauses to remove from this formula
            Set<Integer> clauseIdsToRemove = new HashSet<Integer>();

            // Need a temp variable formulaNum because we will keep dividing formulaNum
            int formulaNum = formulaId;
            // Stores truth value for each atom -> 0:true, 1:false, 2:unknown, we need a 2D list as each 1D list
            // contains atoms of a clause in the formula. Key of map is clauseId, and value is atom signs
            Map<Integer,ArrayList<Integer>> atomSigns = new HashMap<Integer,ArrayList<Integer>>();
            // Keep a boolean variable "clauseExists" which tells if a clause is to be kept or not.
            // Clause can be removed in 3 ways :
            // 1. All atoms are false (to find it, keep a boolean variable allFalseAtoms)
            // 2. Atleast an atom is true (Make clauseExists=false as soon as we find it)
            // 3. Atleast an atom has an empty hypercube (Make clauseExists=false as soon as we find it)
            // After the loop of atoms, if we find allFalseAtoms = true, then make clauseExists = false.

            boolean formulaExists = true;
            boolean allClausesTrue = true;
            // In conditions 1 and 3 above, formula also doesn't exist, so make formulaExists = false also
            // For condition 2, we can make formulaExists false only when all clauses are true i.e. if
            // we haven't added any clause in the formula

            // Heuristic : We will loop over the clauses in order of their numAtoms i.e process shortest clause first
            // so that we can decide removal of formula as quickly as possible.
            ArrayList<Integer> clauseSizes = new ArrayList<Integer>();
            for(WClause clause : formula.clauses)
            {
                clauseSizes.add(clause.atoms.size());
            }
            int clauseIndices[] = new int[clauseSizes.size()];
            for (int i = 0; i < clauseIndices.length; i++) {
                clauseIndices[i] = i;
            }
            ArrayUtils.bubbleSort(clauseSizes,clauseIndices);
            for(int clauseIndex : clauseIndices){
                WClause newClause = newFormula.clauses.get(clauseIndex);
                // Stores truth value for each atom -> 0:true, 1:false, 2:unknown
                ArrayList<Integer> clauseAtomSigns = new ArrayList<>();
                atomSigns.put(clauseIndex, clauseAtomSigns);
                boolean allFalseAtoms = true;
                boolean clauseExists = true;
                // Loop over all the atoms
                for(int atomIndex = 0 ; atomIndex < newClause.atoms.size() ; atomIndex++){
                    Atom atom = newClause.atoms.get(atomIndex);
                    int numPossibleVals = atom.symbol.values.values.size();
                    int rem = formulaNum%(numPossibleVals+1);
                    formulaNum = formulaNum/(numPossibleVals+1);
                    // Definitely unknown
                    if(rem == numPossibleVals){ // last value is unknown so we are checking with numPossibleVals
                        clauseAtomSigns.add(2);
                        allFalseAtoms = false;
                    }

                    // If sign of this atom and equality of rem with valTrue of this atom are opposite,
                    // then this literal is true. Ex : S(x)=5 (so here valTrue is 5), suppose rem is 5, now sign is false and
                    // rem is equal to valTrue(which is 5), and so they are of opposite signs and so atom is true.
                    // Ex : !S(x)=5, suppose rem is 2, then sign is true but equality is false, so again literal becomes true.

                    // Clause Definitely True so don't add into the formula, and break out of atoms loop
                    else if((newClause.sign.get(atomIndex) ^ (rem == newClause.valTrue.get(atomIndex))) == true){
                        clauseExists = false;
                        allFalseAtoms = false;
                        clauseAtomSigns.add(0);
                        // break; there is no break here, we have to set hypercube of every atom
                    }
                    // Literal Definitely False
                    else{
                        // At this point don't remove this atom from the clause, because we have to create union
                        // of hypercubes later.
                        // Ex : S(x) | C(x), suppose hypercube of S(x) is {1,2}, and of C(x) is {2,3},
                        // and suppose C(x) is false, then if we remove C(x) here only, then we are left with
                        // clause S(x) with hypercube of S(x) being {1,2}, which is wrong.
                        // Instead we should have clause S(x) with hypercube being {2}.
                        // so for now, just store which atom is false, so we can remove them later.
                        clauseAtomSigns.add(1);
                    }

                    // Now set atom's hyperCube according to
                    setAtomHyperCubes(atom, predsHyperCubeHashMap.get(atom.symbol).get(rem));
//                    for(HyperCube hyperCube : predsHyperCubeHashMap.get(atom.symbol).get(rem)){
//                        atom.hyperCubes.add(new HyperCube(hyperCube));
//                    }
                    // If there is no hyperCube for this atom, clause will be invalid
                    // This condition is redundant as we have already neglected those formulas in which
                    // empty hypercube atoms occur
                    if(atom.hyperCubes.size() == 0){
                        clauseExists = false;
                        formulaExists = false;
                        break;
                    }

                }
                if(allFalseAtoms == true){
                    clauseExists = false;
                    formulaExists = false;
                }
                // If we already know formula is false, break from the clauses loop
                if(formulaExists == false)
                    break;
                // If this clause doesn't exist, then move to next clause in the formula
                else if(clauseExists == false)
                {
                    clauseIdsToRemove.add(clauseIndex);
                    continue;
                }

                // At this point, this clause definitely exist, so make allClausesTrue false
                allClausesTrue = false;

            }
            // If this formula doesn't exist or allclauses are true, move to next formula
            if(formulaExists == false || allClausesTrue == true)
                continue;

            // At this point, this formula definitely exist
            // The only possibility of it not in existence is due to no common hypercube constant in same term
            // Ex : S(x)=1 ^ C(x)=1, suppose hypercube of S(x) is {1,2}, and of C(x) is {3,4}, then
            // This formula will not exist because there is no common constant for x.
            // Now to get single union of hyperCube of formula, loop over every 2 atoms, and join them accordingly
            // Ex : S(x,y) | C(y,z) ^ R(u,z) will give P(x,y,z,u)
            // atom1 will be formula's built hypercube till now. Initially set it to first atom

            // We definitely know that there is atleast one clause in this formula, because otherwise
            // formulaExists would have been false and we would have continued the loop above.

            Atom atom1 = MLN.create_new_atom(newFormula.clauses.get(0).atoms.get(0));

            // check if atom's hypercube is empty

            atom1.terms = newFormula.clauses.get(0).atoms.get(0).terms;
            // finalAtom will be result of merging two atoms each time
            Atom finalAtom = new Atom();
            // If there is only 1 atom, no need to merge, just return as it is
            if(newFormula.clauses.size() == 1 && newFormula.clauses.get(0).atoms.size()==1){
                finalAtom.terms = atom1.terms;
                finalAtom.hyperCubes = atom1.hyperCubes;
            }
            else
            {
                // Loop over atoms, starting from 2nd atom i.e. merging 1st and 2nd atoms
                for(int clauseId = 0 ; clauseId < newFormula.clauses.size() ; clauseId++)
                {
                    WClause newClause = newFormula.clauses.get(clauseId);
                    for(int atomId = 0 ; atomId < newClause.atoms.size() ; atomId++)
                    {
                        // skip first atom of first clause because that is already in atom1
                        if(clauseId == 0 && atomId == 0)
                            continue;
                        Atom atom = newClause.atoms.get(atomId);
                        int atom1NumTerms = atom1.terms.size();
                        // Arrays to store indices of terms common in both atoms
                        ArrayList<Integer>atom1CommonTermIndices = new ArrayList<Integer>();
                        ArrayList<Integer>atomCommonTermIndices = new ArrayList<Integer>();
                        ///System.out.println("atom1 terms : " + atom1.terms);
                        ///System.out.println("atom terms : " + atom.terms);
                        //Find indices of common terms in both atoms : atom1 and atom
                        for(int atom1TermId = 0 ; atom1TermId < atom1NumTerms ; atom1TermId++){
                            for(int atomTermId = 0 ; atomTermId < atom.terms.size() ; atomTermId++){
                                if(atom1.terms.get(atom1TermId) == (atom.terms.get(atomTermId))){
                                    atom1CommonTermIndices.add(atom1TermId);
                                    atomCommonTermIndices.add(atomTermId);
                                }
                            }
                        }

                        // Finally join two atoms and get result in finalAtom
                        joinAtomsCommonVars(atom1, atom, atom1CommonTermIndices, atomCommonTermIndices, finalAtom);
                        // check if finalatom's hypercube is not empty
                        if(finalAtom.hyperCubes.isEmpty())
                        {
                            formulaExists = false;
                            break;
                        }

                        if(isNormal == false){
                            CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
                            finalAtom.hyperCubes = chcb.mergeHyperCubes(finalAtom.hyperCubes);
                        }
                        // Now for next iteration, atom1 will be result of this iteration
                        atom1 = finalAtom;
                        // If it is last iteration, don't empty out finalAtom because we are setting clause's hypercube acc to it.
                        if(!(clauseId == newFormula.clauses.size()-1 && atomId == newClause.atoms.size()-1))
                            finalAtom = new Atom();
                    }
                    if(formulaExists == false)
                        break;

                }

            }


            if(formulaExists == false)
                continue;

            // Now we can safely remove false atoms from the formula
            // Also due to removal of atoms, some terms may get removed, so create a list of terms which
            // we have to keep
            // Also Remove all clauses from this new formula which do not exist
            // start removing from last so that indexing of starting ones don't get affected
            // Also remove corresponding atomSigns
            Set<Term> termsToStay = new HashSet<Term>();
            for (int clauseId = newFormula.clauses.size()-1; clauseId >= 0; clauseId--) {
                WClause newClause = newFormula.clauses.get(clauseId);
                if(clauseIdsToRemove.contains(clauseId))
                {
                    newFormula.clauses.remove(clauseId);
                    atomSigns.remove(clauseId);
                    continue;
                }
                int numAtoms = newClause.atoms.size();
                for (int atomId = numAtoms-1; atomId >=0 ; atomId--) {
                    if(atomSigns.get(clauseId).get(atomId) == 1)
                    {
                        newClause.removeAtom(atomId);
                    }
                    else
                    {
                        for(Term term : newClause.atoms.get(atomId).terms)
                        {
                            termsToStay.add(term);
                        }
                    }
                }
            }

            newFormula.terms = finalAtom.terms;

            // Now set root's hypercubeList
            for(HyperCube hc : finalAtom.hyperCubes){
                newFormula.root.hyperCubesList.add(hc);
            }

            // Find the termIndices for which we have to remove terms from formulaHyperCubes
            for(int termId = newFormula.terms.size() - 1 ; termId >= 0 ; termId--){
                if(!termsToStay.contains(newFormula.terms.get(termId))){
                    newFormula.terms.remove(termId);
                    for(HyperCube hyperCube : newFormula.root.hyperCubesList)
                    {
                        hyperCube.num_copies *= hyperCube.varConstants.get(termId).size();
//                        if(hyperCube.num_copies > 1)
//                        {
//                            System.out.println("numcopies changed");
//                        }
                        hyperCube.varConstants.remove(termId);
                    }

                    }
            }

            //CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
            //newFormula.root.hyperCubesList = chcb.mergeHyperCubes(newFormula.root.hyperCubesList);

            if(newFormula.root.hyperCubesList.size() > 0){
                newFormula.parentFormulaId = formula.formulaId;
                resultFormulas.add(newFormula);
            }
        }
        return resultFormulas;
    }

    private void setAtomHyperCubes(Atom atom, ArrayList<HyperCube> hyperCubes) {
        for(HyperCube hc : hyperCubes)
        {
            HyperCube newHc = new HyperCube();
            List<Term> termsSeen = new ArrayList<>();
            boolean isEmptyHyperCube = false;
            for(int termIndex = 0 ; termIndex < atom.terms.size() ; termIndex++)
            {
                Term term = atom.terms.get(termIndex);
                Set<Integer> termDomain = new HashSet<>(term.domain);
                termDomain.retainAll(hc.varConstants.get(termIndex));
                if(termDomain.isEmpty())
                {
                    isEmptyHyperCube = true;
                    break;
                }
                int oldTermIndex = termsSeen.indexOf(term);
                if(oldTermIndex == -1)
                {
                    termsSeen.add(term);
                    newHc.varConstants.add(termDomain);
                }
                else
                {
                    newHc.varConstants.get(oldTermIndex).retainAll(termDomain);
                    if(newHc.varConstants.get(oldTermIndex).isEmpty())
                    {
                        isEmptyHyperCube = true;
                        break;
                    }
                }
            }
            if(!isEmptyHyperCube)
            {
                atom.hyperCubes.add(newHc);
            }
        }

    }
    /*public ArrayList<WClause> createClauseHyperCube(WClause clause, HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap, boolean isNormal){
        // Create set of 3^(#atoms) clauses and in each clause, set each atom true, false or unknown, and then find
        // resulting clause's hyperCube. Returns resultant clause
        //System.out.println(predsHyperCubeHashMap.get(clause.atoms.get(0).symbol).get(1));
        int numAtoms = clause.atoms.size();
        int numClauses = (int)Math.pow(3, numAtoms);
        ArrayList<WClause> resultClauses = new ArrayList<WClause>();
        // Iterate over all possible clause and generate its hyperCube
        for(int clauseId = 0 ; clauseId < numClauses ; clauseId++){
            ///System.out.println("clauseId = "+clauseId);
            WClause newClause = MLN.create_new_clause(clause);
            Set<Term> termsToStay = new HashSet<Term>();
            // To get all combinations, convert clauseId to ternary number where left changes fastest. Ex : 8 is 220 for 3 atoms case
            // It means first and second atoms are in unknown state, and third atom is true, although it can become
            //false due to its negation in clause
            int clauseNum = clauseId;
            int rem = clauseNum%3;
            boolean booleanRem = (rem == 1); // booleanRem is true only when we are to pick false hypercubes of pred
            ArrayList<Integer> atomSigns = new ArrayList<Integer>(); // Stores sign for each atom -> 0:True, 1:False, 2:Unknown
            boolean clauseExists = true;
            boolean allFalseAtoms = true;
            for(int atomIndex = 0 ; atomIndex < numAtoms ; atomIndex++){
                Atom atom = newClause.atoms.get(atomIndex);
                // Definitely Unknown
                if(rem == 2){
                    atomSigns.add(2);
                    allFalseAtoms = false;
                    for(Term term : atom.terms)
                        termsToStay.add(term);
                }
                // Definitely True
                else if(booleanRem == newClause.sign.get(atomIndex)){
                    atomSigns.add(0);
                    newClause.satisfied = true;
                }
                // Definitely False
                else{
                    atomSigns.add(1);
                }
                // Now set atom's hyperCube according to
                for(HyperCube hyperCube : predsHyperCubeHashMap.get(atom.symbol).get(rem)){
                    atom.hyperCubes.add(new HyperCube(hyperCube));
                }
                ///System.out.println("Printing " + atom.symbol + "'s hyperCubes : ");
                ///System.out.println(atom.hyperCubes);
                // If there is no hyperCube for this atom, clause will be invalid
                if(atom.hyperCubes.size() == 0){
                    clauseExists = false;
                    break;
                }
                clauseNum = clauseNum/3;
                rem = clauseNum%3;
                booleanRem = (rem == 1);
            }
            if(clauseExists == false || newClause.satisfied || allFalseAtoms == true){
                continue;
            }
            //System.out.println("atomSigns = ");
            //System.out.println(atomSigns);
            // Now to get single union of hyperCube of clause, loop over every 2 atoms, and join them accordingly
            // Ex : S(x,y) | C(y,z) will give P(x,y,z)
            // atom1 will be clause's built hypercube till now. Initially set it to first atom
            Atom atom1 = MLN.create_new_atom(newClause.atoms.get(0));
            atom1.terms = newClause.atoms.get(0).terms;
            // finalAtom will be result of merging two atoms each time
            Atom finalAtom = new Atom();
            // If there is only 1 atom, no need to merge, just return as it is
            if(numAtoms == 1){
                finalAtom.terms = atom1.terms;
                finalAtom.hyperCubes = atom1.hyperCubes;
            }
            // Loop over atoms, starting from 2nd atom i.e. merging 1st and 2nd atoms
            for(int atomId = 1 ; atomId < numAtoms ; atomId++){
                Atom atom = newClause.atoms.get(atomId);
                int atom1NumTerms = atom1.terms.size();
                // Arrays to store indices of terms common in both atoms
                ArrayList<Integer>atom1CommonTermIndices = new ArrayList<Integer>();
                ArrayList<Integer>atomCommonTermIndices = new ArrayList<Integer>();
                ///System.out.println("atom1 terms : " + atom1.terms);
                ///System.out.println("atom terms : " + atom.terms);
                //Find indices of common terms in both atoms : atom1 and atom
                for(int atom1TermId = 0 ; atom1TermId < atom1NumTerms ; atom1TermId++){
                    for(int atomTermId = 0 ; atomTermId < atom.terms.size() ; atomTermId++){
                        if(atom1.terms.get(atom1TermId) == (atom.terms.get(atomTermId))){
                            atom1CommonTermIndices.add(atom1TermId);
                            atomCommonTermIndices.add(atomTermId);
                        }
                    }
                }
                // Finally join two atoms and get result in finalAtom
                joinAtomsCommonVars(atom1, atom, atom1CommonTermIndices, atomCommonTermIndices, finalAtom);
                if(isNormal == false){
                    CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
                    finalAtom.hyperCubes = chcb.mergeHyperCubes(finalAtom.hyperCubes);
                }
                ///System.out.println("After merging first " + (int)(atomId+1) + " atoms, finalAtom becomes : " + finalAtom.hyperCubes.toString());
                // Now for next iteration, atom1 will be result of this iteration
                atom1 = finalAtom;
                // If it is last iteration, don't empty out finalAtom because we are setting clause's hypercube acc to it.
                if(atomId < numAtoms-1)
                    finalAtom = new Atom();
            }
            if(clauseExists == false){
                continue;
            }
            // Remove false atoms from this clause
            for(int atomId = numAtoms - 1 ; atomId >= 0 ; atomId--){
                if(atomSigns.get(atomId) == 1){
                    newClause.removeAtom(atomId);
                }
            }
            newClause.terms = finalAtom.terms;
            for(HyperCube hc : finalAtom.hyperCubes){
                ArrayList<HyperCube> hcList = new ArrayList<HyperCube>();
                hcList.add(hc);
                newClause.root.hyperCubesList.add(hcList);
            }

            // Find the termIndices for which we have to remove terms from clauseHyperCubes
            for(int termId = newClause.terms.size() - 1 ; termId >= 0 ; termId--){
                if(!termsToStay.contains(newClause.terms.get(termId))){
                    newClause.terms.remove(termId);
                    for(int phId = 0 ; phId < newClause.root.hyperCubesList.size() ; phId++){
                        for(HyperCube hyperCube : newClause.root.hyperCubesList.get(phId))
                            hyperCube.varConstants.remove(termId);
                    }

                }
            }
            if(newClause.root.hyperCubesList.size() > 0){
                resultClauses.add(newClause);
            }
            ///System.out.println("result clause : ");
            ///newClause.print();
            ///System.out.println("result Clause's termList : " +newClause.terms);
            ///System.out.println("result Clause's hyperCubes : " + newClause.hyperCubes);
        }
        return resultClauses;
    }*/

	/*
	public ArrayList<WClause> createClauseHyperCubeOld(WClause clause, HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap){
		// Create set of 3^(#atoms) clauses and in each clause, set each atom true, false or unknown, and set its
		// terms domain accordingly
		//System.out.println(predsHyperCubeHashMap.get(clause.atoms.get(0).symbol).get(1));
		int numAtoms = clause.atoms.size();
		int numClauses = (int)Math.pow(3, numAtoms);
		ArrayList<WClause> resultClauses = new ArrayList<WClause>();
		for(int clauseId = 0 ; clauseId < numClauses ; clauseId++){
			System.out.println("clauseId = "+clauseId);
			WClause newClause = MLN.create_new_clause(clause);
			Set<Term> distinctTerms = new HashSet<Term>();
			int clauseNum = clauseId;
			int rem = clauseNum%3;
			boolean booleanRem = (rem == 1);
			ArrayList<Integer> atomSigns = new ArrayList<Integer>();
			boolean clauseExists = true;
			for(int atomIndex = 0 ; atomIndex < numAtoms ; atomIndex++){
				Atom atom = newClause.atoms.get(atomIndex);
				for(Term term : atom.terms){
					distinctTerms.add(term);
				}
				if(rem == 2){
					atomSigns.add(2);

				}
				else if(booleanRem == newClause.sign.get(atomIndex)){
					atomSigns.add(0);
					newClause.satisfied = true;
				}
				else{
					atomSigns.add(1);
				}
				for(HyperCube hyperCube : predsHyperCubeHashMap.get(atom.symbol).get(atomSigns.get(atomIndex))){
					atom.hyperCubes.add(new HyperCube(hyperCube));
				}
				System.out.println("Printing " + atom.symbol + "'s hyperCubes : ");
				System.out.println(atom.hyperCubes);
				if(atom.hyperCubes.size() == 0){
					clauseExists = false;
					break;
				}
				clauseNum = clauseNum/3;
				rem = clauseNum%3;
				booleanRem = (rem == 1);
			}
			if(clauseExists == false){
				continue;
			}
			System.out.println("atomSigns = ");
			System.out.println(atomSigns);
			System.out.println("No. of distinct terms : "+distinctTerms.size());
			int termNum = 0;
			for(Term term : distinctTerms){
				System.out.println("term number : " + (++termNum));
				Set<Integer> intersectionSetTerm = new HashSet<Integer>(term.domain);
				for(int atomIndex = 0 ; atomIndex < numAtoms ; atomIndex++){
					Atom atom = newClause.atoms.get(atomIndex);
					int termIndex = atom.terms.indexOf(term);
					if(termIndex >= 0){
						intersectionSetTerm.retainAll(createSetUnionTerm(predsHyperCubeHashMap.get(atom.symbol).get(atomSigns.get(atomIndex)),termIndex));
					}
				}
				if(intersectionSetTerm.size() == 0){
					clauseExists = false;
					break;
				}
				System.out.println("Term = "+term.toString()+" Its intersectionSetTerm is : ");
				System.out.println(intersectionSetTerm);
				for(int atomIndex = 0 ; atomIndex < numAtoms ; atomIndex++){
					Atom atom = newClause.atoms.get(atomIndex);
					int termIndex = atom.terms.indexOf(term);
					if(termIndex >= 0){
						findIntersectionPredHyperCubeTerm(atom.hyperCubes,intersectionSetTerm, termIndex);
						System.out.println("Now Printing " + atom.symbol + "'s hyperCubes : ");
						System.out.println(atom.hyperCubes);
					}
				}
			}
			if(clauseExists == false){
				continue;
			}
			resultClauses.add(newClause);
			//System.out.println(predsHyperCubeHashMap.get(newClause.atoms.get(2).symbol).get(0));
		}
		return resultClauses;
	}
	*/

    private void joinAtomsCommonVars(Atom atom1, Atom atom2,
                                     ArrayList<Integer> atom1CommonTermIndices,
                                     ArrayList<Integer> atom2CommonTermIndices, Atom finalAtom) {

        // Set resulting terms : the order will be terms appearing in atom1, then terms
        // appearing only in atom2
        for(int termId = 0 ; termId < atom1.terms.size() ; termId++){
            finalAtom.terms.add(atom1.terms.get(termId));
        }

        for(int termId = 0 ; termId < atom2.terms.size() ; termId++){
            if(!atom2CommonTermIndices.contains(termId)){
                finalAtom.terms.add(atom2.terms.get(termId));
            }
        }
        // Loop over atom1's commonTerms segments and for each segment set (of common terms), find the segment sets
        // in atom2's commonterm Segment sets with which it has non zero intersection, and for intersecting set
        // of segments, join atom1 and atom2
        for(int atom1HcId = 0 ; atom1HcId < atom1.hyperCubes.size() ; atom1HcId++){
            for(int atom2HcId = 0 ; atom2HcId < atom2.hyperCubes.size() ; atom2HcId++){
                HyperCube finalHyperCube = new HyperCube(atom1.hyperCubes.get(atom1HcId));
                boolean isEmptyHyperCube = false;
                for(int commonTermIdIndex = 0 ; commonTermIdIndex < atom1CommonTermIndices.size() ; commonTermIdIndex++){
                    int atom1TermId = atom1CommonTermIndices.get(commonTermIdIndex);
                    int atom2TermId = atom2CommonTermIndices.get(commonTermIdIndex);
                    finalHyperCube.varConstants.get(atom1TermId).retainAll(atom2.hyperCubes.get(atom2HcId).varConstants.get(atom2TermId));
                    if(finalHyperCube.varConstants.get(atom1TermId).size() == 0){
                        isEmptyHyperCube = true;
                        break;
                    }
                }
                if(!isEmptyHyperCube){
                    for(int termId = 0 ; termId < atom2.terms.size() ; termId++){
                        if(!atom2CommonTermIndices.contains(termId)){
                            finalHyperCube.varConstants.add(new HashSet<Integer>(atom2.hyperCubes.get(atom2HcId).varConstants.get(termId)));
                        }
                    }
                    finalAtom.hyperCubes.add(finalHyperCube);
                }
            }
        }
        //CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
        //finalAtom.hyperCubes = chcb.mergeHyperCubes(finalAtom.hyperCubes);
    }

    /*
        // Joins two atoms on the basis of common terms and return resultant atom
        // Algo : First take hyperCube list of only common terms in both atoms. Ex : S(x,z,y) | C(y,u,v), then take
        // hypercube list of only variable y in S and C, and take union of both lists.
        // Now make the hyperCubes in this list as disjoint or identical
        // Now for each hyperCube in this disjoint list, see which hypercubes of first atom's y's hypercubes does it
        // intersect with, and make a list of those hyperCubes of S's <x,z>'s hyperCubes
        // Similarly make a list of C's <u,v> hyperCubes
        // Finally take cross product of both these lists of <x,z> and <u,v>. Also while taking cross product, sandwitch
        // y in between them so hypercube becomes <x,z,y,u,v>
        private void joinAtomsCommonVarsOld(Atom atom1, Atom atom2,
                ArrayList<Integer> atom1CommonTermIndices,
                ArrayList<Integer> atom2CommonTermIndices, Atom finalAtom) {

            // Set resulting terms : the order will be terms appearing only in atom1, then common terms, then terms
            // appearing only in atom2
            for(int termId = 0 ; termId < atom1.terms.size() ; termId++){
                if(!atom1CommonTermIndices.contains(termId)){
                    finalAtom.terms.add(atom1.terms.get(termId));
                }
            }

            for(Integer termId : atom2CommonTermIndices){
                finalAtom.terms.add(atom2.terms.get(termId));
            }

            for(int termId = 0 ; termId < atom2.terms.size() ; termId++){
                if(!atom2CommonTermIndices.contains(termId)){
                    finalAtom.terms.add(atom2.terms.get(termId));
                }
            }
            //System.out.println(finalAtom.terms);
            // It will contain union of hypercubes of common term
            Set<HyperCube> commonTermHyperCubes = new HashSet<HyperCube>();
            // It will contain hypercubes of common term in atom1 only
            ArrayList<HyperCube> atom1CommonTermHyperCubes = new ArrayList<HyperCube>();
            // It will contain hyperCubes of nonCommon terms in atom1 only
            ArrayList<HyperCube> atom1NonCommonTermHyperCubes = new ArrayList<HyperCube>();
            int numHyperCubes = atom1.hyperCubes.size();
            // Go over hyperCubes of atom1
            for(int hcId = 0 ; hcId < numHyperCubes ; hcId++){
                HyperCube nonCommonHyperCube = new HyperCube();
                //HyperCube hyperCube = atom1.hyperCubes.get(hcId);
                HyperCube commonTermHyperCube = new HyperCube();
                // Go over common terms and add their hyperCube into commonTermHyperCube
                for(Integer termId : atom1CommonTermIndices){
                    commonTermHyperCube.varConstants.add(new HashSet<Integer>(atom1.hyperCubes.get(hcId).varConstants.get((int)termId)));
                }
                // Add rest terms' hyperCubes into nonCommonTermHyperCube
                for(int termId = 0 ; termId < atom1.terms.size() ; termId++){
                    if(!atom1CommonTermIndices.contains(termId))
                        nonCommonHyperCube.varConstants.add(new HashSet<Integer>(atom1.hyperCubes.get(hcId).varConstants.get(termId)));
                }
                commonTermHyperCubes.add(commonTermHyperCube);
                atom1NonCommonTermHyperCubes.add(nonCommonHyperCube);
                System.out.println("Non common Term hypercube of atom 1 : " + nonCommonHyperCube.toString());
                System.out.println("common Term hypercube of atom 1 : " + commonTermHyperCube.toString());
                atom1CommonTermHyperCubes.add(commonTermHyperCube);
            }
            // Now for atom2
            ArrayList<HyperCube> atom2CommonTermHyperCubes = new ArrayList<HyperCube>();
            ArrayList<HyperCube> atom2NonCommonTermHyperCubes = new ArrayList<HyperCube>();
            numHyperCubes = atom2.hyperCubes.size();
            for(int hcId = 0 ; hcId < numHyperCubes ; hcId++){
                HyperCube nonCommonHyperCube = new HyperCube();
                //HyperCube hyperCube = atom1.hyperCubes.get(hcId);
                HyperCube commonTermHyperCube = new HyperCube();
                for(Integer termId : atom2CommonTermIndices){
                    commonTermHyperCube.varConstants.add(new HashSet<Integer>(atom2.hyperCubes.get(hcId).varConstants.get((int)termId)));
                }
                for(int termId = 0 ; termId < atom2.terms.size() ; termId++){
                    if(!atom2CommonTermIndices.contains(termId))
                        nonCommonHyperCube.varConstants.add(new HashSet<Integer>(atom2.hyperCubes.get(hcId).varConstants.get(termId)));
                }
                commonTermHyperCubes.add(commonTermHyperCube);
                atom2NonCommonTermHyperCubes.add(nonCommonHyperCube);
                System.out.println("Non common Term hypercube of atom 2 : " + nonCommonHyperCube.toString());
                System.out.println("common Term hypercube of atom 2 : " + commonTermHyperCube.toString());
                atom2CommonTermHyperCubes.add(commonTermHyperCube);
            }
            System.out.println("Common term hyperCubes set : ");
            System.out.println(commonTermHyperCubes);
            // Now store this final set into an arrayList
            ArrayList<HyperCube> commonTermHyperCubeList = new ArrayList<HyperCube>(commonTermHyperCubes);
            // Now make this final list as disjoint or identical hyperCubes
            CreateHyperCubeBasic chcb = new CreateHyperCubeBasic();
            chcb.createDisjointHyperCubes(commonTermHyperCubeList);
            System.out.println("Common term disjoint hypercubes : " + commonTermHyperCubeList.toString());
            //For each common HyperCube in commonTermHyperCubeList, find the list of hypercubes in atom1 which
            //intersect with it. Do this for atom2 also. Then take cross product of both lists obtained.
            for(HyperCube commonTermHyperCube : commonTermHyperCubeList){
                numHyperCubes = atom1NonCommonTermHyperCubes.size();
                ArrayList<HyperCube> atom1NonCommonIntersectHyperCubes = new ArrayList<HyperCube>();
                // Go over each hyperCube of atom1
                for(int hcId = 0 ; hcId < numHyperCubes ; hcId++){
                    // If this common hyperCube was part of any hyperCube of atom1's common term hyperCubes
                    if(commonTermHyperCube.hasIntersection(atom1CommonTermHyperCubes.get(hcId))){
                        atom1NonCommonIntersectHyperCubes.add(new HyperCube(atom1NonCommonTermHyperCubes.get(hcId)));
                    }
                }
                // If this hyperCube has no intersection with any hyperCube of atom1, then just continue
                if(atom1NonCommonIntersectHyperCubes.size() == 0)
                    continue;
                // Now do same for atom2
                numHyperCubes = atom2NonCommonTermHyperCubes.size();
                ArrayList<HyperCube> atom2NonCommonIntersectHyperCubes = new ArrayList<HyperCube>();
                for(int hcId = 0 ; hcId < numHyperCubes ; hcId++){
                    if(commonTermHyperCube.hasIntersection(atom2CommonTermHyperCubes.get(hcId))){
                        atom2NonCommonIntersectHyperCubes.add(new HyperCube(atom2NonCommonTermHyperCubes.get(hcId)));
                    }
                }
                if(atom2NonCommonIntersectHyperCubes.size() == 0)
                    continue;
                // Take cross product of both lists
                for(int atom1HcId = 0 ; atom1HcId < atom1NonCommonIntersectHyperCubes.size() ; atom1HcId++){
                    for(int atom2HcId = 0 ; atom2HcId < atom2NonCommonIntersectHyperCubes.size() ; atom2HcId++){
                        // Initialize resulting hyperCube with atom1's nonCommon hyperCube
                        HyperCube hyperCube = new HyperCube(atom1NonCommonIntersectHyperCubes.get(atom1HcId));
                        System.out.println("HyperCube becomes : "+hyperCube.toString());
                        int numCommonVars = atom1CommonTermIndices.size();
                        // Now add hyperCube of common variables
                        for(int varId = 0 ; varId < numCommonVars ; varId++){
                            hyperCube.varConstants.add(new HashSet<Integer>(commonTermHyperCube.varConstants.get(varId)));
                        }
                        System.out.println("HyperCube becomes : "+hyperCube.toString());
                        // Now add hyperCube for atom2's nonCommon hyperCube
                        int numAtom2Vars = atom2NonCommonIntersectHyperCubes.get(0).getVarCount();
                        for(int varId = 0 ; varId < numAtom2Vars ; varId++){
                            hyperCube.varConstants.add(new HashSet<Integer>(atom2NonCommonIntersectHyperCubes.get(atom2HcId).varConstants.get(varId)));
                        }
                        System.out.println("hyperCube becomes : "+hyperCube.toString());
                        finalAtom.hyperCubes.add(hyperCube);
                        System.out.println("final atom hyperCubes becomes : " + finalAtom.hyperCubes);
                    }
                }
            }
        }

        // This method takes a list of set of integers intersectionSetTerm, and finds intersection with the segments
        //in the input list of hypercubes to create new segments and thus modifying the list of hyperCubes.
        // This list of hyperCubes is of a predicate in the clause. Also, intersection is done only on variable termIndex.
        //passed
        private void findIntersectionPredHyperCubeTerm(
                ArrayList<HyperCube> hyperCubes, Set<Integer> intersectionSetTerm, int termIndex) {
            for(int hcId = hyperCubes.size() - 1 ; hcId >= 0 ; hcId--){
                Set<Integer>intersectionSet = new HashSet<Integer>(hyperCubes.get(hcId).varConstants.get(termIndex));
                intersectionSet.retainAll(intersectionSetTerm);
                if(intersectionSet.size() > 0){
                    hyperCubes.get(hcId).setVarConstants(intersectionSet, termIndex);
                }
                else{
                    hyperCubes.remove(hcId);
                }
            }
        }

        private Set<Integer> createSetUnionTerm(ArrayList<HyperCube> inpHyperCubes, int termId) {
            Set<Integer> unionSet = new HashSet<Integer>();
            for(HyperCube hyperCube : inpHyperCubes){
                unionSet.addAll(hyperCube.varConstants.get(termId));
            }
            return unionSet;
        }

        private void setAtomDomain(
                Atom atom,
                HashMap<PredicateSymbol, ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap,
                int sign) {
            PredicateSymbol predSymbol = atom.symbol;
            int numTerms = atom.terms.size();
            int numHyperCubes = predsHyperCubeHashMap.get(predSymbol).get(sign).size();
            for(int hyperCubeIndex = 0 ; hyperCubeIndex < numHyperCubes ; hyperCubeIndex++){
                for(int termIndex = 0 ; termIndex < numTerms ; termIndex++){
                    Term term = atom.terms.get(termIndex);
                    Set<Integer> setToAdd = predsHyperCubeHashMap.get(predSymbol).get(sign).get(hyperCubeIndex).varConstants.get(termIndex);
                    System.out.println(setToAdd);
                    term.segmentDomain.add(new HashSet<Integer>(setToAdd));
                }
            }
        }
        */
    public ArrayList<Set<Integer>> intersectHyperCubeVar(ArrayList<Set<Integer>> list1, ArrayList<Set<Integer>> list2){
        ArrayList<Set<Integer>> result = new ArrayList<Set<Integer>>();
        for(int i = 0 ; i < list1.size() ; i++){
            Set<Integer> intersection = new HashSet<Integer>(list1.get(i));
            intersection.retainAll(list2.get(i));
            result.add(intersection);
        }
        return result;
    }

    private HyperCube createPredDomainHyperCube(MLN mln,
                                                PredicateSymbol predSymbol) {
        HyperCube hyperCube = new HyperCube();
        for(String varType : predSymbol.variable_types)
        {
            Set<Integer> domain = new HashSet<>(mln.varTypeToDomainMap.get(varType));
            hyperCube.varConstants.add(domain);
        }
        //System.out.println(mln.predicateDomainMap.keySet());
//        ArrayList<Integer> domainIndexList = mln.predicateDomainMap.get(predSymbol.id);
//        int varCnt = domainIndexList.size();
//        for(int varId = 0 ; varId < varCnt ; varId++){
//            ArrayList<String> domainValues = (ArrayList<String>) mln.domainList.get(domainIndexList.get(varId)).values;
//            Set<Integer> domain = new HashSet<Integer>();
//            for(String val : domainValues){
//                domain.add(Integer.parseInt(val));
//            }
//            hyperCube.varConstants.add(domain);
//        }
        return hyperCube;
    }

    public static void main(String[] args) throws FileNotFoundException, PredicateNotFound {
        // TODO Auto-generated method stub
        MLN mln = new MLN(false);
        Parser parser = new Parser(mln);
        String filename = new String("smoke/smoke_mln_29.txt");
        //String filename = new String("entity_resolution/er-bnct-eclipse.mln");
        parser.parseInputMLNFile(filename);

        ArrayList<FirstEvidence> evidList = parser.parseInputEvidenceFile("smoke/smoke_evidence.txt");
        //ArrayList<Evidence> evidList = parser.parseInputEvidenceFile("entity_resolution/er-test-eclipse.db");
        MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
        HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList,mln,new HashSet<String>());


        for(PredicateSymbol predSymbol : predsHyperCubeHashMap.keySet()){
            System.out.println(predSymbol.toString()+predsHyperCubeHashMap.get(predSymbol).get(0).toString()+" "+predsHyperCubeHashMap.get(predSymbol).get(1).toString() + " "+predsHyperCubeHashMap.get(predSymbol).get(2).toString());
        }
        int origNumClauses = mln.clauses.size();
        boolean isNormal = false;
//        for(int clauseId = 0 ; clauseId < origNumClauses ; clauseId++){
//            mln.clauses.addAll(mlnToHyperCube.createClauseHyperCube(mln.clauses.get(clauseId), predsHyperCubeHashMap, isNormal));
//        }
        for(int clauseId = origNumClauses-1 ; clauseId >= 0 ; clauseId--){
            mln.clauses.remove(clauseId);
        }
//        ArrayList<Set<Pair>> predEquivalenceClasses = new ArrayList<Set<Pair>>();
//        Decomposer d = new Decomposer();
//        d.findEquivalenceClasses(mln, predEquivalenceClasses);
//        System.out.println("Equivalence Classes : " + predEquivalenceClasses);
//        System.out.println(predEquivalenceClasses);
//        System.out.println(d.isDecomposer(predEquivalenceClasses.get(0),mln));
//        ArrayList<Integer> sizeOfSegments = new ArrayList<>();
//        if(d.isDecomposer(predEquivalenceClasses.get(0),mln))
//            d.initDecomposer(mln, predEquivalenceClasses.get(0), sizeOfSegments);

        //System.out.println(predsHyperCubeHashMap);
        //System.out.println(mlnToHyperCube.createPredDomainHyperCube(mln, mln.symbols.get(2)));
        //System.out.println(mln.domainList);
    }


}
