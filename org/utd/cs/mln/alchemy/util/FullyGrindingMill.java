package org.utd.cs.mln.alchemy.util;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.inference.InferTest;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Created by Happy on 2/28/17.
 */
public class FullyGrindingMill {
    public static boolean fgmdebug = true;
    public static int permutationCount;
    private GroundMLN groundMln;
    private Map<GroundFormula,Integer> groundFormulaToIndexMap = new HashMap<GroundFormula,Integer>();

    private void init(List<GroundPredicate> groundPredicateList, Map<GroundPredicate, Integer> groundPredicateIntegerMap) {
        groundMln = new GroundMLN();
        for (int i = 0; i < groundPredicateList.size(); i++) {
            groundMln.indexToGroundPredMap.put(i, groundPredicateList.get(i));
        }
        groundMln.groundPredToIntegerMap = groundPredicateIntegerMap;
    }

    public GroundMLN ground(MLN mln, List<GroundPredicate> groundPredicateList, Map<GroundPredicate, Integer> groundPredicateIntegerMap) {
        init(groundPredicateList, groundPredicateIntegerMap);
        int formulaNum = 1;
        for(Formula formula : mln.formulas)
        {
            Set<Term> formulaWiseTermToGround = new HashSet<Term>();
            for (WClause clause : formula.clauses) {
                for (Atom atom : clause.atoms) {
                    for (int j = 0; j < atom.terms.size(); j++) {
                        Term term = atom.terms.get(j);
                        formulaWiseTermToGround.add(term);
                    }
                }
            }
            if(fgmdebug) {
                System.out.println("Grounding formula "+formulaNum+" out of " + mln.formulas.size());
            }
            ground(formula, new ArrayList<Term>(formulaWiseTermToGround));
            formulaNum++;
        }

        for(PredicateSymbol symbol : mln.symbols)
        {
            groundMln.symbols.add(new GroundPredicateSymbol(symbol.id, symbol.symbol, symbol.values, symbol.variable_types));
        }
        return groundMln;
    }

    private void ground(Formula formula, ArrayList<Term> terms) {
        int[][] permutations = permute(terms);
        if(fgmdebug) {
            System.out.println("permutations length : " + permutations.length);
        }

        for(int i = 0 ; i < permutations.length ; i++)
        {
            if(fgmdebug) {
                if(i%1000 == 0)
                    System.out.println("\t"+i+" done\n");
            }
            GroundFormula newFormula = new GroundFormula();
            int currentFormulaId = groundMln.groundFormulas.size();
            newFormula.formulaId = currentFormulaId;
            newFormula.parentFormulaId.add(formula.formulaId);
            newFormula.numCopies.add(1);
            newFormula.weight = new LogDouble(formula.weight.getValue(), true);
            newFormula.originalWeight = new LogDouble(formula.weight.getValue(), true);
            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
            for(int c = 0 ; c < formula.clauses.size() ; c++)
            {
                WClause clause = formula.clauses.get(c);
                GroundClause newGroundClause = new GroundClause();
                newGroundClause.formulaId = currentFormulaId;
                newGroundClause.weight = new LogDouble(clause.weight.getValue(), true);
                Map<Integer, BitSet> gpIndexToSatVals = new HashMap<>();
                List<GroundPredicate> newGroundPreds = new ArrayList<>(); // We need this list because once a groundClause is created, we want to
                // update formulaIds info of each groundPred. We can't do it on the go because we don't know whether groundClause will be created
                // or not, since it can be removed due to preprocessing.
                boolean clauseToRemove = false;

                // Iterate over each first order atom, and create ground atom for it.
                for(int j = 0 ; j < clause.atoms.size() ; j++)
                {
                    boolean sign = clause.sign.get(j);
                    Atom oldAtom = clause.atoms.get(j); // first order atom
                    int valTrue = clause.valTrue.get(j);
                    GroundPredicate gp = new GroundPredicate(); // GroundPredicate to create
                    gp.symbol = new GroundPredicateSymbol(oldAtom.symbol.id,oldAtom.symbol.symbol,oldAtom.symbol.values, oldAtom.symbol.variable_types);
                    // Fill in the terms with constants
                    for(Term term : oldAtom.terms)
                    {
                        int termIndex = terms.indexOf(term);
                        gp.terms.add(permutations[i][termIndex]);
                    }

                    // Find the index of this groundPredicate in groundPredicate List
                    MyAssert.assume(groundMln.groundPredToIntegerMap.containsKey(gp));
                    int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
                    gp = groundMln.indexToGroundPredMap.get(gpIndex);

                    // Check if this groundPredicate occurs first time in this ground clause. then update
                    // groundClause's data structures about this groundPredicate.
                    int gpIndexInClause = newGroundClause.groundPredIndices.indexOf(gpIndex);
                    //GroundAtom newGroundAtom = new GroundAtom(gpIndex, gpIndexInClause, valTrue, sign);
                    if(gpIndexInClause == -1)
                    {
                        newGroundPreds.add(gp);
                        newGroundClause.groundPredIndices.add(gpIndex);
                        gpIndexInClause = newGroundClause.groundPredIndices.size()-1;
                        newGroundClause.globalToLocalPredIndex.put(gpIndex,gpIndexInClause);
                        gpIndexToSatVals.put(gpIndexInClause, new BitSet(gp.numPossibleValues));
                    }

                    // Now once we have added new ground Atom, we need to check if ground clause gets satisfied or not.
                    BitSet gpBitSet = new BitSet(gp.numPossibleValues);
                    gpBitSet.set(valTrue);
                    if(sign == true)
                        gpBitSet.flip(0,gp.numPossibleValues);
                    gpBitSet.or(gpIndexToSatVals.get(gpIndexInClause));

                    // If all bits are set for this groundPred, then this clause will always be satisfied and hence,
                    // shouldn't be added into groundformula. Note that, although at this point, we know that
                    // this clause shouldn't be added, but still we shouldn't just break out of this loop, as we
                    // need to add groundPredicates, but we shouldn't add any clauseInfo into groundPredicates appearing
                    // in this clause.
//                    if(gpBitSet.cardinality() == gp.numPossibleValues)
//                        clauseToRemove = true;
                    gpIndexToSatVals.put(gpIndexInClause, gpBitSet);
                }

                // If this clause is to be added, then only update all gp's formulaId's info
                if(clauseToRemove == false)
                {
                    for(int gpId = 0 ; gpId < newGroundPreds.size() ; gpId++)
                    {
                        BitSet b = gpIndexToSatVals.get(gpId);
                        newGroundClause.grounPredBitSet.add(b);
                    }
                    newGroundClauseList.add(newGroundClause);
                    for(GroundPredicate gp : newGroundPreds)
                    {
                        int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
                        newFormula.groundPredIndices.add(gpIndex);
                        if(!gp.groundFormulaIds.containsKey(currentFormulaId))
                        {
                            gp.groundFormulaIds.put(currentFormulaId, new HashSet<Integer>());
                        }
                        gp.groundFormulaIds.get(currentFormulaId).add(newGroundClauseList.size()-1);
                    }
                }
            }
            if(newGroundClauseList.size() > 0)
            {
                newFormula.groundClauses.addAll(newGroundClauseList);
                groundMln.groundFormulas.add(newFormula);
            }
        }
    }

    public GroundMLN groundWithHyperCubes(MLN mln, List<GroundPredicate> groundPredicateList, Map<GroundPredicate, Integer> groundPredicateIntegerMap) {
        init(groundPredicateList, groundPredicateIntegerMap);
        int formulaNum = 1;
        for(Formula formula : mln.formulas)
        {
//            Set<Term> formulaWiseTermToGround = new HashSet<Term>();
//            for (WClause clause : formula.clauses) {
//                for (Atom atom : clause.atoms) {
//                    for (int j = 0; j < atom.terms.size(); j++) {
//                        Term term = atom.terms.get(j);
//                        formulaWiseTermToGround.add(term);
//                    }
//                }
//            }
//            if(fgmdebug) {
//                System.out.println("Grounding formula "+formulaNum+" out of " + mln.formulas.size());
//            }
//            ground(formula, new ArrayList<Term>(formulaWiseTermToGround));
//            ArrayList<Integer> termIndices = new ArrayList<Integer>();
//            for (int i = 0; i < formula.terms.size(); i++) {
//                termIndices.add(i);
//            }


            if(fgmdebug) {
                System.out.println("Grounding formula "+formulaNum+" out of " + mln.formulas.size());
                System.out.println(formula);
            }
            groundFormulaWithHypercubes(formula, formula.terms);
            formulaNum++;
        }

        for(PredicateSymbol symbol : mln.symbols)
        {
            groundMln.symbols.add(new GroundPredicateSymbol(symbol.id, symbol.symbol, symbol.values, symbol.variable_types));
        }

        return groundMln;
    }

    private void groundFormulaWithHypercubes(Formula formula, List<Term> terms)
    {
        for(HyperCube hc : formula.root.hyperCubesList)
        {
            int[][] permutations = permuteHyperCube(hc);
            groundWithPermute(formula, terms, permutations, hc.num_copies);
        }
        System.gc();
    }

    private void groundWithPermute(Formula formula, List<Term> terms, int permutations[][], int numCopies)
    {
        //System.out.println("Number of permutations : "+permutations.length);
        for(int i = 0 ; i < permutations.length ; i++)
        {
            if(fgmdebug) {
                if((i+1)%100000 == 0)
                    System.out.println("\t"+i+" done\n");
            }

            GroundFormula newFormula = new GroundFormula();
            int currentFormulaId = groundMln.groundFormulas.size();
            newFormula.formulaId = currentFormulaId;
            MyAssert.assume(formula.parentFormulaId != -1);
            newFormula.parentFormulaId.add(formula.parentFormulaId);
            newFormula.numCopies.add(numCopies);
            newFormula.weight = new LogDouble(numCopies * formula.weight.getValue(), true);
            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
            for(int c = 0 ; c < formula.clauses.size() ; c++)
            {
                WClause clause = formula.clauses.get(c);
                GroundClause newGroundClause = new GroundClause();
                newGroundClause.formulaId = currentFormulaId;
                newGroundClause.weight = new LogDouble(clause.weight.getValue(), true);
                Map<Integer, BitSet> gpIndexToSatVals = new HashMap<>();
                List<GroundPredicate> allNewGroundPreds = new ArrayList<>();
                List<GroundPredicate> newGroundPreds = new ArrayList<>(); // We need this list because once a groundClause is created, we want to
                // update formulaIds info of each groundPred. We can't do it on the go because we don't know whether groundClause will be created
                // or not, since it can be removed due to preprocessing.
                boolean clauseToRemove = false;

                // Iterate over each first order atom, and create ground atom for it.
                for(int j = 0 ; j < clause.atoms.size() ; j++)
                {
                    boolean sign = clause.sign.get(j);
                    Atom oldAtom = clause.atoms.get(j); // first order atom
                    int valTrue = clause.valTrue.get(j);
                    GroundPredicate gp = new GroundPredicate(); // GroundPredicate to create
                    gp.symbol = new GroundPredicateSymbol(oldAtom.symbol.id,oldAtom.symbol.symbol,oldAtom.symbol.values, oldAtom.symbol.variable_types);
                    // Fill in the terms with constants
                    for(Term term : oldAtom.terms)
                    {
                        int termIndex = terms.indexOf(term);
                        gp.terms.add(permutations[i][termIndex]);
                    }

                    // Find index of this gp in groundPredicateList
                    MyAssert.assume(groundMln.groundPredToIntegerMap.containsKey(gp));
                    int gpIndex = groundMln.groundPredToIntegerMap.get(gp);

                    gp = groundMln.indexToGroundPredMap.get(gpIndex);

                    // Check if this groundPredicate occurs first time in this ground clause. then update
                    // groundClause's data structures about this groundPredicate.
                    int gpIndexInClause = newGroundClause.groundPredIndices.indexOf(gpIndex);
                    allNewGroundPreds.add(gp);
                    //GroundAtom newGroundAtom = new GroundAtom(gpIndex, gpIndexInClause, valTrue, sign);
                    if(gpIndexInClause == -1)
                    {
                        newGroundPreds.add(gp);
                        newGroundClause.groundPredIndices.add(gpIndex);
                        gpIndexInClause = newGroundClause.groundPredIndices.size()-1;
                        newGroundClause.globalToLocalPredIndex.put(gpIndex,gpIndexInClause);
                        gpIndexToSatVals.put(gpIndexInClause, new BitSet(gp.numPossibleValues));
                    }

                    // Now once we have added new ground Atom, we need to check if ground clause gets satisfied or not.
                    BitSet gpBitSet = new BitSet(gp.numPossibleValues);
                    gpBitSet.set(valTrue);
                    if(sign == true)
                        gpBitSet.flip(0,gp.numPossibleValues);
                    gpBitSet.or(gpIndexToSatVals.get(gpIndexInClause));

                    // If all bits are set for this groundPred, then this clause will always be satisfied and hence,
                    // shouldn't be added into groundformula. Note that, although at this point, we know that
                    // this clause shouldn't be added, but still we shouldn't just break out of this loop, as we
                    // need to add groundPredicates, but we shouldn't add any clauseInfo into groundPredicates appearing
                    // in this clause.
                    if(gpBitSet.cardinality() == gp.numPossibleValues)
                        clauseToRemove = true;
                    gpIndexToSatVals.put(gpIndexInClause, gpBitSet);

                }

                // If this clause is to be added, then only update all gp's formulaId's info
                if(clauseToRemove == false)
                {
                    for(int gpId = 0 ; gpId < newGroundPreds.size() ; gpId++)
                    {
                        BitSet b = gpIndexToSatVals.get(gpId);
                        newGroundClause.grounPredBitSet.add(b);
                    }
                    newGroundClauseList.add(newGroundClause);
                    for(GroundPredicate gp : newGroundPreds)
                    {
                        int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
                        newFormula.groundPredIndices.add(gpIndex);
                        if(!gp.groundFormulaIds.containsKey(currentFormulaId))
                        {
                            gp.groundFormulaIds.put(currentFormulaId, new HashSet<Integer>());
                        }
                        gp.groundFormulaIds.get(currentFormulaId).add(newGroundClauseList.size()-1);
                    }
                    for(GroundPredicate gp : allNewGroundPreds)
                    {
                        int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
                        newFormula.allGroundPredIndices.add(gpIndex);
                    }
                }
            }
            if(newGroundClauseList.size() > 0)
            {
                newFormula.groundClauses.addAll(newGroundClauseList);
                // check if this ground formula also exist
                if(groundFormulaToIndexMap.containsKey(newFormula))
                {
                    int oldFormulaIndex = groundFormulaToIndexMap.get(newFormula);

                    for(Integer gpIndex : newFormula.groundPredIndices)
                    {
                        GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpIndex);
                        gp.groundFormulaIds.remove(currentFormulaId);
                    }
                    GroundFormula oldGroundFormula = groundMln.groundFormulas.get(oldFormulaIndex);
                    oldGroundFormula.weight = oldGroundFormula.weight.multiply(newFormula.weight);
                    oldGroundFormula.parentFormulaId.add(formula.formulaId);
                    oldGroundFormula.numCopies.add(numCopies);

                }

                else {
                    groundFormulaToIndexMap.put(newFormula, groundMln.groundFormulas.size());
                    groundMln.groundFormulas.add(newFormula);
                }
            }
        }
    }

//    private void ground(Formula formula, ArrayList<Term> terms) {
//        int[][] permutations = permute(terms);
//        if(fgmdebug) {
//            System.out.println("permutations length : " + permutations.length);
//        }
//
//        groundWithPermute(formula, terms, permutations);
//    }


    public static int[][] permuteHyperCube(HyperCube hyperCube){
        ArrayList<Set<Integer>> sets = hyperCube.varConstants;
        int numTerms = sets.size();
        int permutationSize = 1;
        for (Set<Integer> set : hyperCube.varConstants) {
            permutationSize *= set.size();
        }
        permutationCount += (permutationSize*hyperCube.num_copies);
        int[][] permutations = new int[permutationSize][numTerms];

        int numRepeats = permutationSize;
        // Now fill each element of input sets
        for(int i = 0 ; i < numTerms ; i++){
            numRepeats = numRepeats/sets.get(i).size();
            int numFilledEntries = 0;
            while(numFilledEntries != permutationSize){
                for(Integer elem : sets.get(i)){
                    for(int j = 0 ; j < numRepeats ; j++){
                        permutations[numFilledEntries][i] = elem;
                        numFilledEntries++;
                    }
                }
            }
        }
        return permutations;
    }

    /**
     * Create all possible permutation of a the domains of the terms
     * @param terms
     * @return
     */

    public static int[][] permute(List<Term> terms) {

        int permutaionSize = 1;
        for (Term term : terms) {
            permutaionSize *= term.domain.size();
        }

        int[][] permuations = new int[permutaionSize][terms.size()];

        for (int i = 0; i < permuations.length; i++) {
            int residue = i;
            for (int j = 0; j < terms.size(); j++) {
                int index = residue % terms.get(j).domain.size();
                residue = residue / terms.get(j).domain.size();
                permuations[i][j] = terms.get(j).domain.get(index);
            }
        }

        return permuations;

    }

    public GroundMLN handleEvidence(GroundMLN groundMln, Evidence evidence, Evidence truth, Set<String> evidence_preds, Set<String> query_preds, Set<String> hidden_preds, boolean withEM, boolean queryEvidence) throws CloneNotSupportedException {
        GroundMLN newGroundMln = new GroundMLN();
        newGroundMln.indexToGroundPredMap = groundMln.indexToGroundPredMap;
        newGroundMln.groundPredToIntegerMap = groundMln.groundPredToIntegerMap;
        Map<Integer,Integer> newGpIndexToTrueVal = new HashMap<Integer,Integer>();                  //new truth.predIdVal

        for(GroundPredicate gp : groundMln.indexToGroundPredMap.values())
            gp.groundFormulaIds.clear();
        //new groundMLN.groundPredicates

//        int formulaNum = 0;
//        if(groundHiddenPredList != null){
//            int groundHiddenPredListSize = groundHiddenPredList.size();
//            for(int i = 0; i < groundHiddenPredListSize; ++i){
//
//            }
//        }

        for(GroundFormula gf : groundMln.groundFormulas)
        {
//            formulaNum++;
//            if(formulaNum%10000 == 0)
//                System.out.println("formulaNum : "+formulaNum);
            GroundFormula newGroundFormula = new GroundFormula();
            int currentFormulaId = newGroundMln.groundFormulas.size();
            newGroundFormula.weight = gf.weight;
            newGroundFormula.formulaId = gf.formulaId;
            newGroundFormula.parentFormulaId = gf.parentFormulaId;
            newGroundFormula.originalWeight = new LogDouble(gf.originalWeight);
            boolean keepFormula = true;
            List<GroundClause> newGcList = new ArrayList<>();
            for(GroundClause gc : gf.groundClauses)
            {
                GroundClause newGc = new GroundClause();
                newGc.formulaId = currentFormulaId;
                newGc.weight = new LogDouble(gc.weight.getValue(), true);
                List<GroundPredicate> newGroundPreds = new ArrayList<>();
                boolean clauseToRemove = false;
                for(Integer gpIndex : gc.groundPredIndices)
                {
                    GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpIndex);
                    BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpIndex));

                    // If gp is not in evidence and openworld
                    // and if gp is in truth and queryEvidence
                    boolean toAdd = false;

//                    if(evidence_preds.contains(gp.symbol.symbol) && gp.symbol.world == PredicateSymbol.WorldState.open)
//                        toAdd = true;
//                    else if(query_preds.contains(gp.symbol.symbol) && (truth.predIdVal.containsKey(gpIndex) || !queryEvidence))
//                        toAdd = true;

                    if(!withEM){
                        if(query_preds.contains(gp.symbol.symbol))
                        {
                            if(truth.predIdVal.containsKey(gpIndex) || !queryEvidence)
                            {
                                if(!evidence_preds.contains(gp.symbol.symbol) || !evidence.predIdVal.containsKey(gpIndex))
                                    toAdd = true;
                            }
                        }
                    }
                    else{
                        if(hidden_preds.size() > 0){
                            if (hidden_preds.contains(gp.symbol.symbol))
                                toAdd = true;
                        }
                    }


                    if(toAdd)
                    {
                        GroundPredicate newGp = new GroundPredicate();

                        //TODO : currently, no copy c'tor called
                        newGp.symbol = new GroundPredicateSymbol(gp.symbol.id, gp.symbol.symbol,gp.symbol.values, gp.symbol.variable_types);

                        // Fill in the terms with constants
                        for(Integer term : gp.terms)
                        {
                            newGp.terms.add(term);
                        }

                        //Find index of this gp in groundPredicates list
                        MyAssert.assume(newGroundMln.groundPredToIntegerMap.containsKey(newGp));
                        int newGpIndex = newGroundMln.groundPredToIntegerMap.get(newGp);


                        if(!hidden_preds.contains(gp.symbol.symbol) && truth.predIdVal.containsKey(gpIndex))
                        {
                            int valTrue = truth.predIdVal.get(gpIndex);
                            newGpIndexToTrueVal.put(newGpIndex, valTrue);
                        }


                        newGp = newGroundMln.indexToGroundPredMap.get(newGpIndex);

                        // Check if this groundPredicate occurs first time in this ground clause. then update
                        // groundClause's data structures about this groundPredicate.
                        int newGpIndexInClause = newGc.groundPredIndices.indexOf(newGpIndex);

                        if(newGpIndexInClause == -1)
                        {
                            newGroundPreds.add(newGp);
                            newGc.groundPredIndices.add(newGpIndex);
                            newGpIndexInClause = newGc.groundPredIndices.size()-1;
                            newGc.globalToLocalPredIndex.put(newGpIndex,newGpIndexInClause);
                            newGc.grounPredBitSet.add((BitSet)b.clone());
                        }
                    }

                    else
                    {
                        if(evidence_preds.contains(gp.symbol.symbol))
                        {
                            if(evidence.predIdVal.containsKey(gpIndex)) {
                                if (b.get(evidence.predIdVal.get(gpIndex))) {
                                    clauseToRemove = true;
                                    break;
                                }
                            }
                            else
                            {
                                if(b.get(0)) // If it is closed world and not in evidence, then we assume that its true val is 0.
                                {
                                    clauseToRemove = true;
                                    break;
                                }
                            }
                        }
                        else // If it is not evidence pred, but queryevidence is true
                        {
                            if(b.get(0)) // If it is closed world and not in evidence, then we assume that its true val is 0.
                            {
                                clauseToRemove = true;
                                break;
                            }
                        }

                    }
                }
                if(clauseToRemove == false)
                {
                    if(newGc.groundPredIndices.size() > 0)
                    {
                        //clause add karna h
                        newGcList.add(newGc);
                        for (GroundPredicate gp : newGroundPreds) {
                            int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
                            newGroundFormula.groundPredIndices.add(gpIndex);
                            if (!gp.groundFormulaIds.containsKey(currentFormulaId)) {
                                gp.groundFormulaIds.put(currentFormulaId, new HashSet<Integer>());
                            }
                            gp.groundFormulaIds.get(currentFormulaId).add(newGcList.size() - 1);
                        }
                    }
                    else
                    {
                        keepFormula = false;
                        break;
                    }
                }
            }
            if(newGcList.size() > 0 && keepFormula)
            {
                newGroundFormula.groundClauses.addAll(newGcList);
                newGroundMln.groundFormulas.add(newGroundFormula);
            }
            else
            {
                // remove formulaIds from gp's formulaIdList
                for(GroundClause gc : newGcList)
                {
                    for(Integer gpId : gc.groundPredIndices)
                    {
                        GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpId);
                        gp.groundFormulaIds.remove(currentFormulaId);
                    }
                }
            }

        }
        Set<GroundPredicateSymbol> gpsSet = new HashSet<>();
        for(GroundPredicate gp : groundMln.indexToGroundPredMap.values())
        {
            GroundPredicateSymbol gps = gp.symbol;
            gpsSet.add(new GroundPredicateSymbol(gps.id, gps.symbol, gps.values, gps.variable_types));
        }
        newGroundMln.symbols.addAll(gpsSet);
        if(truth != null) {
            truth.predIdVal = newGpIndexToTrueVal;
        }
        return newGroundMln;
    }

//    public GroundMLN handleEvidenceWithEM(GroundMLN groundMln, Evidence evidence, Evidence truth, List<String> evidence_preds, List<String> query_preds, List<String> hidden_preds, boolean withEM) throws CloneNotSupportedException {
//        GroundMLN newGroundMln = new GroundMLN();
//        Map<Integer,Integer> newGpIndexToTrueVal = new HashMap<Integer,Integer>();
//        List<GroundPredicate> newGpList = new ArrayList<>();
//        Map<GroundPredicate,Integer> newGpToIntegerMap = new HashMap<GroundPredicate,Integer>();
////        int formulaNum = 0;
//        for(GroundFormula gf : groundMln.groundFormulas)
//        {
////            formulaNum++;
////            if(formulaNum%10000 == 0)
////                System.out.println("formulaNum : "+formulaNum);
//            GroundFormula newGroundFormula = new GroundFormula();
//            int currentFormulaId = newGroundMln.groundFormulas.size();
//            newGroundFormula.weight = gf.weight;
//            newGroundFormula.formulaId = gf.formulaId;
//            newGroundFormula.parentFormulaId = gf.parentFormulaId;
//            boolean keepFormula = true;
//            List<GroundClause> newGcList = new ArrayList<>();
//            for(GroundClause gc : gf.groundClauses)
//            {
//                GroundClause newGc = new GroundClause();
//                newGc.formulaId = currentFormulaId;
//                newGc.weight = new LogDouble(gc.weight.getValue(), true);
//                List<GroundPredicate> newGroundPreds = new ArrayList<>();
//                boolean clauseToRemove = false;
//                for(Integer gpIndex : gc.groundPredIndices)
//                {
//                    GroundPredicate gp = groundMln.groundPredicates.get(gpIndex);
//                    BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpIndex));
//
//                    // If gp is not in evidence and openworld
//                    // and if gp is in truth and queryEvidence
//                    boolean toAdd = false;
////                    if(!withEM)
////                    {
//////                        if(evidence_preds.contains(gp.symbol.symbol) && gp.symbol.world == PredicateSymbol.WorldState.open)
//////                            toAdd = true;
//////                        else if(query_preds.contains(gp.symbol.symbol) && (truth.predIdVal.containsKey(gpIndex) || !queryEvidence))
//////                            toAdd = true;
////
////                        if(query_preds.contains(gp.symbol.symbol))
////                        {
////                            if(truth.predIdVal.containsKey(gpIndex) || !queryEvidence)
////                            {
////                                if(!evidence_preds.contains(gp.symbol.symbol) || !evidence.predIdVal.containsKey(gpIndex)) {
////                                    toAdd = true;
////                                }
////                            }
////                        }
////                    }
////
////                    else {
//                        if (hidden_preds.contains(gp.symbol.symbol))
//                            toAdd = true;
////                    }
//
//                    if(toAdd)
//                    {
//                        GroundPredicate newGp = new GroundPredicate();
//
//                        //TODO : currently, no copy c'tor called
//                        newGp.symbol = new GroundPredicateSymbol(gp.symbol.id, gp.symbol.symbol,gp.symbol.values, gp.symbol.world, gp.symbol.variable_types);
//
//                        // Fill in the terms with constants
//                        for(Integer term : gp.terms) {
//                            newGp.terms.add(term);
//                        }
//
//                        // Check if this groundPredicate already exists, if it does not, then add it to groundPredicate List.
//                        // Note that it may happen that this clause gets removed later due to preprocessing, but still,
//                        // we need this groundPredicate, so there is no harm in adding it to groundPredicate List.
//                        int newGpIndex = -1;
//                        if(newGpToIntegerMap.containsKey(newGp))
//                            newGpIndex = newGpToIntegerMap.get(newGp);
//
//                        if(newGpIndex == -1) {
//                            newGpList.add(newGp);
//                            newGp.numPossibleValues = gp.numPossibleValues;
//                            newGpIndex = newGpList.size()-1;
//                            newGpToIntegerMap.put(newGp,newGpIndex);
//                            if(truth.predIdVal.containsKey(gpIndex))
//                            {
//                                int valTrue = truth.predIdVal.get(gpIndex);
//                                newGpIndexToTrueVal.put(newGpIndex, valTrue);
//                            }
//                        }
//
//                        newGp = newGpList.get(newGpIndex);
//
//                        // Check if this groundPredicate occurs first time in this ground clause. then update
//                        // groundClause's data structures about this groundPredicate.
//                        int newGpIndexInClause = newGc.groundPredIndices.indexOf(newGpIndex);
//
//                        if(newGpIndexInClause == -1)
//                        {
//                            newGroundPreds.add(newGp);
//                            newGc.groundPredIndices.add(newGpIndex);
//                            newGpIndexInClause = newGc.groundPredIndices.size()-1;
//                            newGc.globalToLocalPredIndex.put(newGpIndex,newGpIndexInClause);
//                            newGc.grounPredBitSet.add((BitSet)b.clone());
//                        }
//                    }
//
//                    else
//                    {
//                        if(evidence_preds.contains(gp.symbol.symbol))
//                        {
//                            if(evidence.predIdVal.containsKey(gpIndex)) {
//                                if (b.get(evidence.predIdVal.get(gpIndex))) {
//                                    clauseToRemove = true;
//                                    break;
//                                }
//                            }
//                            else
//                            {
//                                if(b.get(0)) // If it is closed world and not in evidence, then we assume that its true val is 0.
//                                {
//                                    clauseToRemove = true;
//                                    break;
//                                }
//                            }
//                        }
//                        else // If it is not evidence pred, but queryevidence is true
//                        {
//                            if(b.get(0)) // If it is closed world and not in evidence, then we assume that its true val is 0.
//                            {
//                                clauseToRemove = true;
//                                break;
//                            }
//                        }
//                    }
//                }
//                if(clauseToRemove == false)
//                {
//                    if(newGc.groundPredIndices.size() > 0)
//                    {
//                        //clause add karna h
//                        newGcList.add(newGc);
//                        for (GroundPredicate gp : newGroundPreds) {
//                            int gpIndex = newGpToIntegerMap.get(gp);
//                            newGroundFormula.groundPredIndices.add(gpIndex);
//                            if (!gp.groundFormulaIds.containsKey(currentFormulaId)) {
//                                gp.groundFormulaIds.put(currentFormulaId, new HashSet<Integer>());
//                            }
//                            gp.groundFormulaIds.get(currentFormulaId).add(newGcList.size() - 1);
//                        }
//                    }
//                    else
//                    {
//                        keepFormula = false;
//                        break;
//                    }
//                }
//            }
//            if(newGcList.size() > 0 && keepFormula)
//            {
//                newGroundFormula.groundClauses.addAll(newGcList);
//                newGroundMln.groundFormulas.add(newGroundFormula);
//            }
//            else
//            {
//                // remove formulaIds from gp's formulaIdList
//                for(GroundClause gc : newGcList)
//                {
//                    for(Integer gpId : gc.groundPredIndices)
//                    {
//                        GroundPredicate gp = newGpList.get(gpId);
//                        gp.groundFormulaIds.remove(currentFormulaId);
//                    }
//                }
//            }
//
//        }
//        newGroundMln.groundPredicates.addAll(newGpList);
//        Set<GroundPredicateSymbol> gpsSet = new HashSet<>();
//        for(GroundPredicate gp : newGpList)
//        {
//            GroundPredicateSymbol gps = gp.symbol;
//            gpsSet.add(new GroundPredicateSymbol(gps.id, gps.symbol, gps.values, gps.world, gps.variable_types));
//        }
//        newGroundMln.symbols.addAll(gpsSet);
//        truth.predIdVal = newGpIndexToTrueVal;
//        return newGroundMln;
//    }

    // Adds softevidence of type st(constant) groundclause to newGroundMLN.
    // Softevidence file contains triplets : <constant,value,weight>. For example. : <0,1,0.4> means
    // there is a soft evidence clause st(0)=1 with weight 0.4.
    // lambda is the constant with which we multiply weight of softevidence.
    // Since there is no corresponding first order formula for these ground formulas, their parentFormulaId list is empty.
    public GroundMLN addSoftEvidence(GroundMLN newGroundMln, String softEvidenceFile, double lambda, String predName) throws FileNotFoundException {
        if (softEvidenceFile == null)
            return newGroundMln;
        Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(softEvidenceFile))));
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine().replaceAll("\\s", "");
            if (line.isEmpty()) {
                continue;
            }
            String words[] = line.split(",");
            int constant = Integer.parseInt(words[0]);
            int value = Integer.parseInt(words[1]);
            double weight = Double.parseDouble(words[2]);
//            weight *= lambda;
            GroundFormula newFormula = new GroundFormula();
            int currentFormulaId = newGroundMln.groundFormulas.size();
            newFormula.formulaId = currentFormulaId;
            newFormula.weight = new LogDouble(weight * lambda, true);
            newFormula.originalWeight = new LogDouble(weight, true);
            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
            GroundClause newGroundClause = new GroundClause();
            newGroundClause.formulaId = currentFormulaId;
            newGroundClause.weight = new LogDouble(weight * lambda, true);
            for (Integer i : newGroundMln.indexToGroundPredMap.keySet()) {
                GroundPredicate gp = newGroundMln.indexToGroundPredMap.get(i);
                if (gp.symbol.symbol.equals(predName)) {
                    if (gp.terms.get(0).equals(constant)) {
                        BitSet b = new BitSet(gp.numPossibleValues);
                        b.set(value);
                        newGroundClause.globalToLocalPredIndex.put(i, 0);
                        newGroundClause.grounPredBitSet.add(b);
                        newGroundClause.groundPredIndices.add(i);
                        newGroundClauseList.add(newGroundClause);
                        newFormula.groundPredIndices.add(i);
                        if(!gp.groundFormulaIds.containsKey(currentFormulaId))
                        {
                            gp.groundFormulaIds.put(currentFormulaId, new HashSet<Integer>());
                        }
                        gp.groundFormulaIds.get(currentFormulaId).add(newGroundClauseList.size()-1);
                        break;
                    }
                }
            }
            newFormula.groundClauses.add(newGroundClause);
            newGroundMln.groundFormulas.add(newFormula);
        }
        return newGroundMln;
    }
    // For each constant, size of feature vector is #firstorderformulas * #positions in each first order formula at which
    // that constant can come
//    public Map<Integer,List<Integer>> getFeatureVectors(GroundMLN groundMln, int numFormulas, Evidence truth, String typeName, Set<Integer> constants, boolean closedWorld) {
//        Map<Integer, List<Integer>> featureVectors = new HashMap<>();
//        for(int constant : constants)
//        {
//            featureVectors.put(constant, new ArrayList<Integer>(Collections.nCopies(numFormulas,0)));
//        }
//        for(GroundFormula gf : groundMln.groundFormulas)
//        {
//            boolean isFormulaSatisfied = true;
//            for(GroundClause gc : gf.groundClauses)
//            {
//                boolean isClauseSatisfied = false;
//                for(Integer gpIndex : gc.groundPredIndices)
//                {
//                    int truthVal = 0;
//                    if(truth.predIdVal.containsKey(gpIndex))
//                    {
//                        truthVal = truth.predIdVal.get(gpIndex);
//                    }
//                    else
//                    {
//                        if(!closedWorld)
//                            continue;
//                    }
//                    int localIndex = gc.globalToLocalPredIndex.get(gpIndex);
//                    isClauseSatisfied = gc.grounPredBitSet.get(localIndex).get(truthVal);
//                    if(isClauseSatisfied)
//                        break;
//                }
//                if(!isClauseSatisfied)
//                {
//                    isFormulaSatisfied = false;
//                    break;
//                }
//            }
//            // If formula is satisfied then update feature vector of each constant appearing in this formula
//            if(isFormulaSatisfied)
//            {
//                for (int i = 0; i < gf.parentFormulaId.size(); i++) {
//                    int pfId = gf.parentFormulaId.get(i);
//                    for (GroundClause gc : gf.groundClauses) {
//                        for (Integer gpIndex : gc.groundPredIndices) {
//                            GroundPredicate gp = groundMln.groundPredicates.get(gpIndex);
//                            for (int j = 0; j < gp.terms.size(); j++) {
//                                if(gp.symbol.variable_types.get(j).equals(typeName))
//                                {
//                                    int constant = gp.terms.get(j);
//                                    int numSatisfied = featureVectors.get(constant).get(pfId);
//                                    featureVectors.get(constant).set(pfId, numSatisfied+gf.numCopies.get(i));
//                                }
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//
//        return featureVectors;
//    }

    /**
     * This function creates a list of {@link GroundPredicate}s from a given MLN
     * @param mln Input MLN
     * @return List of ground predicates
     */
    public List<GroundPredicate> createGroundPredicates(MLN mln, Map<GroundPredicate, Integer> groundPredToIntegerMap, Set<String> queryPreds) {
        List<GroundPredicate> groundPredsList = new ArrayList<GroundPredicate>();
        // For each pred symbol in MLN, create its groundings based on domain of its terms

        for(PredicateSymbol ps : mln.symbols)
        {
            if(!queryPreds.contains(ps.symbol))
                continue;
            // Create set of domains for terms of this ps
            List<Set<Integer>> termDomains = new ArrayList<>();
            for(String varType : ps.variable_types)
            {
                termDomains.add(mln.varTypeToDomainMap.get(varType));
            }
            List<List<Integer>> permutations = OtherUtils.cartesianProd(termDomains);
            for(List<Integer> perm : permutations)
            {
                GroundPredicate gp = new GroundPredicate();
                gp.symbol = new GroundPredicateSymbol(ps.id,ps.symbol,ps.values,ps.variable_types);
                gp.numPossibleValues = gp.symbol.values.values.size();
                gp.terms.addAll(perm);
                groundPredToIntegerMap.put(gp, groundPredsList.size());
                groundPredsList.add(gp);
            }

        }
        return groundPredsList;
    }

    public List<GroundPredicate> createGroundPredicatesFromPredsHashMap(HashMap<PredicateSymbol, ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap, Map<GroundPredicate, Integer> groundPredToIntegerMap) {
        List<GroundPredicate> groundPredsList = new ArrayList<GroundPredicate>();
        for(PredicateSymbol ps : predsHyperCubeHashMap.keySet()){
            int lastIndex = predsHyperCubeHashMap.get(ps).size()-1;
            List<HyperCube> unknownHyperCubes = predsHyperCubeHashMap.get(ps).get(lastIndex);
            if(unknownHyperCubes.size() != 0)
            {
                for(HyperCube hc : unknownHyperCubes)
                {
                    int [][] permutations = permuteHyperCube(hc);
                    for (int i = 0; i < permutations.length; i++) {
                        GroundPredicate gp = new GroundPredicate();
                        gp.symbol = new GroundPredicateSymbol(ps.id,ps.symbol,ps.values,ps.variable_types);
                        gp.numPossibleValues = gp.symbol.values.values.size();
                        for (int j = 0; j < permutations[i].length; j++) {
                            gp.terms.add(permutations[i][j]);
                        }
                        groundPredToIntegerMap.put(gp, groundPredsList.size());
                        groundPredsList.add(gp);
                    }
                }
            }
        }
        return groundPredsList;
    }

    public void removeEvidenceGroundPreds(GroundMLN groundMln, Evidence evidence) {
        for(Integer gpId : evidence.predIdVal.keySet())
        {
            GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpId);
            MyAssert.assume(gp.groundFormulaIds.size() == 0);
            groundMln.indexToGroundPredMap.remove(gpId);
            groundMln.groundPredToIntegerMap.remove(gp);
        }
    }

}