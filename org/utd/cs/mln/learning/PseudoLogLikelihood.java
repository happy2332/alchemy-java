//package org.utd.cs.mln.learning;
//
//import org.utd.cs.mln.alchemy.core.*;
//import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
//import org.utd.cs.mln.alchemy.util.OtherUtils;
//
//import java.util.*;
//
///**
// * Created by Happy on 4/4/18.
// */
//public class PseudoLogLikelihood {
//    private ArrayList<Map<Integer, List<Integer>>> satCounts;
//    State state_;
//    List<List<Double>> wtsPerPredPerVal_;
//    Map<String,Integer> predToNumGndingsMap;
//    public PseudoLogLikelihood(State state) {
//        state_ = state;
//        wtsPerPredPerVal_ = new ArrayList<>();
//        for(GroundPredicate gp : state_.groundMLN.groundPredToIntegerMap.keySet()) {
//            int numVals = gp.numPossibleValues;
//            wtsPerPredPerVal_.add(new ArrayList<Double>(Collections.nCopies(numVals,0.0)));
//        }
//        satCounts = new ArrayList<>();
//        predToNumGndingsMap = new HashMap<>();
//    }
//
//    public double getPseudoLogLikelihood() {
//        countSatVals();
//        createPredToNumGndingsMap();
//        setWtsPerPredPerVal();
//        double pll = 0;
//        for(Integer gpId : state_.groundMLN.indexToGroundPredMap.keySet()) {
//            double d = OtherUtils.getProbabilityDistribution(wtsPerPredPerVal_.get(gpId)).get(state_.truthVals.get(gpId));
//            String predName = state_.groundMLN.indexToGroundPredMap.get(gpId).symbol.symbol;
//            pll += (Math.log(d)/predToNumGndingsMap.get(predName));
//        }
//
//        return pll;
//    }
//
//    private void setWtsPerPredPerVal() {
//        MLN mln = state_.mln;
//        for (Integer gpId : state_.groundMLN.indexToGroundPredMap.keySet()) {
//            Map<Integer, List<Integer>> satCountsPerPred = satCounts.get(gpId);
//            int numPossibleVals = state_.groundMLN.indexToGroundPredMap.get(gpId).numPossibleValues;
//            for(int formulaId : satCountsPerPred.keySet())
//            {
//                double formulaWt = mln.formulas.get(formulaId).weight.getValue();
//                List<Integer> satCountsPerPredPerFormula = satCountsPerPred.get(formulaId);
//                for (int val = 0; val < numPossibleVals; val++) {
//                    int satCount = satCountsPerPredPerFormula.get(val);
//                    double oldVal = wtsPerPredPerVal_.get(gpId).get(val);
//                    wtsPerPredPerVal_.get(gpId).set(val, oldVal + satCount*formulaWt);
//                }
//            }
////            if(priorSoftEvidence)
////            {
////                if(softEvidencePerPredPerVal.get(domainIndex).containsKey(gpId))
////                {
////                    double lambda = weights[numWts-1];
////                    for (int val = 0; val < numPossibleVals; val++) {
////                        double oldVal = state.wtsPerPredPerVal.get(gpId).get(val);
////                        state.wtsPerPredPerVal.get(gpId).set(val, oldVal + lambda * softEvidencePerPredPerVal.get(domainIndex).get(gpId).get(val));
////                    }
////                }
////            }
//
//        }
//    }
//
//    void countSatVals() {
//        for (int i = 0; i < state_.groundMLN.groundPredicates.size(); i++) {
//            satCounts.add(new HashMap<Integer, List<Integer>>());
//        }
//        MLN mln = state_.mln;
//        // for each formula, compute sat counts
//        for (int formulaId = 0; formulaId < mln.formulas.size(); formulaId++) {
//            Formula formula = mln.formulas.get(formulaId);
//            Set<Term> formulaWiseTermToGround = new HashSet<Term>();
//            for (WClause clause : formula.clauses) {
//                for (Atom atom : clause.atoms) {
//                    for (int j = 0; j < atom.terms.size(); j++) {
//                        Term term = atom.terms.get(j);
//                        formulaWiseTermToGround.add(term);
//                    }
//                }
//            }
//            countSatValsPerFormula(formulaId, satCounts, new ArrayList<Term>(formulaWiseTermToGround));
//        }
//    }
//
//    private void countSatValsPerFormula(int formulaId, List<Map<Integer, List<Integer>>> satCounts, ArrayList<Term> terms) {
//        MLN mln = state_.mln;
//        int[][] permutations = FullyGrindingMill.permute(terms);
//        List<GroundPredicate> groundPredicates = state_.groundMLN.groundPredicates;
//        Map<GroundPredicate, Integer> groundPredToIntegerMap = state_.groundMLN.groundPredToIntegerMap;
//        Formula formula = mln.formulas.get(formulaId);
//        for (int i = 0; i < permutations.length; i++) {
//            GroundFormula newFormula = new GroundFormula();
//            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
//            Set<Integer> formulaGpIndices = new HashSet<>();
//            for (int c = 0; c < formula.clauses.size(); c++) {
//                WClause clause = formula.clauses.get(c);
//                GroundClause newGroundClause = new GroundClause();
//                Map<Integer, BitSet> gpIndexToSatVals = new HashMap<>();
//                // Iterate over each first order atom, and create ground atom for it.
//                for (int j = 0; j < clause.atoms.size(); j++) {
//                    boolean sign = clause.sign.get(j);
//                    Atom oldAtom = clause.atoms.get(j); // first order atom
//                    int valTrue = clause.valTrue.get(j);
//                    GroundPredicate gp = new GroundPredicate(); // GroundPredicate to create
//                    gp.symbol = new GroundPredicateSymbol(oldAtom.symbol.id, oldAtom.symbol.symbol, oldAtom.symbol.values, oldAtom.symbol.variable_types);
//                    // Fill in the terms with constants
//                    for (Term term : oldAtom.terms) {
//                        int termIndex = terms.indexOf(term);
//                        gp.terms.add(permutations[i][termIndex]);
//                    }
//
//                    int numPossibleValues = oldAtom.symbol.values.values.size();
//                    gp.numPossibleValues = numPossibleValues;
//                    int gpIndex = groundPredToIntegerMap.get(gp);
//                    formulaGpIndices.add(gpIndex);
//
//                    // Check if this groundPredicate occurs first time in this ground clause. then update
//                    // groundClause's data structures about this groundPredicate.
//                    int gpIndexInClause = newGroundClause.groundPredIndices.indexOf(gpIndex);
//                    //GroundAtom newGroundAtom = new GroundAtom(gpIndex, gpIndexInClause, valTrue, sign);
//                    if (gpIndexInClause == -1) {
//                        newGroundClause.groundPredIndices.add(gpIndex);
//                        gpIndexInClause = newGroundClause.groundPredIndices.size() - 1;
//                        newGroundClause.globalToLocalPredIndex.put(gpIndex, gpIndexInClause);
//                        gpIndexToSatVals.put(gpIndexInClause, new BitSet(gp.numPossibleValues));
//                    }
//
//                    // Now once we have added new ground Atom, we need to check if ground clause gets satisfied or not.
//                    BitSet gpBitSet = new BitSet(gp.numPossibleValues);
//                    gpBitSet.set(valTrue);
//                    if (sign == true)
//                        gpBitSet.flip(0, gp.numPossibleValues);
//                    gpBitSet.or(gpIndexToSatVals.get(gpIndexInClause));
//
//                    // If all bits are set for this groundPred, then this clause will always be satisfied and hence,
//                    // shouldn't be added into groundformula. Note that, although at this point, we know that
//                    // this clause shouldn't be added, but still we shouldn't just break out of this loop, as we
//                    // need to add groundPredicates, but we shouldn't add any clauseInfo into groundPredicates appearing
//                    // in this clause.
//
//                    gpIndexToSatVals.put(gpIndexInClause, gpBitSet);
//                }
//                for (int gpId = 0; gpId < newGroundClause.groundPredIndices.size(); gpId++) {
//                    BitSet b = gpIndexToSatVals.get(gpId);
//                    newGroundClause.grounPredBitSet.add(b);
//                }
//                newGroundClauseList.add(newGroundClause);
//            }
//            newFormula.groundClauses.addAll(newGroundClauseList);
//
//            for (Integer gpIndex : formulaGpIndices) {
//                int numPossibleValues = groundPredicates.get(gpIndex).numPossibleValues;
//                int tempTrueValue = state_.truthVals.get(gpIndex);
//                int ftcPerValue[] = new int[numPossibleValues];
//
//                for (int predValue = 0; predValue < numPossibleValues; ++predValue) {
//                    state_.truthVals.set(gpIndex, predValue);
//                    boolean isFormulaSatisfied = true;
//                    for (GroundClause gc : newFormula.groundClauses) {
//                        boolean isClauseSatisfied = false;
//                        for (int gpId : gc.groundPredIndices) {
//                            BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
//                            isClauseSatisfied = b.get(state_.truthVals.get(gpId));
//                            if (isClauseSatisfied) {
//                                break;
//                            }
//                        }
//                        isFormulaSatisfied = isClauseSatisfied;
//                        if (!isFormulaSatisfied) {
//                            break;
//                        }
//                    }
//                    if (isFormulaSatisfied) {
//
//                        ++ftcPerValue[predValue];
//
//                    }
//                }
//                Map<Integer, List<Integer>> satCountsMap = satCounts.get(gpIndex);
//                if (!satCountsMap.containsKey(formulaId)) {
//                    satCountsMap.put(formulaId, new ArrayList<>(Collections.nCopies(numPossibleValues, 0)));
//                }
//                for (int j = 0; j < numPossibleValues; j++) {
//                    satCountsMap.get(formulaId).set(j, satCountsMap.get(formulaId).get(j) + ftcPerValue[j]);
//                }
//                state_.truthVals.set(gpIndex, tempTrueValue);
//            }
//        }
//    }
//
//    private void createPredToNumGndingsMap() {
//        GroundMLN groundMln = state_.groundMLN;
//        List<GroundPredicate> groundPreds = groundMln.groundPredicates;
//        for(GroundPredicate gp : groundPreds)
//        {
//            String predName = gp.symbol.symbol;
//            if(!predToNumGndingsMap.containsKey(predName))
//            {
//                predToNumGndingsMap.put(predName,0);
//            }
//            int oldVal = predToNumGndingsMap.get(predName);
//            predToNumGndingsMap.put(predName,oldVal+1);
//        }
//
//    }
//
//
//}
