package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;

import java.util.*;

/**
 * Created by Happy on 2/23/17.
 */
public class State {
    public MLN mln;
    public GroundMLN groundMLN;
    public Map<Integer, Integer> truthVals = new HashMap<>(); // For each groundPredicate in mln.groundPredicates, stores its truthval
    public List<Set<Integer>> falseClausesSet = new ArrayList<>(); // for each groundformula, stores set of groundClauseIds which are false in this state
    public List<List<Integer>> numTrueLiterals = new ArrayList<>(); // for each groundformula, for each clauseId, stores numSatLiterals in that clause
    public Map<Integer, List<Double>> wtsPerPredPerVal = new HashMap<>(); // For each GroundPred, stores sat wts for each value
    public ArrayList<Integer> groundedGfIndicesList = new ArrayList<>(); // Contains indices of those gfs in groundMLN, which were not present in first order mln, but directly added during code. Thses gfs have parentFormulaId as -1.

    public State(MLN mln, GroundMLN groundMLN) {
        this.mln = mln;
        this.groundMLN  = groundMLN;
        //int numGroundPreds = groundMLN.groundPredicates.size();
        for(GroundPredicate gp : groundMLN.groundPredToIntegerMap.keySet())
        {
            int gpId = groundMLN.groundPredToIntegerMap.get(gp);
            truthVals.put(gpId, 0);
            wtsPerPredPerVal.put(gpId, new ArrayList<>(Collections.nCopies(gp.numPossibleValues,0.0)));
        }
        int numGroundFormulas = groundMLN.groundFormulas.size();
        for(int i = 0 ; i < numGroundFormulas ; i++)
        {
            falseClausesSet.add(new HashSet<Integer>());
            int numGroundClauses = groundMLN.groundFormulas.get(i).groundClauses.size();
            numTrueLiterals.add(new ArrayList<Integer>());
            for(int j = 0 ; j < numGroundClauses ; j++)
            {
                numTrueLiterals.get(i).add(0);
            }
            if(groundMLN.groundFormulas.get(i).parentFormulaId.isEmpty())
                groundedGfIndicesList.add(i);
        }
    }

    public State(MLN mln, GroundMLN groundMLN, Evidence truthVals)
    {
        this.mln = mln;
        this.groundMLN  = groundMLN;
        //int numGroundPreds = groundMLN.groundPredicates.size();
        for(GroundPredicate gp : groundMLN.groundPredToIntegerMap.keySet())
        {
            int gpId = groundMLN.groundPredToIntegerMap.get(gp);
            if(!truthVals.predIdVal.containsKey(gpId)) {
                this.truthVals.put(gpId, 0);
            }
            else
            {
                this.truthVals.put(gpId, truthVals.predIdVal.get(gpId));
            }

            wtsPerPredPerVal.put(gpId, new ArrayList<>(Collections.nCopies(gp.numPossibleValues,0.0)));
        }
        int numGroundFormulas = groundMLN.groundFormulas.size();
        for(int i = 0 ; i < numGroundFormulas ; i++)
        {
            falseClausesSet.add(new HashSet<Integer>());
            int numGroundClauses = groundMLN.groundFormulas.get(i).groundClauses.size();
            numTrueLiterals.add(new ArrayList<Integer>());
            for(int j = 0 ; j < numGroundClauses ; j++)
            {
                numTrueLiterals.get(i).add(0);
            }
            if(groundMLN.groundFormulas.get(i).parentFormulaId.isEmpty())
                groundedGfIndicesList.add(i);
        }
    }

//    public void setGroundFormulaWtsToParentWts(MLN mln) {
//        for(GroundFormula gf : groundMLN.groundFormulas)
//        {
//            int parentFormulaId = gf.parentFormulaId;
//            gf.weight = new LogDouble(mln.formulas.get(parentFormulaId).weight.getValue(), true);
//        }
//    }
//
//    public void setGroundFormulaWtsToParentWtsSoftEvidence(MLN mln, double lambda) {
//        for(GroundFormula gf : groundMLN.groundFormulas)
//        {
//            int parentFormulaId = gf.parentFormulaId;
//            if(parentFormulaId != -1)
//                gf.weight = new LogDouble(mln.formulas.get(parentFormulaId).weight.getValue(), true);
//            else
//                gf.weight = new LogDouble(lambda*gf.originalWeight.getValue(),true);
//        }
//    }

    public double[] getNumTrueGndings(int numWts) {
        double []numTrueGndings = new double[numWts];
        for(GroundFormula gf : groundMLN.groundFormulas)
        {
            boolean isFormulaSatisfied = true;
            for(GroundClause gc : gf.groundClauses)
            {
                boolean isClauseSatisfied = false;
                for(Integer gpId : gc.groundPredIndices)
                {
                    BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                    int trueVal = truthVals.get(gpId);
                    isClauseSatisfied |= b.get(trueVal);
                    if(isClauseSatisfied)
                        break;
                }
                isFormulaSatisfied &= isClauseSatisfied;
                if(!isFormulaSatisfied)
                    break;
            }
            if(isFormulaSatisfied)
            {
                List<Integer> parentFormulaId = gf.parentFormulaId;
                List<Integer> numCopies = gf.numCopies;
                if(!parentFormulaId.isEmpty()) {
                    for (int i = 0; i < parentFormulaId.size(); i++) {
                        numTrueGndings[parentFormulaId.get(i)] += numCopies.get(i);
                    }
                }
                else
                {
                    numTrueGndings[numWts-1] += gf.originalWeight.getValue();
                }
            }
        }
        return numTrueGndings;
    }

    public void setTruthVals(Evidence truthVals) {
        for(int gpId : truthVals.predIdVal.keySet())
        {
            this.truthVals.put(gpId, truthVals.predIdVal.get(gpId));
        }
    }
}