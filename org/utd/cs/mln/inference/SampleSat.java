package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.*;

import java.util.*;

import static org.utd.cs.mln.alchemy.util.OtherUtils.getUniformAssignment;

/**
 * This class implements SampleSat algorithm, which samples a satisfying assignment uninformly (given a list of
 * clauses)
 * The algo works as follows :
 * Start with a random state
 * while(not all clauses satisfied):
 *  with prob p, make walksat move, otherwise make Simulated Annealing move (SA Move) (we will keep p=0.8)
 * return state
 * Walksat move :
 *  Pick an unsatisfied clause randomly, and flip one of its literal (which makes max number of clauses satisfiable)
 * SA move :
 *  Pick a literal (ground pred) randomly. Find its delta, where delta is difference of satisfied clauses if we flip
 *  this literal and satisfied clauses if don't flip this literal.
 *  If delta > 0:
 *      Flip the literal with prob 1.0
 *  else:
 *      Flip the literal with prob exp(delta/temp) (we set temp=14)
 *
 * Created by Happy on 5/28/18.
 */
public class SampleSat {


    private final State state;
    private final Set<Integer> groundFormulaIds;
    private Map<Integer, Integer> truthValues;
    private Map<Integer, Integer> numSatLiterals; // Key : groundFormulaId, value : numSatLiterals
    private int numSatFormulas;
    private List<Integer> unsatFormulaIds;
    private final double p = 0.8; // prob by which we make walksat move

    /***
     * @param state_ State in which all infor about groundFormulas are present
     * @param groundFormulaIds_ Indices of groundformulas which are required to be satisfied
     * @param truthValues_ Map in which generated sample will be stored. key : gpId, value : truthvalue
     */
    public SampleSat(State state_, Set<Integer> groundFormulaIds_, Map<Integer, Integer> truthValues_) {
        state = state_;
        groundFormulaIds = groundFormulaIds_;
        truthValues = truthValues_;
        numSatLiterals = new HashMap<>();
        unsatFormulaIds = new ArrayList<>();
        for(Integer gfId : groundFormulaIds)
        {
            numSatLiterals.put(gfId, 0);
        }
    }

    /**
     * This function performs one samplesatstep.
     */
    public void performSampleSatStep() {
        randomInitializeTruthVals();
        initializeNumSatValues();
        int numFormulas = groundFormulaIds.size();
        while(numSatFormulas != numFormulas)
        {
            Random rand = new Random();
            if(rand.nextDouble() < p)
            {
                // Make a walksat move
                performWalkSatStep();
            }
            else
            {
                // Make SA move
            }
        }
    }

    private void performWalkSatStep() {
        // Pick an unsatisfied formula randomly
        int formulaId = unsatFormulaIds.get((int)(Math.random()*unsatFormulaIds.size()));
        GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
        // If the weight is negative, then randomly pick any unsatisfied literal and flip it
        if(gf.weight.getValue() < 0)
        {
            // TODO : This is potentially slow operation
            Map<Integer, BitSet> unsatLitsMap = new HashMap<>();
            List<Integer> unsatLits = new ArrayList<>();
            GroundClause gc = gf.groundClauses.get(0);
            int numPreds = gc.groundPredIndices.size();
            for(int k = 0 ; k < numPreds ; k++) {
                int globalPredIndex = gc.groundPredIndices.get(k);
                int currentAssignment = truthValues.get(globalPredIndex);
                if (!gc.grounPredBitSet.get(k).get(currentAssignment)) {
                    unsatLitsMap.put(globalPredIndex, gc.grounPredBitSet.get(k));
                    unsatLits.add(globalPredIndex);
                }
            }
            // unsatLits contains all unsat Lits, now randomly pick one of the lits and flip it.
            int predId = unsatLits.get((int)(Math.random()*unsatLits.size()));
            GroundPredicate gp = state.groundMLN.indexToGroundPredMap.get(predId);
            // Now find the value to which this predId should be flipped to. The value which satisfied most number
            // of formulas should be the answer
            int numPossibleVals = gp.numPossibleValues;
            BitSet b = unsatLitsMap.get(predId);
            int maxFormulas = 0;
            List<Integer> numSatFormulasPerVal = getNumSatFormulasPerVal(predId);
        }
    }

    private List<Integer> getNumSatFormulasPerVal(int gpId) {
        GroundPredicate gp = state.groundMLN.indexToGroundPredMap.get(gpId);
        // Loop over all formulas in which this gpId occurs
        Map<Integer, Set<Integer>> formulaIds = gp.groundFormulaIds;
        for(Integer formulaId : formulaIds.keySet())
        {
            GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
            // HACK : assume clausal form
            for(Integer clauseId : formulaIds.get(formulaId))
            {
                BitSet clauseBitSet = new BitSet(gp.numPossibleValues);
            }
        }
        return new ArrayList<>();
    }

    /**
     * Randomly initialize truth values of all ground predicates
     */
    void randomInitializeTruthVals() {
        for(GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet())
        {
            int numPossibleVals = gp.numPossibleValues;
            int assignment = getUniformAssignment(numPossibleVals);
            int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
            truthValues.put(gpId,assignment);
        }
    }

    // According to present truthVals, set numSatLiterals
    // HACK : assumes clausal form i.e. each formula contains only one clause
    void initializeNumSatValues() {
        for(Integer formulaId : groundFormulaIds) {
            GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
            GroundClause gc = gf.groundClauses.get(0);
            int numPreds = gc.groundPredIndices.size();
            int numSatLiteralsLocal = 0;
            for(int k = 0 ; k < numPreds ; k++) {
                int globalPredIndex = gc.groundPredIndices.get(k);
                int currentAssignment = truthValues.get(globalPredIndex);
                if (gc.grounPredBitSet.get(k).get(currentAssignment)) {
                    numSatLiteralsLocal++;
                }
            }
            numSatLiterals.put(formulaId, numSatLiteralsLocal);

            if(gf.weight.getValue() >= 0)
            {
                if (numSatLiteralsLocal > 0)
                    numSatFormulas++;
                if(numSatLiteralsLocal == 0)
                    unsatFormulaIds.add(formulaId);
            }
            else
            {
                if(numSatLiteralsLocal == 0)
                    numSatFormulas++;
                if (numSatLiteralsLocal > 0)
                    unsatFormulaIds.add(formulaId);
            }

        }// end of formula
    }

}
