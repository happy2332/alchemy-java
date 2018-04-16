package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.OtherUtils;

import java.io.PrintWriter;
import java.util.*;

/**
 * Abstract class from which all MCMC algorithms are sub classed.
 * @author Happy
 * @since 03/28/18
 * @see MCMCParams,GibbsSampler_v3
 */
public abstract class MCMC extends Inference {
    // No. of chains which MCMC will use
    int numChains;
    // Min. no. of burn-in steps MCMC will take per chain
    int burnMinSteps;
    // Max. no. of burn-in steps MCMC will take per chain
    int burnMaxSteps;
    // Min. no. of sampling steps MCMC will take per chain
    int minSteps;
    // Max. no. of sampling steps MCMC will take per chain
    int maxSteps;
    // Max. no. of seconds MCMC should run
    int maxSeconds;

    boolean mcmcdebug = true;

    // Truth values in each chain for each ground predicate (truthValues[c][p])
    Map<Integer, Integer> truthValues [];

    // Num. of satisfying literals in each chain for each ground formula for each ground clause
    // numSatLiterals[chain][formula][clause]
    List<List<List<Integer>>> numSatLiterals;

    // For each chain, For each formula, set of clauseIds which are false in this truthVals configuration
    // falseClausesSet[chain][formula][
    List<List<Set<Integer>>> falseClausesSet;

    // For each chain, Holds number of times a ground predicate is set to a value
    // numValPerPred[c][g][val]
    List<Map<Integer, List<Double>>> numValPerChainPerPred_;

    // For each chain, for each ground predicate, for each possible value, stores satweight
    // wtsPerPredVal[c][g][val]
    List<Map<Integer, List<Double>>> wtsPerPredPerVal;

    /**
     * Constructor: User-set parameters are set.
     * @see Inference#Inference(State, int, boolean, boolean)
     */
    public MCMC(State state, int seed, boolean trackFormulaTrueCnts, boolean priorSoftEvidence, MCMCParams params) {
        super(state, seed, trackFormulaTrueCnts, priorSoftEvidence);
        this.numChains = params.numChains;
        this.burnMinSteps = params.burnMinSteps;
        this.burnMaxSteps = params.burnMaxSteps;
        this.minSteps = params.minSteps;
        this.maxSteps = params.maxSteps;
        this.maxSeconds = params.maxSeconds;
    }

    /**
     * Allocates memory to all truthValues, numSatLiterals, falseClauseSet and wtsPerPredPerVal
     */
    void initTruthValsAndSatTrueLits() {
        truthValues = new Map[numChains];
        numSatLiterals = new ArrayList<>();
        wtsPerPredPerVal = new ArrayList<>();
        falseClausesSet = new ArrayList<>();
        int numFormulas = state.groundMLN.groundFormulas.size();
        for (int c = 0; c < numChains; c++) {
            numSatLiterals.add(new ArrayList<List<Integer>>());
            wtsPerPredPerVal.add(new HashMap<Integer, List<Double>>());
            falseClausesSet.add(new ArrayList<Set<Integer>>());
            for (int formulaId = 0; formulaId < numFormulas; formulaId++) {
                int numClauses = state.groundMLN.groundFormulas.get(formulaId).groundClauses.size();
                numSatLiterals.get(c).add(new ArrayList<>(Collections.nCopies(numClauses,0)));
                falseClausesSet.get(c).add(new HashSet<Integer>());
            }
            truthValues[c] = new HashMap<>();
            for (GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet()) {
                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
                truthValues[c].put(gpId, 0);
                wtsPerPredPerVal.get(c).put(gpId, new ArrayList<Double>(Collections.nCopies(gp.numPossibleValues,0.0)));
            }
        }
    }

    /**
     * Randomly initialize truth values of all ground predicates for each chain
     */
    void randomInitializeTruthVals() {
        for (int c = 0; c < numChains; c++) {
            for(GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet())
            {
                int numPossibleVals = gp.numPossibleValues;
                int assignment = getUniformAssignment(numPossibleVals);
                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
                truthValues[c].put(gpId,assignment);
            }
        }
    }

    /**
     * Generates a random integer from 0 to n-1
     */
    int getUniformAssignment(int n) {
        int r = rand.nextInt(n);
        return r;
    }

    // According to present truthVals, set falseClausesSet and numSatLiterals
    void initializeNumSatValues() {
        GroundMLN gm = state.groundMLN;
        int numFormulas = state.groundMLN.groundFormulas.size();
        for (int c = 0; c < numChains; c++) {
            for(int formulaId = 0 ; formulaId < numFormulas ; formulaId++) {
                GroundFormula gf = gm.groundFormulas.get(formulaId);
                Set<Integer> falseClauseIds = new HashSet<>();
                int numClauses = gf.groundClauses.size();

                for(int clauseId = 0 ; clauseId < numClauses ; clauseId++) {
                    GroundClause gc = gf.groundClauses.get(clauseId);
                    int numPreds = gc.groundPredIndices.size();
                    int numSatLiteralsLocal = 0;
                    for(int k = 0 ; k < numPreds ; k++) {
                        int globalPredIndex = gc.groundPredIndices.get(k);
                        int currentAssignment = truthValues[c].get(globalPredIndex);
                        if (gc.grounPredBitSet.get(k).get(currentAssignment)) {
                            numSatLiteralsLocal++;
                        }
                    }
                    numSatLiterals.get(c).get(formulaId).set(clauseId, numSatLiteralsLocal);

                    if (numSatLiteralsLocal == 0) {
                        falseClauseIds.add(clauseId);
                    }
                } // end of clause
                falseClausesSet.get(c).set(formulaId, falseClauseIds);
            }// end of formula
        }
    }

    /**
     * Initializes structure for holding number of times a predicate was set
     * to each value.
     */
    void initNumValPerChainPerPred()
    {
        numValPerChainPerPred_ = new ArrayList<>();
        for (int chainIdx = 0; chainIdx < numChains; chainIdx++) {
            numValPerChainPerPred_.add(new HashMap<Integer, List<Double>>());
            for (GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet()) {
                int numPossibleVals = gp.numPossibleValues;
                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
                numValPerChainPerPred_.get(chainIdx).put(gpId, new ArrayList<>(Collections.nCopies(numPossibleVals,0.0)));
            }
        }

    }

    public void updateWtsForGndPreds(List<Integer> affectedGndPredIndices, int chainIndex) {
        Map<Integer, GroundPredicate> indexToGroundPredMap = state.groundMLN.indexToGroundPredMap;
        for (Integer i : affectedGndPredIndices) {
            GroundPredicate gp = indexToGroundPredMap.get(i);
            int numPossibleVals = gp.numPossibleValues;
            double wtsPerVal[] = new double[numPossibleVals];
            Map<Integer, Set<Integer>> formulaIds = gp.groundFormulaIds;

            for (Integer formulaId : formulaIds.keySet()) {
                GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
                double wt = gf.weight.getValue();
                Set<Integer> tempSet = new HashSet<Integer>();
                tempSet.addAll(falseClausesSet.get(chainIndex).get(formulaId));
                tempSet.removeAll(formulaIds.get(formulaId));
                BitSet formulaBitSet = new BitSet(numPossibleVals);
                formulaBitSet.flip(0, numPossibleVals);

                // If there is a clause which is false, and doesn't contain gp, then this formula is always false
                if (tempSet.size() == 0) {
                    for (Integer cid : formulaIds.get(formulaId)) {
                        BitSet clauseBitSet = new BitSet(numPossibleVals);
                        GroundClause gc = gf.groundClauses.get(cid);
                        int localPredIndex = gc.globalToLocalPredIndex.get(i);
                        int numSatLiteralsLocal = numSatLiterals.get(chainIndex).get(formulaId).get(cid);
                        if (numSatLiteralsLocal > 1)
                            clauseBitSet.flip(0, numPossibleVals);
                        else if (numSatLiteralsLocal == 1) {
                            BitSet b = gc.grounPredBitSet.get(localPredIndex);
                            if (b.get(truthValues[chainIndex].get(i))) {
                                clauseBitSet = (BitSet) b.clone();
                            } else {
                                clauseBitSet.flip(0, numPossibleVals);
                            }

                        } else {
                            BitSet b = gc.grounPredBitSet.get(localPredIndex);
                            clauseBitSet = (BitSet) b.clone();
                        }
                        formulaBitSet.and(clauseBitSet);
                    }// end clauses loop

                    int startIndex = 0;
                    while (startIndex < numPossibleVals) {
                        int index = formulaBitSet.nextSetBit(startIndex);
                        if (index == -1)
                            break;
                        wtsPerVal[index] += wt;
                        startIndex = index + 1;
                    }
                }// end if condition
            } //end formulas loops

            List<Double> tempWts = new ArrayList<>();
            for (double d : wtsPerVal) {
                tempWts.add(d);
            }
            wtsPerPredPerVal.get(chainIndex).put(i, tempWts);
        }
    }

    public int get_probabilistic_assignment(List<Double> satWeights) {
        List<Double> probabilities = OtherUtils.getProbabilityDistribution(satWeights);
        return getRandomAssignment(probabilities);
    }

    private int getRandomAssignment(List<Double> probabilities) {
//        double p = Math.random();
//        Random rand = new Random(LearnTest.getSeed());
//        Random rand = new Random();
        double p = rand.nextDouble();
        double cumulativeProbability = 0.0;
        for (int i = 0 ; i < probabilities.size() ; i++) {
            cumulativeProbability += probabilities.get(i);
            if (p <= cumulativeProbability) {
                return i;
            }
        }
        return probabilities.size()-1;
    }

    void updateSatCounts(int chainIdx, int gpId, int assignment, int prev_assignment) {
        GroundPredicate gp = state.groundMLN.indexToGroundPredMap.get(gpId);
        for(int formulaId : gp.groundFormulaIds.keySet())
        {
            for(int cid : gp.groundFormulaIds.get(formulaId))
            {
                GroundClause gc = state.groundMLN.groundFormulas.get(formulaId).groundClauses.get(cid);
                int localPredindex = gc.globalToLocalPredIndex.get(gpId);
                BitSet b = gc.grounPredBitSet.get(localPredindex);
                int cur_val = b.get(assignment) ? 1 : 0;
                int prev_val = b.get(prev_assignment) ? 1 : 0;
                int satDifference = cur_val - prev_val; // numsatLiterals according to new assignment - old assignment
                int curSatLiterals = numSatLiterals.get(chainIdx).get(formulaId).get(cid);
                curSatLiterals += satDifference;
                numSatLiterals.get(chainIdx).get(formulaId).set(cid, curSatLiterals);

                if(curSatLiterals > 0)
                {
                    falseClausesSet.get(chainIdx).get(formulaId).remove(cid);
                }
                else
                {
                    falseClausesSet.get(chainIdx).get(formulaId).add(cid);
                }
            }
        }
    }

    public void updateWtsForGndPred(int chainIdx, int gpId) {
        GroundPredicate gp = state.groundMLN.indexToGroundPredMap.get(gpId);
        int numPossibleVals = gp.numPossibleValues;
        double wtsPerVal[] = new double[numPossibleVals];
        Map<Integer, Set<Integer>> formulaIds = gp.groundFormulaIds;

        for(Integer formulaId : formulaIds.keySet())
        {
            GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
            double wt = gf.weight.getValue();
            Set<Integer> tempSet = new HashSet<Integer>();
            tempSet.addAll(falseClausesSet.get(chainIdx).get(formulaId));
            tempSet.removeAll(formulaIds.get(formulaId));
            BitSet formulaBitSet = new BitSet(numPossibleVals);
            formulaBitSet.flip(0,numPossibleVals);

            // If there is a clause which is false, and doesn't contain gp, then this formula is always false
            if(tempSet.size() == 0)
            {

                for(Integer cid : formulaIds.get(formulaId))
                {
                    BitSet clauseBitSet = new BitSet(numPossibleVals);
                    GroundClause gc = gf.groundClauses.get(cid);
                    if(!gc.globalToLocalPredIndex.containsKey(gpId))
                    {
                        System.out.println("problem....");
                    }
                    int localPredIndex = gc.globalToLocalPredIndex.get(gpId);
                    int numSatLiteralsLocal = numSatLiterals.get(chainIdx).get(formulaId).get(cid);
                    if(numSatLiteralsLocal > 1)
                        clauseBitSet.flip(0,numPossibleVals);
                    else if(numSatLiteralsLocal == 1)
                    {
                        BitSet b = gc.grounPredBitSet.get(localPredIndex);
                        if(b.get(truthValues[chainIdx].get(gpId)))
                        {
                            clauseBitSet = (BitSet) b.clone();
                        }
                        else
                        {
                            clauseBitSet.flip(0,numPossibleVals);
                        }
                    }
                    else {
                        BitSet b = gc.grounPredBitSet.get(localPredIndex);
                        clauseBitSet = (BitSet) b.clone();
                    }
                    formulaBitSet.and(clauseBitSet);
                }// end clauses loop

                int startIndex = 0;
                while(startIndex < numPossibleVals)
                {
                    int index = formulaBitSet.nextSetBit(startIndex);
                    if(index == -1)
                        break;
                    wtsPerVal[index] += wt;
                    startIndex = index+1;
                }
            }// end if condition

        } //end formulas loops

        List<Double> tempWts = new ArrayList<>();
        for(double d : wtsPerVal)
        {
            tempWts.add(d);
        }
        wtsPerPredPerVal.get(chainIdx).put(gpId, tempWts);
    }

    /**
     * Update trueCounts such as formulaTrueCnts, formulaTrueSqCnts, and allFormulaTrueCnts for this chain
     * @param chainIdx Chain for which this updation needs to be done.
     */
    void updateTrueCnts(int chainIdx)
    {
        int numCounts = formulaTrueCnts.length;
        double []numTrueGndings = getNumTrueGndingsFromTruthVals(truthValues[chainIdx]);
        if(saveAllCounts)
        {
            allFormulaTrueCnts.add(new ArrayList<Double>(Collections.nCopies(numCounts,0.0)));
        }
        for (int i = 0; i < numCounts; i++) {
            if(saveAllCounts)
            {
                allFormulaTrueCnts.get(numSamples).set(i, numTrueGndings[i]);
            }
            formulaTrueCnts[i] += numTrueGndings[i];
            formulaTrueSqCnts[i] += numTrueGndings[i]*numTrueGndings[i];
        }
        numSamples++;
    }

    private double[] getNumTrueGndingsFromTruthVals(Map<Integer, Integer> truthValues)
    {
        double numTrueGndings[] = new double[numWts];
        for(GroundFormula gf : state.groundMLN.groundFormulas)
        {
            boolean isFormulaSatisfied = true;
            for(GroundClause gc : gf.groundClauses)
            {
                boolean isClauseSatisfied = false;
                for(Integer gpId : gc.groundPredIndices)
                {
                    BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                    int trueVal = truthValues.get(gpId);
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
                    numTrueGndings[numWts] += gf.originalWeight.getValue();
                }
            }
        }
        return numTrueGndings;
    }

    public void writeProbs(PrintWriter writer)
    {
        for(Integer g : state.groundMLN.indexToGroundPredMap.keySet())
        {
            int numVals = state.groundMLN.indexToGroundPredMap.get(g).numPossibleValues;
            for(int val = 0 ; val < numVals ; val++)
            {
                double marginal = 0.0;
                for (int c = 0; c < numChains; c++) {
                    marginal += numValPerChainPerPred_.get(c).get(g).get(val);
                }
                double marginal_prob = marginal/numChains;
                writer.printf(state.groundMLN.indexToGroundPredMap.get(g) + "=" + val + " %.4f\n",marginal_prob);
            }
        }
    }

}
