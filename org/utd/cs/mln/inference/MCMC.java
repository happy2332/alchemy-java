package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.OtherUtils;

import java.io.PrintWriter;
import java.util.*;

import static org.utd.cs.mln.alchemy.util.OtherUtils.getUniformAssignment;

/**
 * Abstract class from which all MCMC algorithms are sub classed.
 * @author Happy
 * @since 03/28/18
 * @see MCMCParams,GibbsSampler_v3
 */
public abstract class MCMC extends Inference {
    // No. of chains which MCMC will use
    public int numChains;
    // Min. no. of burn-in steps MCMC will take per chain
    public int burnMinSteps;
    // Max. no. of burn-in steps MCMC will take per chain
    public int burnMaxSteps;
    // Min. no. of sampling steps MCMC will take per chain
    public int minSteps;
    // Max. no. of sampling steps MCMC will take per chain
    public int maxSteps;
    // Max. no. of seconds MCMC should run
    public int maxSeconds;

    boolean mcmcdebug = true;
    boolean trackState = false;

    // Truth values in each chain for each ground predicate (truthValues[c][p])
    public Map<Integer, Integer> truthValues [];

    // Num. of satisfying literals in each chain for each ground formula for each ground clause
    // numSatLiterals[chain][formula][clause]
    List<List<List<Integer>>> numSatLiterals;

    List<Set<Integer>> gpIdsToBeChanged;

    // For each chain, For each formula, set of clauseIds which are false in this truthVals configuration
    // falseClausesSet[chain][formula][
    List<List<Set<Integer>>> falseClausesSet;

    // For each chain, Holds number of times a ground predicate is set to a value
    // numValPerPred[c][g][val]
    List<Map<Integer, List<Double>>> numValPerChainPerPred_;

    // For each chain, for each ground predicate, for each possible value, stores satweight
    // wtsPerPredVal[c][g][val]
    List<Map<Integer, List<Double>>> wtsPerPredPerVal;

    // For each chain c, allFormulaTrueCntsPerChain[c][i][j] is the number of true groundings of jth formula in ith sample
    List<List<List<Double>>> allFormulaTrueCntsPerChain;

    // allFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample
    // oldAllFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample in the previous iter of learning
    public List<List<Double>> allFormulaTrueCnts, oldAllFormulaTrueCnts;

    // Whether to save all counts for each sample
    boolean saveAllCounts;



    /**
     * Constructor: User-set parameters are set.
     * @see Inference#Inference(State, int, boolean, boolean, boolean)
     */
    public MCMC(State state, int seed, boolean trackFormulaTrueCnts, boolean saveAllCounts, boolean priorSoftEvidence, MCMCParams params, boolean calculateMarginal) {
        super(state, seed, trackFormulaTrueCnts, priorSoftEvidence, calculateMarginal);
        this.numChains = params.numChains;
        this.burnMinSteps = params.burnMinSteps;
        this.burnMaxSteps = params.burnMaxSteps;
        this.minSteps = params.minSteps;
        this.maxSteps = params.maxSteps;
        this.maxSeconds = params.maxSeconds;
        this.saveAllCounts = saveAllCounts;

        truthValues = new Map[numChains];
        numSatLiterals = new ArrayList<>();
        gpIdsToBeChanged = new ArrayList<>();
        falseClausesSet = new ArrayList<>();
        numValPerChainPerPred_ = new ArrayList<>();
        wtsPerPredPerVal = new ArrayList<>();
        allFormulaTrueCntsPerChain = new ArrayList<>();
        allFormulaTrueCnts = new ArrayList<>();
        oldAllFormulaTrueCnts = new ArrayList<>();

        int numFormulas = state.groundMLN.groundFormulas.size();
        for (int c = 0; c < numChains; c++) {
            numSatLiterals.add(new ArrayList<List<Integer>>());
            gpIdsToBeChanged.add(new HashSet<Integer>());
            wtsPerPredPerVal.add(new HashMap<Integer, List<Double>>());
            falseClausesSet.add(new ArrayList<Set<Integer>>());
            numValPerChainPerPred_.add(new HashMap<Integer, List<Double>>());
            allFormulaTrueCntsPerChain.add(new ArrayList<List<Double>>());
            for (int formulaId = 0; formulaId < numFormulas; formulaId++) {
                int numClauses = state.groundMLN.groundFormulas.get(formulaId).groundClauses.size();
                numSatLiterals.get(c).add(new ArrayList<>(Collections.nCopies(numClauses,0)));
                falseClausesSet.get(c).add(new HashSet<Integer>());
            }
            truthValues[c] = new HashMap<>();
            for (GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet()) {
                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
                truthValues[c].put(gpId, 0);
                int numPossibleVals = gp.numPossibleValues;
                numValPerChainPerPred_.get(c).put(gpId, new ArrayList<>(Collections.nCopies(numPossibleVals,0.0)));
                wtsPerPredPerVal.get(c).put(gpId, new ArrayList<Double>(Collections.nCopies(numPossibleVals,0.0)));
            }
        }

        if(trackState)
        {

        }

    }

//    /**
//     * Allocates memory to all datamembers of MCMC : truthValues, numSatLiterals, falseClauseSet, numValPerChainPerPred_,
//     * wtsPerPredPerVal, formulaTrueCntsPerChain, formulaTrueSqCntsPerChain, and allFormulaTrueCntsPerChain.
//     */
//    void init() {
//        truthValues = new Map[numChains];
//        numSatLiterals = new ArrayList<>();
//        falseClausesSet = new ArrayList<>();
//        numValPerChainPerPred_ = new ArrayList<>();
//        wtsPerPredPerVal = new ArrayList<>();
//        allFormulaTrueCntsPerChain = new ArrayList<>();
//        allFormulaTrueCnts = new ArrayList<>();
//        oldAllFormulaTrueCnts = new ArrayList<>();
//        formulaTrueCntsPerChain = new double[numChains][numWts];
//        formulaTrueSqCntsPerChain = new double[numChains][numWts];
//
//        int numFormulas = state.groundMLN.groundFormulas.size();
//        for (int c = 0; c < numChains; c++) {
//            numSatLiterals.add(new ArrayList<List<Integer>>());
//            wtsPerPredPerVal.add(new HashMap<Integer, List<Double>>());
//            falseClausesSet.add(new ArrayList<Set<Integer>>());
//            numValPerChainPerPred_.add(new HashMap<Integer, List<Double>>());
//            allFormulaTrueCntsPerChain.add(new ArrayList<List<Double>>());
//            for (int formulaId = 0; formulaId < numFormulas; formulaId++) {
//                int numClauses = state.groundMLN.groundFormulas.get(formulaId).groundClauses.size();
//                numSatLiterals.get(c).add(new ArrayList<>(Collections.nCopies(numClauses,0)));
//                falseClausesSet.get(c).add(new HashSet<Integer>());
//            }
//            truthValues[c] = new HashMap<>();
//            for (GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet()) {
//                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
//                truthValues[c].put(gpId, 0);
//                int numPossibleVals = gp.numPossibleValues;
//                numValPerChainPerPred_.get(c).put(gpId, new ArrayList<>(Collections.nCopies(numPossibleVals,0.0)));
//                wtsPerPredPerVal.get(c).put(gpId, new ArrayList<Double>(Collections.nCopies(numPossibleVals,0.0)));
//            }
//        }
//    }

    /**
     * Computes Hessian vector product. Hessian is approximated from the matrix allFormulaTrueCnts.
     * @param v vector to which Hessian is to be multiplied
     * @return resultant vector
     */
    public double[] getHessianVectorProduct(double[] v) {
        // For minimizing the negative log likelihood,
        // the ith element of H v is:
        //   E[n_i * vn] - E[n_i] E[vn]
        // where n is the vector of all clause counts
        // and vn is the dot product of v and n.

        double sumVN = 0;
        double []sumN = new double[numWts];
        double []sumNiVN = new double[numWts];

        // Get sufficient statistics from each sample,
        // so we can compute expectations
        int numSamples = allFormulaTrueCnts.size();
        for (int s = 0; s < numSamples; s++)
        {
            List<Double> n1 = allFormulaTrueCnts.get(s);

            // Compute v * n
            double vn = 0;

            for (int i = 0; i < numWts; i++)
                vn += v[i] * n1.get(i);

            // Tally v*n, n_i, and n_i v*n
            sumVN += vn;
            for (int i = 0; i < numWts; i++)
            {
                sumN[i]    += n1.get(i);
                sumNiVN[i] += n1.get(i) * vn;
            }
        }

        // Compute actual product from the sufficient stats
        double []product = new double[numWts];
        for (int formulano = 0; formulano < numWts; formulano++)
        {
            double E_vn = sumVN/numSamples;
            double E_ni = sumN[formulano]/numSamples;
            double E_nivn = sumNiVN[formulano]/numSamples;
            product[formulano] = E_nivn - E_ni * E_vn;
        }

        return product;
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
        double []numTrueGndings = getNumTrueGndingsFromTruthVals(truthValues[chainIdx]);
        int numWts = numTrueGndings.length;
        if(saveAllCounts)
        {
            allFormulaTrueCntsPerChain.get(chainIdx).add(new ArrayList<Double>(Collections.nCopies(numWts,0.0)));
        }
        for (int i = 0; i < numWts; i++) {
            int size = allFormulaTrueCntsPerChain.get(chainIdx).size();
            if(saveAllCounts)
            {
                allFormulaTrueCntsPerChain.get(chainIdx).get(size-1).set(i, numTrueGndings[i]);
            }
        }
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
                writer.printf(state.groundMLN.indexToGroundPredMap.get(g) + "=" + val + " %.4f\n",numValPerPred_.get(g).get(val));
            }
        }
    }

    @Override
    /**
     * Reset all data structures before starting inference
     */
    public void resetCnts(){
        super.resetCnts();
        truthValues = new Map[numChains];
        numSatLiterals = new ArrayList<>();
        falseClausesSet = new ArrayList<>();
        numValPerChainPerPred_ = new ArrayList<>();
        wtsPerPredPerVal = new ArrayList<>();
        allFormulaTrueCntsPerChain = new ArrayList<>();
        allFormulaTrueCnts = new ArrayList<>();

        int numFormulas = state.groundMLN.groundFormulas.size();
        for (int c = 0; c < numChains; c++) {
            numSatLiterals.add(new ArrayList<List<Integer>>());
            wtsPerPredPerVal.add(new HashMap<Integer, List<Double>>());
            falseClausesSet.add(new ArrayList<Set<Integer>>());
            numValPerChainPerPred_.add(new HashMap<Integer, List<Double>>());
            allFormulaTrueCntsPerChain.add(new ArrayList<List<Double>>());
            for (int formulaId = 0; formulaId < numFormulas; formulaId++) {
                int numClauses = state.groundMLN.groundFormulas.get(formulaId).groundClauses.size();
                numSatLiterals.get(c).add(new ArrayList<>(Collections.nCopies(numClauses,0)));
                falseClausesSet.get(c).add(new HashSet<Integer>());
            }
            truthValues[c] = new HashMap<>();
            for (GroundPredicate gp : state.groundMLN.groundPredToIntegerMap.keySet()) {
                int gpId = state.groundMLN.groundPredToIntegerMap.get(gp);
                truthValues[c].put(gpId, 0);
                int numPossibleVals = gp.numPossibleValues;
                numValPerChainPerPred_.get(c).put(gpId, new ArrayList<>(Collections.nCopies(numPossibleVals,0.0)));
                wtsPerPredPerVal.get(c).put(gpId, new ArrayList<Double>(Collections.nCopies(numPossibleVals,0.0)));
            }
        }
    }

    /**
     * If we are storing all counts, then restore all counts from old counts
     */
    public void restoreCnts() {
        if (!saveAllCounts)
            return;

        //resetCnts();
        for (int i = 0; i < allFormulaTrueCnts.size(); i++) {
            allFormulaTrueCnts.get(i).clear();
        }
        allFormulaTrueCnts.clear();
        System.out.println("Restoring counts...");
        if(trackFormulaTrueCnts){
            for (int i = 0; i < oldAllFormulaTrueCnts.size(); i++)
            {
                allFormulaTrueCnts.add(new ArrayList<Double>());
                int numcounts = oldAllFormulaTrueCnts.get(i).size();
                for (int j = 0; j < numcounts; j++)
                {
                    allFormulaTrueCnts.get(i).add(oldAllFormulaTrueCnts.get(i).get(j));
//                    formulaTrueCnts[j] += allFormulaTrueCnts.get(i).get(j);
//                    formulaTrueSqCnts[j] += allFormulaTrueCnts.get(i).get(j) * allFormulaTrueCnts.get(i).get(j);
                }
            }
        }
    }

    /**
     * If we are storing all counts, save current all counts to old all counts
     */
    public void saveToOldCnts() {
        if (!saveAllCounts)
            return;
        for (int i = 0; i < oldAllFormulaTrueCnts.size(); i++) {
            oldAllFormulaTrueCnts.get(i).clear();
        }
        oldAllFormulaTrueCnts.clear();
        if(trackFormulaTrueCnts){
            for (int i = 0; i < allFormulaTrueCnts.size(); i++)
            {
                oldAllFormulaTrueCnts.add(new ArrayList<Double>());
                int numcounts = allFormulaTrueCnts.get(i).size();
                for (int j = 0; j < numcounts; j++)
                {
                    oldAllFormulaTrueCnts.get(i).add(allFormulaTrueCnts.get(i).get(j));
                }
            }
        }
    }
}
