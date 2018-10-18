package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.SimpleGroundClause;
import org.utd.cs.mln.alchemy.core.SimpleGroundMln;
import org.utd.cs.mln.alchemy.core.SimpleState;
import org.utd.cs.mln.alchemy.core.State;

import java.io.PrintWriter;
import java.util.*;

import static org.utd.cs.mln.alchemy.util.OtherUtils.getUniformAssignment;

/**
 * Created by Happy on 8/16/18.
 */
public class SimpleMCSAT {

    public SimpleGroundMln simpleGroundMln;
    private final int DEFAULT_SEED = 2350877;
    // Number of samples after which print info
    private int samplesPerTest;

    // whether we want to calculate marginal of each predicate or not
    boolean calculateMarginal;

    int seed;

    // Indicates if need to store formulaTrueCnts and formulaTrueSqCnts or not
    boolean trackFormulaTrueCnts;
    boolean saveAllCounts;
    MCMCParams params;

    // Expected true counts and true squared counts of first order formulas.
    // If softEvidence, then size is +1
    // Note that we are not taking them to be List because we want to increment the counts regularly, and doing that
    // is cumbersome in Lists.
    public double [] formulaTrueCnts, formulaTrueSqCnts;

    // For each ground pred g, holds number of times g is set to true
    List<Double> numTruePerPred_;

    Random rand;

    //Indicates whether softEvidence is given or not
    boolean priorSoftEvidence;

    // numWts is equal to number of first order formulas (+1 if softEvidence)
    int numWts;

    /**
     * Constructor: User-set parameters are set.
     *
     * @param sgmln
     * @param seed
     * @param trackFormulaTrueCnts
     * @param saveAllCounts
     * @param priorSoftEvidence
     * @param params
     * @param calculateMarginal
     * @see Inference#Inference(State, int, boolean, boolean, boolean)
     */
    public SimpleMCSAT(SimpleGroundMln sgmln, int seed, boolean trackFormulaTrueCnts, boolean saveAllCounts, boolean priorSoftEvidence, MCMCParams params, boolean calculateMarginal) {
        this.simpleGroundMln = sgmln;
        this.seed = seed;
        if(seed == -1)
            this.seed = DEFAULT_SEED;
        this.trackFormulaTrueCnts = trackFormulaTrueCnts;
        this.saveAllCounts = saveAllCounts;
        this.priorSoftEvidence = priorSoftEvidence;
        this.params = params;
        this.calculateMarginal = calculateMarginal;
        this.rand = new Random();
        numWts = simpleGroundMln.mln.formulas.size();
        if(priorSoftEvidence)
            numWts++;

        if(trackFormulaTrueCnts)
        {
            formulaTrueCnts = new double[numWts];
            formulaTrueSqCnts = new double[numWts];
        }
        numTruePerPred_ = new ArrayList<>(Collections.nCopies(simpleGroundMln.groundPreds.size(),0.0));
    }

    public void infer() {
        // First randomly initialize state
        int numGroundPreds = simpleGroundMln.groundPreds.size();

        boolean burningIn = true;
        long time = System.currentTimeMillis();
        double secondsElapsed = 0;

        // This stores indices of clauses which are selected for samplesat step in the previous state.
        Set<Integer> samplesatClauses = new HashSet<>();

        SimpleState prevState = new SimpleState(simpleGroundMln);
        for(int i = 1 ; i < numGroundPreds ; i++)
        {
            int assignment = getUniformAssignment(2);
            prevState.world.set(i,assignment);
        }
        for(int numSamples = 1 ; numSamples <= params.burnMaxSteps ; numSamples++)
        {
            if(numSamples % samplesPerTest == 0)
            {
                secondsElapsed = (System.currentTimeMillis() - time) / 1000.0;
                System.out.println("Sample : " + numSamples + ", Elapsed Time : " + org.utd.cs.gm.utility.Timer.time(secondsElapsed));
            }
            SimpleState newState = performMCSatStep(prevState, samplesatClauses);
            prevState = newState;
        }
    }

    SimpleState performMCSatStep(SimpleState prevState, Set<Integer> samplesatClauses)
    {
        // List of positive weighted clauses selected for samplesat
        List<SimpleGroundClause> M = new ArrayList<>();
        // For each ground formula, check if it is satisfied or not. It it is, add it to list M with some prob p = (1-e^-w)
        for(int clauseId = 0 ; clauseId < simpleGroundMln.groundClauses.size() ; clauseId++)
        {
            // If this ground formula is already present in samplesatClauses, then
            // no need to check for satisfiability, just add it with prob p (after making the weight positive), because in this state,
            // all the formulas in samplesatClauses are already satisfied since this state
            // comes from last sample sat step.
            if(samplesatClauses.contains(clauseId))
            {
                addClauseIntoM(clauseId, M, samplesatClauses);

            }
            // else check for the satisfiability of this clause, and then add into M with prob p
            boolean isSatisfied = false;
            SimpleGroundClause sgc = simpleGroundMln.groundClauses.get(clauseId);
            for(Integer atomIndex : sgc.groundAtomIndices)
            {
                int sign = 1;
                if(atomIndex < 0)
                {
                    sign = -1;
                }
                int assignment = prevState.world.get(Math.abs(atomIndex));
                if((sign == 1 && assignment == 1) || (sign == -1 && assignment == 0))
                {
                    isSatisfied = true;
                    break;
                }
            }
            if(isSatisfied)
            {
                addClauseIntoM(clauseId, M, samplesatClauses);
            }
        }
        SimpleState newState = performSampleSatStep(M, prevState);
        return newState;
    }

    /**
     *
     * @param M Set of clauses (unweighted) which are required to be satisfied.
     * @param prevState Previous state of MC-SAT. If we can't find a new state in samplesat, then return this state.
     * @return A new uniformly sampled state
     *
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
     */
    private SimpleState performSampleSatStep(List<SimpleGroundClause> M, SimpleState prevState) {
        SimpleState resultState = new SimpleState(prevState.simpleGroundMln);
        resultState.randomize();

        // Find initial number of satisfied clauses
        int numSatClauses = 0;

        return resultState;
    }

    private void addClauseIntoM(int clauseId, List<SimpleGroundClause> M, Set<Integer> samplesatClauses) {
        SimpleGroundClause sgc = simpleGroundMln.groundClauses.get(clauseId);
        double p = (1-Math.exp(-Math.abs(sgc.weight)));
        if(rand.nextDouble() < p)
        {
            // If weight is positive, then add this clause into M as it is
            if(sgc.weight >= 0)
            {
                SimpleGroundClause newsgc = new SimpleGroundClause();
                newsgc.groundAtomIndices = new ArrayList<>(sgc.groundAtomIndices);
                M.add(newsgc);
            }
            // else separate the atoms and add them negated
            else
            {
                for(Integer atomIndex : sgc.groundAtomIndices)
                {
                    SimpleGroundClause newsgc = new SimpleGroundClause();
                    newsgc.groundAtomIndices.add(-1*atomIndex);
                    M.add(newsgc);
                }
            }

        }
        // else we don't put this clause into M, and hence remove this from samplesatclauses
        else
        {
            samplesatClauses.remove(clauseId);
        }
    }

    public void writeProbs(PrintWriter writer) {
    }
}
