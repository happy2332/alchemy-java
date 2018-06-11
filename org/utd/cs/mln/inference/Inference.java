package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.State;

import java.io.PrintWriter;
import java.util.*;

/**
 * Abstract class from which all inference algorithms are derived.
 * At least one method is abstract making this an abstract class
 * (it can not be instantiated).
 * @author Happy
 * @since 03/28/18
 * @see MCMC
 */
public abstract class Inference {
    private final int DEFAULT_SEED = 2350877;

    // whether we want to calculate marginal of each predicate or not
    boolean calculateMarginal;

    public State state;

    int seed;

    // Indicates if need to store formulaTrueCnts and formulaTrueSqCnts or not
    boolean trackFormulaTrueCnts;

    // Expected true counts and true squared counts of first order formulas.
    // If softEvidence, then size is +1
    // Note that we are not taking them to be List because we want to increment the counts regularly, and doing that
    // is cumbersome in Lists.
    public double [] formulaTrueCnts, formulaTrueSqCnts;

    // For each ground pred g, holds number of times g is set to a value val
    // numValPerPred[g][val]
    Map<Integer, List<Double>> numValPerPred_;

    Random rand;

    //Indicates whether softEvidence is given or not
    boolean priorSoftEvidence;

    // numWts is equal to number of first order formulas (+1 if softEvidence)
    int numWts;

    /**
     * Constructor: Every inference algorithm is required to have a VariableState
     * representing the state of variables and clauses and a seed for any
     * randomization in the algorithm. If there is no randomization, seed is not
     * used.
     *
     * @param state State of the variables and clauses of the inference.
     * @param seed Seed used to initialize randomization in the algorithm.
     * @param trackFormulaTrueCnts Indicates if need to store true counts for each first-order formula
     */
    public Inference(State state, int seed, boolean trackFormulaTrueCnts, boolean priorSoftEvidence, boolean calculateMarginal)
    {
        this.state = state;
        this.seed = seed;
        if(seed == -1)
            this.seed = DEFAULT_SEED;
        this.trackFormulaTrueCnts = trackFormulaTrueCnts;
        this.priorSoftEvidence = priorSoftEvidence;
        this.calculateMarginal = calculateMarginal;
        this.rand = new Random();
        numWts = state.mln.formulas.size();
        if(priorSoftEvidence)
            numWts++;

        if(trackFormulaTrueCnts)
        {
            formulaTrueCnts = new double[numWts];
            formulaTrueSqCnts = new double[numWts];
        }

        numValPerPred_ = new HashMap<>();
        for(Integer g : state.groundMLN.indexToGroundPredMap.keySet())
        {
            int numVals = state.groundMLN.indexToGroundPredMap.get(g).numPossibleValues;
            numValPerPred_.put(g, new ArrayList<Double>(Collections.nCopies(numVals, 0.0)));
        }
    }

    /**
     * Initializes the inference algorithm.
     */
    abstract void init();

    /**
     * Performs the inference algorithm.
     */
    abstract public void infer();

    /**
     * Prints the probabilities of each predicate.
     */
    abstract void writeProbs(PrintWriter writer);

    /**
     * Resets all counts
     */
    public void resetCnts() {
        if(trackFormulaTrueCnts){
            Arrays.fill(formulaTrueCnts,0.0);
            Arrays.fill(formulaTrueSqCnts,0.0);
        }
        if(calculateMarginal)
        {
            for(int gpId : numValPerPred_.keySet()) {
                int numPossibleVals = numValPerPred_.get(gpId).size();
                numValPerPred_.put(gpId, new ArrayList<>(Collections.nCopies(numPossibleVals, 0.0)));
            }
        }

    }


    public void writeNetwork(PrintWriter writer) {
        writer.write(state.groundMLN.writableString());
    }
}
