package org.utd.cs.mln.inference;

import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.State;
import org.utd.cs.mln.alchemy.util.ConvergenceTest;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Implements Gibbs Sampling
 * @author Happy
 * @since 03/28/18
 */
public class GibbsSampler_v3 extends MCMC{

    // Gamma used by convergence test
    double gamma;
    // Epsilon used by convergence test
    double epsilonError;
    // Fraction of samples needed to converge
    double fracConverged;
    // If true, test for convergence, otherwise do not test
    boolean testConvergence;
    // Number of samples between checking for convergence
    int samplesPerTest;

    // Convergence test for burning in
    private GelmanConvergenceTest[] burnConvergenceTests;

    // Convergence test for sampling
    private ConvergenceTest[] gibbsConvergenceTests;

    /**
     * Constructor: User-set parameters are set.
     * @see MCMC#MCMC(State, int, boolean, boolean, MCMCParams)
     */
    public GibbsSampler_v3(State state, int seed, boolean trackFormulaTrueCnts, boolean priorSoftEvidence, GibbsParams gibbsParams) {
        super(state, seed, trackFormulaTrueCnts, priorSoftEvidence, gibbsParams);
        // User-set parameters
        this.gamma = gibbsParams.gamma;
        this.epsilonError = gibbsParams.epsilonError;
        this.fracConverged = gibbsParams.fracConverged;
        this.testConvergence = gibbsParams.testConvergence;
        this.samplesPerTest = gibbsParams.samplesPerTest;
    }

    @Override
    void init() {
        initTruthValsAndSatTrueLits();
        if(mcmcdebug)
            System.out.println("initializing Gibbs sampling randomly");
        randomInitializeTruthVals();
        initializeNumSatValues();
        // Initialize convergence test
        int numGndPreds = state.groundMLN.groundPredicates.size();
        initConvergenceTests(gamma, epsilonError, numGndPreds, numChains);
    }

    /**
     * Initializes convergence tests for burning in and sampling.
     */
    void initConvergenceTests(double gamma, double epsilonFrac,
                              int numGndPreds, int numChains)
    {
        burnConvergenceTests = new GelmanConvergenceTest[numGndPreds];
        gibbsConvergenceTests = new ConvergenceTest[numGndPreds];
        for (int i = 0; i < numGndPreds; i++)
        {
            burnConvergenceTests[i]  = new GelmanConvergenceTest(numChains);
            gibbsConvergenceTests[i] = new ConvergenceTest(numChains, gamma, epsilonFrac);
        }
    }

    @Override
    void infer() {
        initNumValPerChainPerPred();
        // Burn-in only if burnMaxSteps positive
        boolean burningIn = (burnMaxSteps > 0) ? true : false;

        // If keeping track of true clause groundings, then init to zero
        if(trackFormulaTrueCnts)
        {
            for (int clauseId = 0; clauseId < formulaTrueCnts.length; clauseId++) {
                formulaTrueCnts[clauseId] = 0;
            }
        }
        ArrayList<Integer> affectedGndPredIndices = new ArrayList<>();
        int numGndPreds = state.groundMLN.groundPredicates.size();

        for(int g = 0 ; g < numGndPreds ; g++) {
            affectedGndPredIndices.add(g);
        }

        for (int c = 0; c < numChains; c++) {
            updateWtsForGndPreds(affectedGndPredIndices, c);
        }
        affectedGndPredIndices.clear();

        System.out.println("Running Gibbs Sampling...");

        // Sampling loop
        int sample = 0;
        int numSamplesPerChainPerPred = 0;
        boolean done = false;
        long time = System.currentTimeMillis();
        double secondsElapsed = 0;
        while(!done)
        {
            ++sample;
            if(sample % samplesPerTest == 0)
            {
                secondsElapsed = (System.currentTimeMillis() - time) / 1000.0;
                System.out.println("Sample (per pred per chain) : " + sample + ", Elapsed Time : " + Timer.time(secondsElapsed));
            }

            GibbsPerChain gibbsPerChain;
//            Thread t[] = new Thread[numChains];
//            Thread.currentThread().setPriority(10);
            // For each chain, for each node, generate the node's new truth value
            for (int c = 0; c < numChains; c++) {
                gibbsPerChain = new GibbsPerChain(c, burningIn, trackFormulaTrueCnts);
                gibbsPerChain.run();
//                t[c] = new Thread(gibbsPerChain);
//                t[c].setPriority(1);
            }
//            try {
//                for (int c = 0; c < numChains; c++) {
//                    t[c].start();
//                }
//                //System.out.println("Waiting for numChains threads to finish.");
//                Thread.currentThread().setPriority(1);
//                for (int c = 0; c < numChains; ++c) {
//                    t[c].join();
//                }
//                Thread.currentThread().setPriority(10);
//            }
//            catch (InterruptedException e) {
//                System.out.println("Main thread Interrupted");
//            }
            if (!burningIn) numSamplesPerChainPerPred++;


            // Add current truth values to the convergence testers
            for (int g = 0; g < state.groundMLN.groundPredicates.size(); g++)
            {
                double vals[] = new double[numChains];
                for (int c = 0; c < numChains; c++) {
                    vals[c] = truthValues[c][g];
                }
                if (burningIn)
                    burnConvergenceTests[g].appendNewValues(vals);
                else
                    gibbsConvergenceTests[g].appendNewValues(vals);
            }

            if (sample % samplesPerTest != 0) continue;

            if (burningIn)
            {
                // Use convergence criteria stated in "Probability and Statistics",
                // DeGroot and Schervish
                boolean burnConverged = false;

                if (testConvergence)
                    burnConverged =
                            GelmanConvergenceTest.checkConvergenceOfAll(burnConvergenceTests, state.groundMLN.groundPredicates.size(), true);
                if ( (sample >= burnMinSteps && burnConverged)
                        || (burnMaxSteps >= 0 && sample >= burnMaxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds))
                {
                    System.out.println("Done burning. " + sample + " samples per pred per chain");
                    if (testConvergence)
                    {
                        System.out.println(" (" + (burnConverged? "converged":"didn't converge")
                                +" at total of " + numChains*sample + " samples per pred)");
                    }
                    burningIn = false;
                    sample = 0;
                }
            }

            else
            {  // Doing actual gibbs sampling
                boolean gibbsConverged = false;

                if (testConvergence)
                    gibbsConverged = ConvergenceTest.checkConvergenceOfAtLeast(gibbsConvergenceTests, state.groundMLN.groundPredicates.size(), sample, fracConverged, true);

                if (   (sample >= minSteps && gibbsConverged)
                        || (maxSteps >= 0 && sample >= maxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds))
                {
                    System.out.println("Done Gibbs sampling. " + sample
                            + " samples per pred per chain");
                    if (testConvergence)
                    {
                        System.out.println(" (" + (gibbsConverged? "converged":"didn't converge")
                                +" at total of " + numSamplesPerChainPerPred + " samples per chain per pred)");
                    }
                    done = true;
                }
            }
        }
        System.out.println("Time taken for Gibbs sampling : "+Timer.time((System.currentTimeMillis() - time) / 1000.0));
        // update gndPreds probability
        for (int c = 0; c < numChains; c++) {
            for (int g = 0; g < numValPerChainPerPred_.get(c).size(); g++) {
                for(int val = 0 ; val < numValPerChainPerPred_.get(c).get(g).size() ; val++)
                {
                    double p = numValPerChainPerPred_.get(c).get(g).get(val)/numSamplesPerChainPerPred;
                    numValPerChainPerPred_.get(c).get(g).set(val,p);
                }
            }
        }

    }


    /**
     * Performs one step of Gibbs sampling in one chain.
     *
     * @param chainIdx Index of chain in which the Gibbs step is performed.
     * @param burningIn If true, burning-in is occuring. Otherwise, false.
     */
    public void performGibbsStep(int chainIdx, boolean burningIn)
    {
        int numGndPreds = state.groundMLN.groundPredicates.size();
        for(int gpId = 0; gpId < numGndPreds; gpId++){
            int assignment = performGibbsStep(chainIdx,gpId, (gpId+1)%numGndPreds);
            // If in actual gibbs sampling phase, update numValPerPred
            if(!burningIn)
                numValPerChainPerPred_.get(chainIdx).get(gpId).set(assignment, numValPerChainPerPred_.get(chainIdx).get(gpId).get(assignment)+1);
        }
    }

    /**
     * Performs one step of Gibbs sampling for predicate at gpId.
     * @param chainIdx Index of chain in which the Gibbs step is performed.
     * @param gpId PredicateId for which we try to flip
     * @param nextGpId PredicateId which will be affected if we flip gpId
     * @return assignment by gibbs step of gpId
     */
    int performGibbsStep(int chainIdx, int gpId, int nextGpId)
    {
        int assignment = get_probabilistic_assignment(wtsPerPredPerVal.get(chainIdx).get(gpId));
        int prev_assignment = truthValues[chainIdx][gpId];
        truthValues[chainIdx][gpId] = assignment;
        if(assignment != prev_assignment)
        {
            updateSatCounts(chainIdx, gpId, assignment, prev_assignment);
            updateWtsForGndPred(chainIdx,nextGpId);
        }
        return assignment;
    }

    private class GibbsPerChain implements Runnable{
        private final int chainIdx;
        private final boolean burningIn, trackFormulaTrueCnts;
        public GibbsPerChain(int chainIdx, boolean burningIn, boolean trackFormulaTrueCnts) {
            this.chainIdx = chainIdx;
            this.burningIn = burningIn;
            this.trackFormulaTrueCnts = trackFormulaTrueCnts;
        }

        @Override
        public void run() {
            GibbsSampler_v3.this.performGibbsStep(chainIdx, burningIn);
            if(!burningIn && trackFormulaTrueCnts)
            {
                updateTrueCnts(chainIdx);
            }
        }
    }
}