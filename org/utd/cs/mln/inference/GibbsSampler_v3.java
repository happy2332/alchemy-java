package org.utd.cs.mln.inference;

import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.GroundFormula;
import org.utd.cs.mln.alchemy.core.GroundMLN;
import org.utd.cs.mln.alchemy.core.GroundPredicate;
import org.utd.cs.mln.alchemy.core.State;
import org.utd.cs.mln.alchemy.util.ConvergenceTest;

import java.util.ArrayList;
import java.util.Arrays;
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
     * @see MCMC#MCMC(State, int, boolean, boolean, boolean, MCMCParams, boolean)
     */
    public GibbsSampler_v3(State state, int seed, boolean trackFormulaTrueCnts, boolean saveAllCounts, boolean priorSoftEvidence, GibbsParams gibbsParams, boolean calculateMarginal) {
        super(state, seed, trackFormulaTrueCnts, saveAllCounts, priorSoftEvidence, gibbsParams, calculateMarginal);
        // User-set parameters
        this.gamma = gibbsParams.gamma;
        this.epsilonError = gibbsParams.epsilonError;
        this.fracConverged = gibbsParams.fracConverged;
        this.testConvergence = gibbsParams.testConvergence;
        this.samplesPerTest = gibbsParams.samplesPerTest;
        init();
    }

    @Override
    void init() {
        if(mcmcdebug)
            System.out.println("initializing Gibbs sampling randomly");
        randomInitializeTruthVals();
        initializeNumSatValues();
        //findMarkovBlankets();
        // Initialize convergence test
        int numGndPreds = state.groundMLN.indexToGroundPredMap.size();
        initConvergenceTests(gamma, epsilonError, numGndPreds, numChains);
    }

//    private void findMarkovBlankets() {
//        GroundMLN groundMln = state.groundMLN;
//        for(GroundFormula gf : groundMln.groundFormulas)
//        {
//            for(int gpId : gf.groundPredIndices)
//            {
//                GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpId);
//                gp.markovBlanketSet.addAll(gf.groundPredIndices);
//            }
//        }
//    }

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
    public void infer() {
        // Burn-in only if burnMaxSteps positive
        boolean burningIn = (burnMaxSteps > 0) ? true : false;

//        // If keeping track of true clause groundings, then init to zero
//        if(trackFormulaTrueCnts)
//        {
//            for (int clauseId = 0; clauseId < formulaTrueCnts.length; clauseId++) {
//                formulaTrueCnts[clauseId] = 0;
//            }
//        }
        ArrayList<Integer> affectedGndPredIndices = new ArrayList<>();
        affectedGndPredIndices.addAll(state.groundMLN.indexToGroundPredMap.keySet());

        for (int c = 0; c < numChains; c++) {
            updateWtsForGndPreds(affectedGndPredIndices, c);
        }
        affectedGndPredIndices.clear();

        System.out.println("Running Gibbs Sampling...");

        // Sampling loop
        int numSamplesPerChainPerPred = 0;
        boolean done = false;
        long time = System.currentTimeMillis();
        double secondsElapsed = 0;
        while(!done)
        {
            ++numSamplesPerChainPerPred;
            if(numSamplesPerChainPerPred % samplesPerTest == 0)
            {
                secondsElapsed = (System.currentTimeMillis() - time) / 1000.0;
                System.out.println("Sample (per pred per chain) : " + numSamplesPerChainPerPred + ", Elapsed Time : " + Timer.time(secondsElapsed));
            }
            if(numSamplesPerChainPerPred == 4000)
            {
                System.out.println("hello...");
            }
            performGibbsStepForAllChains(burningIn);


            // Add current truth values to the convergence testers
            int g = 0;
            for (Integer gpId : state.groundMLN.indexToGroundPredMap.keySet())
            {
                double vals[] = new double[numChains];
                for (int c = 0; c < numChains; c++) {
                    vals[c] = truthValues[c].get(gpId);
                }
                if (burningIn)
                    burnConvergenceTests[g].appendNewValues(vals);
                else
                    gibbsConvergenceTests[g].appendNewValues(vals);
                g++;
            }

            if (numSamplesPerChainPerPred % samplesPerTest != 0) continue;

            if (burningIn)
            {
                // Use convergence criteria stated in "Probability and Statistics",
                // DeGroot and Schervish
                boolean burnConverged = false;

                if (testConvergence)
                    burnConverged =
                            GelmanConvergenceTest.checkConvergenceOfAll(burnConvergenceTests, state.groundMLN.indexToGroundPredMap.size(), true);
                if ( (numSamplesPerChainPerPred >= burnMinSteps && burnConverged)
                        || (burnMaxSteps >= 0 && numSamplesPerChainPerPred >= burnMaxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds))
                {
                    System.out.println("Done burning. " + numSamplesPerChainPerPred + " samples per pred per chain");
                    if (testConvergence)
                    {
                        System.out.println(" (" + (burnConverged? "converged":"didn't converge")
                                +" at total of " + numChains*numSamplesPerChainPerPred + " samples per pred)");
                    }
                    burningIn = false;
                    numSamplesPerChainPerPred = 0;
                }
            }

            else
            {  // Convergence test for actual gibbs sampling
                boolean gibbsConverged = false;

                if (testConvergence)
                    gibbsConverged = ConvergenceTest.checkConvergenceOfAtLeast(gibbsConvergenceTests, state.groundMLN.indexToGroundPredMap.size(), numSamplesPerChainPerPred, fracConverged, true);

                if (   (numSamplesPerChainPerPred >= minSteps && gibbsConverged)
                        || (maxSteps >= 0 && numSamplesPerChainPerPred >= maxSteps)
                        || (maxSeconds > 0 && secondsElapsed >= maxSeconds))
                {
                    System.out.println("Done Gibbs sampling. " + numSamplesPerChainPerPred
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

        if(calculateMarginal)
        {
            // update gndPreds probability
            for (Integer g : state.groundMLN.indexToGroundPredMap.keySet()) {
                int numVals = state.groundMLN.indexToGroundPredMap.get(g).numPossibleValues;
                for(int val = 0 ; val < numVals; val++)
                {
                    double marginal = 0.0;
                    for (int c = 0; c < numChains; c++) {
                        marginal += numValPerChainPerPred_.get(c).get(g).get(val)/numSamplesPerChainPerPred;
                    }
                    marginal /= numChains;
                    numValPerPred_.get(g).set(val,marginal);
                }
            }
        }
        if(trackFormulaTrueCnts)
        {
            int totalSamples = allFormulaTrueCntsPerChain.get(0).size();
            if(saveAllCounts)
            {
                int numCounts = formulaTrueCnts.length;
                int startSample = allFormulaTrueCnts.size();
                for (int sample = startSample; sample < totalSamples; sample++) {
                    allFormulaTrueCnts.add(new ArrayList<>(Collections.nCopies(numCounts,0.0)));
                    int size = allFormulaTrueCnts.size();
                    for (int formulaId = 0; formulaId < numCounts; formulaId++) {
                        double totalCount = 0.0;
                        for (int chainIdx = 0; chainIdx < numChains; chainIdx++) {
                            totalCount += allFormulaTrueCntsPerChain.get(chainIdx).get(sample).get(formulaId)/numChains;
                        }
                        allFormulaTrueCnts.get(size-1).set(formulaId, totalCount);
                    }
                }
            }

            Arrays.fill(formulaTrueCnts, 0.0);
            Arrays.fill(formulaTrueSqCnts, 0.0);

            for (int sampleNum = 0; sampleNum < allFormulaTrueCnts.size() ; sampleNum++) {
                int numWts = allFormulaTrueCnts.get(sampleNum).size();
                for (int formulaId = 0; formulaId < numWts; formulaId++) {
                    double trueCntsInASample = allFormulaTrueCnts.get(sampleNum).get(formulaId);
                    formulaTrueCnts[formulaId] += trueCntsInASample/totalSamples;
                    formulaTrueSqCnts[formulaId] += (trueCntsInASample * trueCntsInASample)/totalSamples;
                }
            }

        }
    }

    private void performGibbsStepForAllChains(boolean burningIn) {
        GibbsPerChain gibbsPerChain;
        boolean withThread = true;
        if(withThread)
        {
            Thread t[] = new Thread[numChains];
            Thread.currentThread().setPriority(10);
            // For each chain, for each node, generate the node's new truth value
            for (int c = 0; c < numChains; c++) {
                gibbsPerChain = new GibbsPerChain(c, burningIn, trackFormulaTrueCnts);
//                gibbsPerChain.run();
                t[c] = new Thread(gibbsPerChain);
                t[c].setPriority(1);
            }
            try {
                for (int c = 0; c < numChains; c++) {
                    t[c].start();
                }
                //System.out.println("Waiting for numChains threads to finish.");
                Thread.currentThread().setPriority(1);
                for (int c = 0; c < numChains; ++c) {
                    t[c].join();
                }
                Thread.currentThread().setPriority(10);
            }
            catch (InterruptedException e) {
                System.out.println("Main thread Interrupted");
            }
        }
        else
        {
            for (int c = 0; c < numChains; c++) {
                gibbsPerChain = new GibbsPerChain(c, burningIn, trackFormulaTrueCnts);
                gibbsPerChain.run();

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
        for(Integer gpId : state.groundMLN.indexToGroundPredMap.keySet()){

            int assignment = performGibbsStep(chainIdx,gpId);
            // If in actual gibbs sampling phase, update numValPerPred
            if(!burningIn)
                numValPerChainPerPred_.get(chainIdx).get(gpId).set(assignment, numValPerChainPerPred_.get(chainIdx).get(gpId).get(assignment)+1);
        }
    }

    /**
     * Performs one step of Gibbs sampling for predicate at gpId.
     * @param chainIdx Index of chain in which the Gibbs step is performed.
     * @param gpId PredicateId for which we try to flip
     * @return assignment by gibbs step of gpId
     */
    int performGibbsStep(int chainIdx, int gpId)
    {
//        if(gpIdsToBeChanged.get(chainIdx).contains(gpId)) {
//            updateWtsForGndPred(chainIdx, gpId);
//            gpIdsToBeChanged.get(chainIdx).remove(gpId);
//        }
        updateWtsForGndPred(chainIdx, gpId);
        int assignment = get_probabilistic_assignment(wtsPerPredPerVal.get(chainIdx).get(gpId));
        int prev_assignment = truthValues[chainIdx].get(gpId);
        truthValues[chainIdx].put(gpId, assignment);
        if(assignment != prev_assignment)
        {
            updateSatCounts(chainIdx, gpId, assignment, prev_assignment);
            //updateToBeChanged(chainIdx, gpId);
        }
        return assignment;
    }

//    private void updateToBeChanged(int chainIdx, int gpId) {
//        GroundPredicate gp = state.groundMLN.indexToGroundPredMap.get(gpId);
//        gpIdsToBeChanged.get(chainIdx).addAll(gp.markovBlanketSet);
//    }

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

    @Override
    public void resetCnts() {
        super.resetCnts();
        randomInitializeTruthVals();
        initializeNumSatValues();
        // Initialize convergence test
        int numGndPreds = state.groundMLN.indexToGroundPredMap.size();
        initConvergenceTests(gamma, epsilonError, numGndPreds, numChains);
    }
}