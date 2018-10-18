package org.utd.cs.mln.inference;

import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.GroundFormula;
import org.utd.cs.mln.alchemy.core.State;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by Happy on 5/28/18.
 */
public class MCSAT extends MCMC {
    // Number of samples after which print info
    private int samplesPerTest;
    // For each chain, stores indices of ground formulas which are required to be satisfied in SampleSat
    private List<Set<Integer>> toBeSatisfiedFormulasPerChain;
    // For each chain, stores

    /**
     * Constructor: User-set parameters are set.
     *
     * @param state
     * @param seed
     * @param trackFormulaTrueCnts
     * @param saveAllCounts
     * @param priorSoftEvidence
     * @param params
     * @param calculateMarginal
     * @see Inference#Inference(State, int, boolean, boolean, boolean)
     */
    public MCSAT(State state, int seed, boolean trackFormulaTrueCnts, boolean saveAllCounts, boolean priorSoftEvidence, MCMCParams params, boolean calculateMarginal) {
        super(state, seed, trackFormulaTrueCnts, saveAllCounts, priorSoftEvidence, params, calculateMarginal);
        toBeSatisfiedFormulasPerChain = new ArrayList<>();
        for (int i = 0; i < params.numChains; i++) {
            toBeSatisfiedFormulasPerChain.add(new HashSet<Integer>());
        }
        init();
    }

    @Override
    void init() {
        randomInitializeTruthVals();
    }

    @Override
    public void infer() {
        // Burn-in only if burnMaxSteps positive
        boolean burningIn = (burnMaxSteps > 0) ? true : false;
        // Sampling loop
        int numSamplesPerChain = 0;
        boolean done = false;
        long time = System.currentTimeMillis();
        double secondsElapsed = 0;
        while(!done)
        {
            ++numSamplesPerChain;
            if(numSamplesPerChain % samplesPerTest == 0)
            {
                secondsElapsed = (System.currentTimeMillis() - time) / 1000.0;
                System.out.println("Sample (per pred per chain) : " + numSamplesPerChain + ", Elapsed Time : " + Timer.time(secondsElapsed));
            }
            performMCSatStepForAllChains(burningIn);
        }
    }

    private void performMCSatStepForAllChains(boolean burningIn) {
        MCSatPerChain mcsatPerChain;
        boolean withThread = false;
        if(withThread)
        {
            Thread t[] = new Thread[numChains];
            Thread.currentThread().setPriority(10);
            // For each chain, for each node, generate the node's new truth value
            for (int c = 0; c < numChains; c++) {
                mcsatPerChain = new MCSatPerChain(c, burningIn, trackFormulaTrueCnts);
//                gibbsPerChain.run();
                t[c] = new Thread(mcsatPerChain);
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
                mcsatPerChain = new MCSatPerChain(c, burningIn, trackFormulaTrueCnts);
                mcsatPerChain.run();

            }
        }
    }

    private class MCSatPerChain implements Runnable{
        private final int chainIdx;
        private final boolean burningIn, trackFormulaTrueCnts;
        public MCSatPerChain(int c, boolean burningIn, boolean trackFormulaTrueCnts) {
            this.chainIdx = c;
            this.burningIn = burningIn;
            this.trackFormulaTrueCnts = trackFormulaTrueCnts;
        }

        @Override
        public void run() {
            MCSAT.this.performMCSatStep(chainIdx, burningIn);
            if(!burningIn && trackFormulaTrueCnts)
            {
                updateTrueCnts(chainIdx);
            }
        }
    }

    /**
     * Performs one step of MC-SAT sampling in one chain.
     *
     * @param chainIdx Index of chain in which the Gibbs step is performed.
     * @param burningIn If true, burning-in is occuring. Otherwise, false.
     */
    private void performMCSatStep(int chainIdx, boolean burningIn) {
        // For each ground formula, check if it is satisfied or not. It it is, add it to
        // toBeSatisfiedFormulasPerChain list with some prob p = (1-e^-w)
        int numGroundFormulas = state.groundMLN.groundFormulas.size();
        for (int i = 0; i < numGroundFormulas; i++) {
            GroundFormula gf = state.groundMLN.groundFormulas.get(i);
            // If this ground formula is already present in toBeSatisfiedFormulasPerChain, then
            // no need to check for satisfiability, just add it with prob p, because in this state,
            // all the formulas in toBeSatisfiedFormulasPerChain are already satisfied since this state
            // comes from last sample sat step.
            if(toBeSatisfiedFormulasPerChain.get(chainIdx).contains(i))
            {
                double p = (1-Math.exp(-Math.abs(gf.weight.getValue())));
                if(rand.nextDouble() < p)
                    toBeSatisfiedFormulasPerChain.get(chainIdx).add(i);
            }
            // Else check the satisfiability of this ground formula and if satisfied, add it
            // with prob p
            else if(gf.isSatisfied(truthValues[chainIdx]))
            {
                double p = (1-Math.exp(-Math.abs(gf.weight.getValue())));
                if(rand.nextDouble() < p)
                    toBeSatisfiedFormulasPerChain.get(chainIdx).add(i);
            }
            SampleSat ss = new SampleSat(state, toBeSatisfiedFormulasPerChain.get(chainIdx), truthValues[chainIdx]);
            ss.performSampleSatStep();
        }
    }
}
