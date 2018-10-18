package org.utd.cs.mln.learning;

import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.ConditionalLogLikelihood;
import org.utd.cs.mln.alchemy.util.OtherUtils;
import org.utd.cs.mln.inference.GibbsParams;
import org.utd.cs.mln.inference.GibbsSampler_v3;
import org.utd.cs.mln.inference.Inference;
import org.utd.cs.mln.inference.MCMC;

import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by Happy on 4/11/18.
 */
public class DiscLearn extends WeightLearner {
    private boolean isPll = false;
    private String method = "cg";
    public List<Inference> inferences, inferencesEM;
    private long time;

    // for each domain d, for each first order formula f, stores number of true groundings according to domain d.
    // formulaTrainCnts[d][f]. If priorSoftEvidence, last entry is trainCnt for lambda, which is sum(w_i), where
    // i indexes over all softevidence ground formulas which are satisfied in the training data.
    // This is used for finding gradient.
    public double[][] formulaTrainCnts;

//    // for each domain d, for each first order formula f, stores number of total groundings according to domain d.
//    // formulaNumGndings[d][f]
//    // This is used for setting lambda of prior. lambda for prior determines weightage of the prior. Higher the number of
//    // groundings of a formula,
//    private double[][] formulaNumGroundings;

    private boolean dldebug = false;

    /**
     * Constructor : Initializes the fields
     *
     * @param mlnsParam List of MLNs of size number of databases. All MLNs are same except for domains for predicates.
     * @param lArgs
     */
    public DiscLearn(List<MLN> mlnsParam, List<GroundMLN> groundMlnsParam, List<GroundMLN> groundMlnsEMParam,
                     List<Evidence> truthsParam, List<Evidence> truthsEMParam, LearnArgs lArgs, List<Map<Integer,Integer>> groundHiddenPredMapMToEStep) throws FileNotFoundException {
        super(mlnsParam, groundMlnsParam, groundMlnsEMParam, truthsParam, truthsEMParam, lArgs, groundHiddenPredMapMToEStep);
        if(lArgs.debug)
            dldebug = true;
        isPll = lArgs.pll;
        method = lArgs.method;
        init();
    }

    private void init()
    {
        //assign memory to all data members here
        inferences = new ArrayList<>();
        if(withEM)
            inferencesEM = new ArrayList<>();
        formulaTrainCnts = new double[domain_cnt][numWts];
        // Initialize required data members here
        for (int i = 0; i < domain_cnt; i++) {
            GibbsParams gibbsparams = new GibbsParams(learnArgs.gibbsParam);
            Inference inference = new GibbsSampler_v3(states.get(i),-1,true, true, priorSoftEvidence, gibbsparams, false);
            inferences.add(inference);
            if(withEM)
            {
                GibbsParams gibbsparamsEM = new GibbsParams(learnArgs.gibbsParam);
                gibbsparamsEM.maxSteps = learnArgs.numEMSamples;
                gibbsparamsEM.numChains = learnArgs.numEMSamples;
                gibbsparamsEM.samplesPerTest = learnArgs.numEMSamples;
                Inference inferenceEM = new GibbsSampler_v3(statesEM.get(i), -1, false, true, priorSoftEvidence, gibbsparamsEM, false);
                inferencesEM.add(inferenceEM);
            }
        }
    }

    @Override
    public void learnWeights() {
        time = System.currentTimeMillis();

        if(method.equals("cg"))
        {
            learnByCG();
        }
        else
        {
            //TODO : learnByLBFGS
        }
    }

    private void learnByCGOld() {
        CGParams cgParams = new CGParams(learnArgs.cgParam);
        ConjugateGradient cg = new ConjugateGradient(cgParams, numWts);
        cg.cgdebug = learnArgs.debug;
        if(!isPll)
        {
            ConditionalLogLikelihood cllLoss = new ConditionalLogLikelihood(inferences, inferencesEM, withEM, learnArgs.usePrior, learnArgs.agg);
            cllLoss.cllDebug = learnArgs.debug;
            if(!withEM)
            {
                findFormulaTrainCnts();
                cllLoss.formulaTrainCnts = formulaTrainCnts;
            }
            time = System.currentTimeMillis();
            // Start learning weights
            for (int iter = 1; iter <= cg.numIter; iter++) {
                System.out.println("ITER : " + iter);
                // At the start of every iteration, set MLN weight to weights learned till now.
                // We need this because later we want to set weights of ground formulas according to weights learned
                // of first order formulas, which are set according to MLN weights. In case of first iteration,
                // this will be initial weights which are assigned in weights array.
                setMLNWeights();

                // Also set weight for every ground formula according to new MLN weights set
                setGroundMLNWeights();

                // Now find gradient and fill gradient field.
                double gradient[] = cllLoss.getGradient(weights);
                int status = cg.updateWts(weights, gradient, inferences, inferencesEM, iter, withEM, cllLoss);
                if(status == -1)
                    break;
                if (cg.backtracked) {
                    iter--;

                }
                else{
                    // If we are backtracking, then we DON'T want to reset formulaCounts as we want more samples to add on.
                    for (int i = 0; i < domain_cnt; i++) {
                        inferences.get(i).resetCnts();
                        if(withEM)
                            inferencesEM.get(i).resetCnts();
                    }
                }
                System.out.println("ELAPSED TIME : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
            }
        }
    }

    private void learnByCG() {
        CGParams cgParams = new CGParams(learnArgs.cgParam);
        ConjugateGradient cg = new ConjugateGradient(cgParams, numWts);
        cg.cgdebug = learnArgs.debug;
        if (!isPll) {
            ConditionalLogLikelihood cllLoss = new ConditionalLogLikelihood(inferences, inferencesEM, withEM, learnArgs.usePrior, learnArgs.agg);
            cllLoss.cllDebug = learnArgs.debug;
            if (!withEM) {
                findFormulaTrainCnts();
                cllLoss.formulaTrainCnts = formulaTrainCnts;
            }
            time = System.currentTimeMillis();
            // Start learning weights
            for (int iter = 1; iter <= cg.numIter; iter++) {
                System.out.println("ITER : " + iter);
                // At the start of every iteration, set MLN weight to weights learned till now.
                // We need this because later we want to set weights of ground formulas according to weights learned
                // of first order formulas, which are set according to MLN weights. In case of first iteration,
                // this will be initial weights which are assigned in weights array.
                setMLNWeights();

                // Also set weight for every ground formula according to new MLN weights set
                setGroundMLNWeights();
                double gradient[] = new double[numWts];
                if (withEM) {
                    cllLoss.infer(inferencesEM);
                    for (int i = 0; i < learnArgs.numEMSamples; i++) {
                        cllLoss.formulaTrainCnts = findFormulaTrainCntsWithEM(i);
                        // Now find gradient and fill gradient field.
                        double gradientPerSample[] = cllLoss.getGradient(weights);
                        for (int j = 0; j < numWts; j++) {
                            gradient[j] += gradientPerSample[j];
                        }
                    }
                    for (int j = 0; j < numWts; j++) {
                        gradient[j] /= learnArgs.numEMSamples;
                    }
                }
                else
                {
                    // Now find gradient and fill gradient field.
                    gradient = cllLoss.getGradient(weights);
                }

                int status = cg.updateWts(weights, gradient, inferences, inferencesEM, iter, withEM, cllLoss);
                if (status == -1)
                    break;
                if (cg.backtracked) {
                    iter--;

                }
                else{
                    // If we are backtracking, then we DON'T want to reset formulaCounts as we want more samples to add on.
                    for (int i = 0; i < domain_cnt; i++) {
                        inferences.get(i).resetCnts();
//                        if(withEM)
//                            inferencesEM.get(i).resetCnts();
                    }
                }
                System.out.println("ELAPSED TIME : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
            }
        }
    }

    private double[][] findFormulaTrainCntsWithEM(int sampleIndex) {
        double formulaTrainCnts[][] = new double[domain_cnt][numWts];
        for (int i = 0; i < domain_cnt; i++) {
            GroundMLN gm = inferences.get(i).state.groundMLN;
            for(GroundFormula gf : gm.groundFormulas)
            {
                List<Integer> parentFormulaId = gf.parentFormulaId;
                List<Integer> numCopies = gf.numCopies;
                boolean isFormulaSatisfied = true;
                for(GroundClause gc : gf.groundClauses)
                {
                    boolean isClauseSatisfied = false;
                    for(int gpId : gc.groundPredIndices)
                    {
                        int trueVal = 0;
                        if(groundHiddenPredMapMToEStep.get(i).containsKey(gpId))
                        {
                            if(groundHiddenPredMapMToEStep.get(i).get(gpId) != -1)
                            {
                                int gpIdEStep = groundHiddenPredMapMToEStep.get(i).get(gpId);
                                trueVal = ((MCMC)inferencesEM.get(i)).truthValues[sampleIndex].get(gpIdEStep);
                            }
                            else
                            {
                                int numPossibleVals = gm.indexToGroundPredMap.get(gpId).numPossibleValues;
                                trueVal = OtherUtils.getUniformAssignment(numPossibleVals);
                            }
                        }
                        else
                        {
                            trueVal = inferences.get(i).state.truthVals.get(gpId);
                        }
                        BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                        isClauseSatisfied |= b.get(trueVal);
                        if(isClauseSatisfied) {
                            break;
                        }
                    }
                    isFormulaSatisfied &= isClauseSatisfied;
                    if(!isFormulaSatisfied) {
                        break;
                    }
                }

                for (int j = 0; j < parentFormulaId.size(); j++) {
                    if(isFormulaSatisfied) {
                        formulaTrainCnts[i][parentFormulaId.get(j)] += numCopies.get(j);
                    }
                    //formulaNumGroundings[i][parentFormulaId.get(j)] += numCopies.get(j);
                }
                if(isFormulaSatisfied && parentFormulaId.isEmpty())
                {
                    formulaTrainCnts[i][numWts-1] += gf.originalWeight.getValue();
                }
            }
        }
        return formulaTrainCnts;
    }

    /**
     * This method fills in {@link DiscLearn#formulaTrainCnts}.
     */

    public void findFormulaTrainCnts() {
        for (int i = 0; i < domain_cnt; i++) {
            GroundMLN gm = inferences.get(i).state.groundMLN;
            for(GroundFormula gf : gm.groundFormulas)
            {
                List<Integer> parentFormulaId = gf.parentFormulaId;
                List<Integer> numCopies = gf.numCopies;
                boolean isFormulaSatisfied = true;
                for(GroundClause gc : gf.groundClauses)
                {
                    boolean isClauseSatisfied = false;
                    for(int gpId : gc.groundPredIndices)
                    {
                        int trueVal = inferences.get(i).state.truthVals.get(gpId);
                        BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                        isClauseSatisfied |= b.get(trueVal);
                        if(isClauseSatisfied) {
                            break;
                        }
                    }
                    isFormulaSatisfied &= isClauseSatisfied;
                    if(!isFormulaSatisfied) {
                        break;
                    }
                }

                for (int j = 0; j < parentFormulaId.size(); j++) {
                    if(isFormulaSatisfied) {
                        formulaTrainCnts[i][parentFormulaId.get(j)] += numCopies.get(j);
                    }
                    //formulaNumGroundings[i][parentFormulaId.get(j)] += numCopies.get(j);
                }
                if(isFormulaSatisfied && parentFormulaId.isEmpty())
                {
                    formulaTrainCnts[i][numWts-1] += gf.originalWeight.getValue();
                }
            }
        }
    }


    private void setGroundMLNWeights() {
        for (int i = 0; i < domain_cnt; i++) {
            State state = inferences.get(i).state;
            State stateEM = null;
            if(withEM)
            {
                stateEM = inferencesEM.get(i).state;
            }
            if(learnArgs.agg)
            {
                //state.groundMLN.setNumConnections();
                state.groundMLN.setEffWts(state.mln);
                if(withEM)
                    stateEM.groundMLN.setEffWts(stateEM.mln);
            }
            else
            {
                state.groundMLN.setGroundFormulaWtsToSumOfParentWts(state.mln);
                if(withEM)
                    stateEM.groundMLN.setGroundFormulaWtsToSumOfParentWts(stateEM.mln);
            }
        }
    }

}
