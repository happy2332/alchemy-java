package org.utd.cs.mln.learning;

import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.inference.GibbsParams;
import org.utd.cs.mln.inference.GibbsSampler_v3;
import org.utd.cs.mln.inference.Inference;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 * Created by Happy on 4/11/18.
 */
public class DiscLearn extends WeightLearner {
    private List<Inference> inferences;
    private double[] gradient;
    private long time;

    // for each domain d, for each first order formula f, stores number of true groundings according to domain d.
    // formulaTrainCnts[d][f]. If priorSoftEvidence, last entry is trainCnt for lambda
    private double[][] formulaTrainCnts;

    // for each domain d, for each first order formula f, stores expected true groundings according to domain d.
    // formulaInferedCnts[d][f]. If priorSoftEvidence, last entry is for lambda
    private double[][] formulaInferedCnts;

    // for each domain d, for each first order formula f, stores number of total groundings according to domain d.
    // formulaNumGndings[d][f]
    private double[][] formulaNumGroundings;

    //max number of learning iterations
    private int numIter;

    /**
     * Constructor : Initializes the fields
     *
     * @param mlnsParam List of MLNs of size number of databases. All MLNs are same except for domains for predicates.
     * @param lArgs
     */
    public DiscLearn(List<MLN> mlnsParam, List<GroundMLN> groundMlnsParam, List<GroundMLN> groundMlnsEMParam,
                     List<Evidence> truthsParam, List<Evidence> truthsEMParam, LearnArgs lArgs) throws FileNotFoundException {
        super(mlnsParam, groundMlnsParam, groundMlnsEMParam, truthsParam, truthsEMParam, lArgs);
        init(groundMlnsParam, truthsParam, groundMlnsEMParam, truthsEMParam, lArgs);
    }

    private void init(List<GroundMLN> groundMlnsParam, List<Evidence> truthsParam, List<GroundMLN> groundMlnsEMParam, List<Evidence> truthsEMParam, LearnArgs lArgs)
    {
        //assign memory to all data members here
        inferences = new ArrayList<>();
        gradient = new double[numWts];
        formulaTrainCnts = new double[domain_cnt][numWts];
        formulaInferedCnts = new double[domain_cnt][numWts];
        formulaNumGroundings = new double[domain_cnt][numWts];

        // Initialize required data members here
        for (int i = 0; i < domain_cnt; i++) {
            GibbsParams gibbsparams = new GibbsParams();
            Inference inference = new GibbsSampler_v3(states.get(i),-1,true, priorSoftEvidence, gibbsparams);
            inferences.add(inference);
        }

        numIter = lArgs.dNumIter;
    }

    @Override
    public void learnWeights() {
        time = System.currentTimeMillis();
        if(!withEM)
            findFormulaTrainCnts();

        // Start learning weights
        for (int iter = 0; iter < numIter; iter++) {
            // At the start of every iteration, set MLN weight to weights learned till now. In case of first iteration,
            // this will be initial weights which are assigned in weights array.
            setMLNWeights();

            // Also set weight for every ground formula according to new MLN weights set
            setGroundMLNWeights();

            // Now find gradient
            findAndSetGradient();
        }

    }

    private void setGroundMLNWeights() {
        for (int i = 0; i < domain_cnt; i++) {
            State state = inferences.get(i).state;
            state.groundMLN.setGroundFormulaWtsToSumOfParentWts(state.mln);
        }
    }

    /**
     * This method computes gradient of log likelihood of data and fills in with gradient field. Note that we want
     * to minimize neg log likelihood, so gradient of that will be negative of this computed gradient.
     * g[i] = (n_i - E[n_i]) - reg[i], where n_i is number of true groundings of ith formula in training data,
     * E[n_i] is expected number of true groundings of ith formula, where expectation is over P(y|x),
     * and reg[i] comes from regularization, which is priorLambda[i]*(weights[i] - priorMeans[i])/priorStdDev[i].
     * n_i is stored in formulaTrainCnts[i]
     * To find E[n_i], run inference, and store results in formulaInferedCnts i.e. E[n_i] is formulaInferedCnts[i].
     */
    private void findAndSetGradient()
    {
        // For finding gradient, we need to do inference first so that we can calculate expected true counts.
        //infer(inferences, formulaInferedCnts);

        // for each domain, calculate gradient, and then add them up.
        for (int domainId = 0; domainId < domain_cnt; domainId++) {

        }
    }

//    /**
//     * This method runs inference for all domains in parallel, and then fills in formulaInferedCntsParam
//     * @param inferencesParam List of inference objects over which inference is to be done
//     * @param formulaInferedCntsParam Fill infered counts in this 2D array
//     */
//    private void infer(List<Inference> inferencesParam, double [][] formulaInferedCntsParam) {
//        InferPerDomain inferPerDomain;
//        Thread t[] = new Thread[domain_cnt];
//        MLN mln = inferencesParam.get(0).state.mln;
//        Thread.currentThread().setPriority(10);
//        for (int i = 0; i < domain_cnt; i++) {
//            State state = inferences.get(i).state;
//            if(priorSoftEvidence) {
//                state.groundMLN.setGroundFormulaWtsToParentWts(mln, weights[numFormulas]);
//            } else {
//                state.groundMLN.setGroundFormulaWtsToParentWts(mln);
//            }
//            inferences.get(i).updateWtsForNextGndPred(0);
//            ipmln = new InferPerMLN(inferences.get(i), burningIn, isInit);
//            t[i] = new Thread(ipmln);
//            t[i].setPriority(1);
//        }
//
//        try {
//            for (int i = 0; i < domain_cnt; ++i){
//                System.out.println("Doing inference for domain " + i);
//                t[i].start();
//            }
//            System.out.println("Waiting for threads to finish.");
//            Thread.currentThread().setPriority(1);
//            for (int i = 0; i < domain_cnt; ++i){
//                t[i].join();
//            }
//            Thread.currentThread().setPriority(10);
//        } catch (InterruptedException e) {
//            System.out.println("Main thread Interrupted");
//        }
//    }

    private void findFormulaTrainCnts() {
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
                    formulaNumGroundings[i][parentFormulaId.get(j)] += numCopies.get(j);
                }
                if(isFormulaSatisfied && parentFormulaId.isEmpty())
                {
                    formulaTrainCnts[i][numWts-1] += gf.originalWeight.getValue();
                }

            }
        }
    }
}
