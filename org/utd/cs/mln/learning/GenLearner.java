package org.utd.cs.mln.learning;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.Evidence;
import org.utd.cs.mln.alchemy.core.GroundMLN;
import org.utd.cs.mln.alchemy.core.MLN;
import org.utd.cs.mln.alchemy.util.Loss;
import org.utd.cs.mln.alchemy.util.PseudoLogLikelihood;
import util.LBFGS;
import util.U;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * This class implements generative learning of MLN by maximizing pseudologlikelihood of data.
 * The algorithm used is LBFGS, for which we use @see <a href="https://github.com/brendano/myutil">this library</a>.
 * @author Happy
 * @since 27/03/18
 */
public class GenLearner extends WeightLearner{

    private boolean genLearnDebug = false;
    private long time;
    private boolean isPll = false;
    private String method = "cg";
    Loss loss;

    /**
     * Constructor : Initializes fields
     * @param mlnsParam : List of MLNs of size number of databases. All MLNs are same except for domains for predicates
     * @param groundMlnsParam : List of groundMLNs of size number of databases.
     * @param groundMlnsEMParam : List of groundMLNs of size number of databases (for EM)
     * @param truthsParam : For each database, Truth values of every ground predicate in corresponding groundMLN.
     * @param truthsEMParam : For each database, Truth values of every ground predicate in corresponding groundMLN (for EM)
     * @param lArgs : Learning Arguments
     * @throws FileNotFoundException
     */
    public GenLearner(List<MLN> mlnsParam, List<GroundMLN> groundMlnsParam, List<GroundMLN> groundMlnsEMParam,
                      List<Evidence> truthsParam, List<Evidence> truthsEMParam, LearnArgs lArgs, List<Map<Integer,Integer>> groundHiddenPredMapMToEStep) throws FileNotFoundException {
        super(mlnsParam, groundMlnsParam, groundMlnsEMParam, truthsParam, truthsEMParam, lArgs, groundHiddenPredMapMToEStep);
        if(lArgs.debug)
            genLearnDebug = true;
        isPll = lArgs.pll;
        method = lArgs.method;
        //init(groundMlnsParam, truthsParam, lArgs);
    }

    public void learnWeights() {
        time = System.currentTimeMillis();
        setMLNWeights();
        if (isPll) {
            loss = new PseudoLogLikelihood(states, learnArgs.usePrior);
        } else {
            //TODO : get loss for likelihood
        }
        if (method.equals("lbfgs")) {
            learnByLBFGS();
        }

//        if(withEM) {
//            inferEM(true, true);
//            subtypeTrueVal = getStTrueVal();
//            findFormulaTrainCnts();
//            updateMLN(true, 0);
//        }
//        else
        {
//            if(!withEM){
            //updateMLN(true);

        }
//        prior = setPrior();
//        System.out.println("Prior: " + Arrays.toString(prior));
    }

    void learnByLBFGS()
    {
        LBFGS.Function f = new LBFGS.Function() {
            @Override
            public double evaluate(double[] beta, double[] g, int n, double step) {
                double logLoss = 0;
                Arrays.fill(g, 0);
                if(learnArgs.agg)
                {
                    for (int i = 0; i < states.get(0).mln.formulas.size(); i++) {
                        weights[i] /= LearnTest.firstOrderNumConnections.get(i);
                    }
                }
                setMLNWeights();

//                if(withEM){
//                    int numIter = inferencesEM.get(0).numIter;
//                    inferEM(true, true);
//                    for (int iter = 0; iter < numIter; ++iter){
//                        updateMLN(false, iter);
//
//                        for (int domainID = 0; domainID < domain_cnt; ++domainID){
//                            inferences.get(domainID).findFtcPerPredicate();
////                                inferences.get(domainID).debugFTC(10);
//                        }
//
//                        logLoss -= getPseudoLogLikelihood(prior) / numIter;
//                        double []tempGrad = getPseudoGradient(prior);
//                        for (int i = 0; i < g.length; i++) {
//                            g[i] -= tempGrad[i] / numIter;
//                        }
//                    }
//                }
//                else
                {
                    //updateMLN(false);

                    logLoss += loss.getLossValue();
                    double []tempGrad = loss.getGradient();
                    for (int i = 0; i < g.length; i++) {
                        g[i] += tempGrad[i];
                    }
                }
                if(learnArgs.agg)
                {
                    for (int i = 0; i < states.get(0).mln.formulas.size(); i++) {
                        weights[i] *= LearnTest.firstOrderNumConnections.get(i);
                    }
                }

//                System.out.println("Weights: "+Arrays.toString(beta));
//                System.out.println("Log Loss (-Log Likelihood): "+logLoss + "n : "+n);
//                System.out.println("Gradient: "+Arrays.toString(g));
                return logLoss;
            }
        };
        LBFGS.Params p = new LBFGS.Params();
        p.m = 5;
        p.epsilon = 1E-5D;
        p.max_iterations = 100;
        System.out.println("epsion : " + p.epsilon);
//            p.past = 5;
//            p.delta = 1.0E-7D;

        LBFGS.ProgressCallback cb = new LBFGS.ProgressCallback() {
            @Override
            public int apply(double[] x, double[] g, double fx, double xnorm,
                             double gnorm, double step, int n, int k, LBFGS.Status ls) {
                if(k % 1 == 0) {
                    U.pf("ITER %d obj=%g |proj g|=%g\n", k, fx, gnorm);
                    if(genLearnDebug)
                    {
                        System.out.println("Weights: "+Arrays.toString(x));
                        System.out.println("Gradient: "+Arrays.toString(g));
                    }

                    System.out.println("Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
                }
                return 0;
            }
        };
//            weights = new double[0];
        LBFGS.Result r = LBFGS.lbfgs(weights, f, cb, p);
        U.p(r);
    }
}

