package org.utd.cs.mln.learning;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.Parser;
import org.utd.cs.mln.inference.GibbsSampler_v2;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

import util.Arr;
import util.LBFGS;
import util.LBFGS.Function;
import util.U;

/**
 * Created by Happy on 3/11/17.
 */
public class DiscLearner {

    private double lambdaSoftEvidence, oldLambdaSoftEvidence = 1.0;
    private boolean priorSoftEvidence = false;
    private List<GibbsSampler_v2> inferences = new ArrayList<>();
    private List<GibbsSampler_v2> inferencesEM = new ArrayList<>();
    private double[] lambdaTrainCnts;
    private double[] weights, oldWeights, averageWeights, gradient, old_gradient, d, oldd,
            delta_pred, priorMeans, priorStdDevs;
    private double[][] formulaTrainCnts, formulaNumGroundings;
    private int num_iter, domain_cnt, backtrackCount, maxBacktracks;
    private double cg_lambda, cg_max_lambda, alpha, min_ll_change, lambda_grad;
    private boolean withEM = false, backtracked = false, preConditionCG, dldebug=true, usePrior=false;
    private LearnTest.Method method;
    private int numFormulas, numWts;
    private Map<Integer, Integer>[] subtypeMaps, subtypeTrueVal;        //subtypeTrueVal: for prior
    private double[] prior;
    double time;

    public DiscLearner(List<GibbsSampler_v2> inferences, List<GibbsSampler_v2> inferencesEM, Map<Integer, Integer>[] subtypeMaps, int num_iter, double lambda, double min_ll_change, double max_lambda,
                       boolean withEM, boolean predConditionCG, boolean usePrior, boolean priorSoftEvidence, double lambdaSoftEvidence, LearnTest.Method method)
    {
        this.inferences = inferences;
        this.inferencesEM = inferencesEM;
        this.num_iter = num_iter;
        this.backtrackCount = 0;
        this.maxBacktracks = 1000;
        this.cg_lambda = lambda;
        this.min_ll_change = min_ll_change;
        this.cg_max_lambda = max_lambda;
        this.alpha = 1;
        this.withEM = withEM;
        this.usePrior = usePrior;
        this.preConditionCG = predConditionCG;
        this.domain_cnt = inferences.size();
        this.priorSoftEvidence = priorSoftEvidence;
        this.numFormulas = inferences.get(0).mln.formulas.size();
        this.subtypeMaps = subtypeMaps;
        this.method = method;
//        this.likelihood = Likelihood.LOG;

        if(priorSoftEvidence){
            this.lambdaSoftEvidence = lambdaSoftEvidence;
            this.numWts = numFormulas + 1;
        }
        else
            this.numWts = numFormulas;

        weights = new double[numWts];
        oldWeights = new double[numWts];
        averageWeights = new double[numWts];
        gradient = new double[numWts];
        old_gradient = new double[numWts];
        d = new double[numWts];
        oldd = new double[numWts];
        lambdaTrainCnts = new double[domain_cnt];
        delta_pred = null;
        formulaTrainCnts = new double[domain_cnt][numFormulas];
        formulaNumGroundings = new double[domain_cnt][numFormulas];
        priorMeans = new double[numWts];
        priorStdDevs = new double[numWts];
        prior = new double[numWts];
        Arrays.fill(priorStdDevs, 4);
        Arrays.fill(priorMeans, 0);
        Arrays.fill(prior, 0);
    }

    public double[] setPrior() {
        double[] prior = new double[numWts];
        for (int i = 0; i < numWts; ++i){
            int temp = 0;
            for (int j = 0; j < domain_cnt; ++j){
                if(i == numFormulas){
//                    temp += inferences.get(j).state.groundMLN.sePredCount;
                } else {
                    temp += formulaTrainCnts[j][i];
//                    temp += formulaNumGroundings[j][i];
                }
            }
//            temp /= domain_cnt;
//            temp = 1;
            prior[i] = temp / (priorStdDevs[i] * priorStdDevs[i]);
        }
        return prior;
    }

    public double[] learnWeights() {
        initWeights();
        if(!withEM){
            findFormulaTrainCnts();
            if(priorSoftEvidence)
                findFormulaTrainCountsLambda();
        }

        for (int i = 0; i < domain_cnt; i++) {
            inferences.get(i).saveAllCounts(true);
            if(withEM)
                inferencesEM.get(i).saveAllCounts(true);
        }

        time = System.currentTimeMillis();
        boolean burningIn = true, isInit = true;

//        System.out.println("Prior: " + Arrays.toString(prior));
        System.out.println("Variance: " + priorStdDevs[0]);

        if(method == LearnTest.Method.LBFGS){
            setMLNWeights();
            if(withEM) {
                inferEM(true, true);
                subtypeTrueVal = getStTrueVal();
                findFormulaTrainCnts();
                updateMLN(true, 0);
            } else{
//            if(!withEM){
                updateMLN(true);
                for (int domainID = 0; domainID < domain_cnt; ++domainID){
                    inferences.get(domainID).findFtcPerPredicate();
                }
            }
            prior = setPrior();
            System.out.println("Prior: " + Arrays.toString(prior));

            Function f = new Function() {
                @Override
                public double evaluate(double[] beta, double[] g, int n, double step) {
                    double logLoss = 0;
                    Arrays.fill(g, 0);
                    setMLNWeights();

                    if(withEM){
                        int numIter = inferencesEM.get(0).numIter;
                        inferEM(true, true);
                        for (int iter = 0; iter < numIter; ++iter){
                            updateMLN(false, iter);

                            for (int domainID = 0; domainID < domain_cnt; ++domainID){
                                inferences.get(domainID).findFtcPerPredicate();
//                                inferences.get(domainID).debugFTC(10);
                            }

                            logLoss -= getPseudoLogLikelihood(prior) / numIter;
                            double []tempGrad = getPseudoGradient(prior);
                            for (int i = 0; i < g.length; i++) {
                                g[i] -= tempGrad[i] / numIter;
                            }
                        }
                    } else {
                        updateMLN(false);
                        logLoss -= getPseudoLogLikelihood(prior);
                        double []tempGrad = getPseudoGradient(prior);
                        for (int i = 0; i < g.length; i++) {
                            g[i] -= tempGrad[i];
                        }
                    }

                    System.out.println("Weights: "+Arrays.toString(beta));
                    System.out.println("Log Loss (-Log Likelihood): "+logLoss);
                    System.out.println("Gradient: "+Arrays.toString(g));
                    return logLoss;
                }
            };
            LBFGS.Params p = new LBFGS.Params();
//            p.m = 10;
//            p.epsilon = 1.0E-7D;
//            p.past = 5;
//            p.delta = 1.0E-7D;

            LBFGS.ProgressCallback cb = new LBFGS.ProgressCallback() {
                @Override
                public int apply(double[] x, double[] g, double fx, double xnorm,
                                 double gnorm, double step, int n, int k, LBFGS.Status ls) {
                    U.pf("ITER %d obj=%g sol=%.6g\n", k, fx, x[0]);
                    System.out.println("Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
                    return 0;
                }
            };
//            weights = new double[0];
            LBFGS.Result r = LBFGS.lbfgs(weights, f, cb, p);
            U.p(r);
        } else {
            prior = setPrior();
            System.out.println("Prior: " + Arrays.toString(prior));
            for(int iter = 1 ; iter <= num_iter ; iter++) {
                System.out.println("iter : " + iter + ", Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
                setMLNWeights();
                System.out.println("Running inference...");
                if(backtracked) {
                    burningIn = true;
                    isInit = true;
                }
                infer(burningIn, isInit);
                if(withEM)
                    inferEM(burningIn, isInit);
                burningIn = false;
                isInit = false;
                System.out.println("Done Inference");
                System.out.println("Getting gradient...");

                findGradient(prior);

                if(dldebug)
                    System.out.println("gradient = " + Arrays.toString(gradient));
                if(method == LearnTest.Method.CG) {
                    int status = updateWtsByCG(iter, prior);
                    if(status == -1)
                        break;
                    if (backtracked)
                        iter--;
                } else if(method == LearnTest.Method.VP){
                    updateWtsByVP(iter);
                }
            }
        }

        System.out.println("Learning done...");
        System.out.println("Final weights : " + Arrays.toString(weights));
        return weights;
    }

    private Map<Integer,Integer>[] getStTrueVal() {
        Map<Integer, Integer>[] subtypeTrueVal = new HashMap[domain_cnt];
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            subtypeTrueVal[domainID] = new HashMap<>();
            for (Integer predID : subtypeMaps[domainID].keySet()){
                int predID_EM = subtypeMaps[domainID].get(predID);
                Map<Integer, Integer> occurrence = new HashMap<>();
                for (int iter = 0; iter < inferencesEM.get(0).numIter; ++iter){
                    int predVal = inferencesEM.get(domainID).mlnStatePerIter[iter].get(predID_EM);
                    if(occurrence.containsKey(predVal)){
                        occurrence.put(predVal, occurrence.get(predVal) + 1);
                    } else {
                        occurrence.put(predVal, 1);
                    }
                }
                int stVal = 0;
                for (Integer key : occurrence.keySet()){
                    if(!occurrence.containsKey(stVal) || occurrence.get(stVal) < occurrence.get(key)){
                        stVal = key;
                    }
                }
                subtypeTrueVal[domainID].put(predID, stVal);
            }
        }
        return subtypeTrueVal;
    }

    //    without EM
    private void updateMLN(boolean init) {
        MLN mln = inferences.get(0).mln;
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            State state = inferences.get(domainID).state;
            if(priorSoftEvidence) {
                state.groundMLN.setGroundFormulaWtsToParentWtsSoftEvidence(mln, weights[numFormulas]);
            } else {
                state.groundMLN.setGroundFormulaWtsToParentWts(mln);
            }

            inferences.get(domainID).updateMLN(init);
        }
    }

//    with EM
    private void updateMLN(boolean init, int iter) {
        MLN mln = inferences.get(0).mln;
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            State state = inferences.get(domainID).state;
            if(priorSoftEvidence) {
                state.groundMLN.setGroundFormulaWtsToParentWtsSoftEvidence(mln, weights[numFormulas]);
            } else {
                state.groundMLN.setGroundFormulaWtsToParentWts(mln);
            }

            inferences.get(domainID).updateMLN(init, inferencesEM.get(domainID).mlnStatePerIter[iter], subtypeMaps[domainID]);
        }
    }

    private double getPseudoLogLikelihood(double[] prior) {
        double likelihood = 0;
        for (int domainIndex = 0; domainIndex < domain_cnt; ++domainIndex){
            likelihood += inferences.get(domainIndex).getPseudoLogLikelihood();
        }

        for (int formulaID = 0; formulaID < numWts; ++formulaID){
            double temp = weights[formulaID] - priorMeans[formulaID];
            likelihood -= temp * temp * prior[formulaID] / 2;
        }
        return likelihood;
    }

    private double[] getPseudoGradient(double[] prior) {
        double[] gradientForDomain;
        Arrays.fill(gradient,0.0);
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            gradientForDomain = inferences.get(domainID).getPseudoGradient(priorSoftEvidence);
            for (int formulaID = 0; formulaID < numWts; ++formulaID){
                gradient[formulaID] += gradientForDomain[formulaID];
            }
        }

        for (int formulaID = 0; formulaID < numWts; ++formulaID){
            gradient[formulaID] -= (weights[formulaID] - priorMeans[formulaID]) * prior[formulaID];
        }
        return gradient;
    }

    private void updateWtsByVP(int iter) {
        for (int i = 0; i < numFormulas; ++i) {
            double avgNi = 0;
            for (int j = 0; j < domain_cnt; ++j) {
                avgNi += formulaTrainCnts[j][i];
            }
            weights[i] -= ((1 / avgNi) * gradient[i]);
            averageWeights[i] = ((iter - 1) * averageWeights[i] + weights[i]) / iter;
        }
        for (int i = 0; i < domain_cnt; ++i){
            inferences.get(i).resetCnts();
        }
        System.out.println("weights = " + Arrays.toString(weights));
        System.out.println("Avg weights = " + Arrays.toString(averageWeights));
    }

    private void initWeights() {
//        MLN mln = inferences.get(0).mln;
//        if(usePrior){
//            for (int i = 0; i < numFormulas; i++) {
//                priorMeans[i] = mln.formulas.get(i).weight.getValue();
//            }
////            Random rand = new Random();
////            for (int i = 0; i < mln.formulas.size(); i++) {
////                double p = rand.nextGaussian()*priorStdDevs[i]+priorMeans[i];
////                weights[i] = p;
////            }
//        }
//        for (int i = 0; i < numFormulas; ++i) {
//            weights[i] = priorMeans[i];
//            averageWeights[i] = priorMeans[i];
//        }

        Arrays.fill(weights, 0);
        Arrays.fill(averageWeights, 0);

        if(priorSoftEvidence) {
            weights[numFormulas] = lambdaSoftEvidence;
            averageWeights[numFormulas] = lambdaSoftEvidence;
        }

        System.out.println("Initial Weights: " + Arrays.toString(weights));
    }

    private int updateWtsByCG(int iter, double[] prior) {
        double realdist = 1.0;
        double preddist = 1.0;
        if (iter > 1 && delta_pred != null && !backtracked) {
            double []dist = new double[numWts];
            for (int i = 0; i < numWts; i++) {
                dist[i] = weights[i] - oldWeights[i];
            }
            double []avgPred = new double[numWts];
            for (int i = 0; i < numWts; i++) {
                avgPred[i] = old_gradient[i] + delta_pred[i] / 2.0;
            }
            preddist = dotprod(avgPred, dist, numWts);

            // Real change is lower bound on actual change
            realdist = dotprod(gradient, dist, numWts);
            System.out.println("pred*dist = " + preddist);
            System.out.println("real*dist = " + realdist);
        }
        if(iter > 1) {
            double delta = realdist / preddist;

            if (!backtracked && preddist == 0)
                cg_lambda /= 4;

            if (!backtracked && preddist != 0 && delta > 0.75)
                cg_lambda /= 2;
            // if (delta < 0.25)   // Update schedule from (Fletcher, 1987)
            if (delta < 0.0)       // Gentler update schedule, to handle noise
            {
                if (cg_lambda * 4 > cg_max_lambda)
                    cg_lambda = cg_max_lambda;
                else
                    cg_lambda *= 4;
            }
            if (delta < 0.0 && backtrackCount < maxBacktracks) {
                System.out.println("Backtracking...");
                for (int i = 0; i < numWts; i++)
                    weights[i] = oldWeights[i];

                for (int i = 0; i < domain_cnt; i++)
                {
                    inferences.get(i).restoreCnts();
                    if(withEM)
                        inferencesEM.get(i).restoreCnts();
                }

                backtracked = true;
                backtrackCount++;
            } else {
                backtracked = false;
                backtrackCount = 0;
            }
        }

        if (!backtracked)
        {
            for (int i = 0; i < domain_cnt; i++)
            {
                inferences.get(i).saveCnts();
                if(withEM)
                    inferencesEM.get(i).saveCnts();
            }

            double preCond[] = new double[numWts];
            Arrays.fill(preCond, 1.0);
            if(preConditionCG)
            {
                double variance[] = getVariance(prior);
                if(withEM)
                {
                    double varianceEM[] = getVarianceEM(prior);
                    for (int i = 0; i < variance.length; i++) {
                        variance[i] -= varianceEM[i];
                    }
                }

                if(dldebug)
                    System.out.println("variance : " + Arrays.toString(variance));
                for (int formulaNum = 0; formulaNum < numWts; formulaNum++)
                    preCond[formulaNum] = 1.0/(variance[formulaNum] + 0.00001);
            }
            double beta = 0.0;

            // Compute beta using Polak-Ribiere form:
            //   beta = g_j+1 (g_j+1 - g_j) / (g_j g_j)
            // Preconditioned:
            //   beta = g_j+1 M-1 (g_j+1 - g_j) / (g_j M-1 g_j)
            double beta_num = 0.0;
            double beta_denom = 0.0;
            if (iter > 1)
            {
                for (int i = 0; i < numWts; i++)
                {
                    beta_num   += gradient[i] * preCond[i] * (gradient[i] - old_gradient[i]);
                    beta_denom += old_gradient[i] * preCond[i] * old_gradient[i];
                }
                beta = beta_num/(beta_denom + 0.00001);
            }
            else
                beta = 0.0;
            
            if(dldebug)
                System.out.println("beta = " + beta);

            // Compute new direction
            for (int w = 0; w < numWts; w++) {
                d[w] = -preCond[w] * gradient[w] + beta * oldd[w];
            }
            System.out.println("Direction D = " + Arrays.toString(d));

            double Hd[] = getHessianVectorProduct(d, prior);
            if(withEM)
            {
                double []HdEM = getHessianVectorProductEM(d, prior);
                for (int i = 0; i < Hd.length; i++) {
                    Hd[i] -= HdEM[i];
                }
            }
            alpha = computeQuadraticStepLength(Hd);
            if (alpha < 0.0)
            {
                for (int w = 0; w < numWts; w++)
                    d[w] = -preCond[w]*gradient[w];
                Hd = getHessianVectorProduct(d, prior);
                if(withEM)
                {
                    double []HdEM = getHessianVectorProductEM(d, prior);
                    for (int i = 0; i < Hd.length; i++) {
                        Hd[i] -= HdEM[i];
                    }
                }
                alpha = computeQuadraticStepLength(Hd);
            }
        }

        if (!backtracked && alpha < 0.0)
        {
            // If alpha is negative, then either the direction or the
            // Hessian is in error.  We call this a backtrack so that
            // we can gather more samples while keeping the old samples.
            backtracked = true;
        }
        if (!backtracked) {
            // Compute total weight change
            double wchange[] = new double[numWts];
            for (int w = 0; w < numWts; w++) {
                //wchange[w] = d[w] * alpha + (weights[w] - oldWeights[w]) * momentum;
                // above line was present in alchemy, but momentum is 0 always for disc learning
                wchange[w] = d[w] * alpha;
            }

            // Convergence criteria for 2nd order methods:
            // Stop when the maximum predicted improvement in log likelihood
            // is very small.
            double maxchange = -dotprod(gradient, wchange, numWts);
            System.out.println("Maximum Estimated Improvement = " + maxchange);
            if ((method == LearnTest.Method.CG) && maxchange < min_ll_change) {
                System.out.println("Upper bound is less than " + min_ll_change + ", halting learning.");
                return -1;
            }

            // Save weights, gradient, and direction and adjust the weights
            for (int w = 0; w < numWts; w++) {
                oldWeights[w] = weights[w];
                oldd[w] = d[w];
                old_gradient[w] = gradient[w];

                if(w < numFormulas){
                    if(!inferences.get(0).mln.formulas.get(w).tautology)
                        weights[w] += wchange[w];
                } else {
                    weights[w] += wchange[w];
                }
                averageWeights[w] = ((iter - 1) * averageWeights[w] + weights[w]) / iter;
            }
            delta_pred = getHessianVectorProduct(wchange, prior);
            //TODO : delta pred for EM ?
            for (int i = 0; i < domain_cnt; i++) {
                inferences.get(i).resetCnts();
                if(withEM)
                    inferencesEM.get(i).resetCnts();
            }
            System.out.println("weights = " + Arrays.toString(weights));
            System.out.println("Avg weights = " + Arrays.toString(averageWeights));
        }
        return 0;
    }

    private double computeQuadraticStepLength(double[] hd) {
        int numWeights = d.length;
        // Compute step length using trust region approach
        double dHd = dotprod(d, hd, numWeights);
        double dd = dotprod(d, d, numWeights);
        double dg = dotprod(gradient, d, numWeights);
        double alpha_denom = (dHd + cg_lambda * dd);
        if(alpha_denom == 0.0)
            return 0.0;
        double alpha = -dg/alpha_denom;

        if(dldebug)
        {
            System.out.println("dHd = " + dHd);
            System.out.println("dd = " + dd);
            System.out.println("dg = " + dg);
            System.out.println("alpha = " + alpha);
        }

        // Because the problem is convex, the Hessian should always
        // be positive definite, and so alpha should always be non-negative.
        // Since our Hessian is approximate, we may end up with a negative
        // alpha.  In these cases, we used to make alpha zero (i.e., take a
        // step of size zero), but now we return the negative alpha and
        // let the caller deal with it.
        if (alpha < 0.0)
        {
            System.out.println("Alpha < 0!  Bad direction or Hessian.");
        }

        return alpha;
    }

    private double[] getHessianVectorProduct(double[] d, double[] prior) {
        double Hd[] = inferences.get(0).getHessianVectorProduct(d);
        for (int i = 1; i < domain_cnt; i++) {
            double Hd_i[] = inferences.get(i).getHessianVectorProduct(d);
            for (int j = 0; j < Hd.length; j++) {
                Hd[j] += Hd_i[j];
            }
        }

        for(int i = 0; i < prior.length; ++i){
            Hd[i] += d[i] * prior[i];
        }

        return Hd;
    }

    private double[] getHessianVectorProductEM(double[] d, double[] prior) {
        double Hd[] = inferencesEM.get(0).getHessianVectorProduct(d);
        for (int i = 1; i < domain_cnt; i++) {
            double Hd_i[] = inferencesEM.get(i).getHessianVectorProduct(d);
            for (int j = 0; j < Hd.length; j++) {
                Hd[j] += Hd_i[j];
            }
        }

        for(int i = 0; i < Hd.length; ++i){
            Hd[i] += d[i] * prior[i];
        }

        return Hd;
    }

    // get variance of inferred counts for each formula
    private double[] getVariance(double[] prior) {
        double []variance = new double[numWts];
        Arrays.fill(variance, 0);
        for (int i = 0; i < domain_cnt; i++) {
            double trueCnts[] = inferences.get(i).numFormulaTrueCnts;
            double trueSqCnts[] = inferences.get(i).numFormulaTrueSqCnts;
            double lambdaTrueCnts = inferences.get(i).numLambdaTrueCnts;
            double lambdaTrueSqCnts = inferences.get(i).numLambdaTrueSqCnts;
            int numSamples = inferences.get(i).numIter;

            for (int formulaNum = 0; formulaNum < numFormulas; formulaNum++) {
                double x   = trueCnts[formulaNum];
                double xsq = trueSqCnts[formulaNum];

                // Add variance for this domain
                variance[formulaNum] += xsq/numSamples - (x/numSamples)*(x/numSamples);
            }
            if(priorSoftEvidence)
                variance[numFormulas] += (lambdaTrueSqCnts / numSamples) - ((lambdaTrueCnts / numSamples) * (lambdaTrueCnts / numSamples));
        }
        for (int formulaNum = 0; formulaNum < numFormulas; formulaNum++) {
            variance[formulaNum] += prior[formulaNum];
        }
        return variance;
    }

    private double[] getVarianceEM(double[] prior) {
        double []variance = new double[numWts];
        for (int i = 0; i < domain_cnt; i++) {
            double trueCnts[] = inferencesEM.get(i).numFormulaTrueCnts;
            double trueSqCnts[] = inferencesEM.get(i).numFormulaTrueSqCnts;
            double lambdaTrueCnts = inferencesEM.get(i).numLambdaTrueCnts;
            double lambdaTrueSqCnts = inferencesEM.get(i).numLambdaTrueSqCnts;
            int numSamples = inferencesEM.get(i).numIter;

            for (int formulaNum = 0; formulaNum < numFormulas; formulaNum++) {
                double x   = trueCnts[formulaNum];
                double xsq = trueSqCnts[formulaNum];
                // Add variance for this domain
                variance[formulaNum] += xsq/numSamples - (x/numSamples)*(x/numSamples);
            }
            if(priorSoftEvidence) {
                variance[numFormulas] += (lambdaTrueSqCnts / numSamples) - ((lambdaTrueCnts / numSamples) * (lambdaTrueCnts / numSamples));
            }
        }
        for (int formulaNum = 0; formulaNum < numFormulas; formulaNum++) {
            variance[formulaNum] += prior[formulaNum];
        }
        return variance;
    }

    private double dotprod(double[] v1, double[] v2, int numWts) {
        double total = 0.0;
        for (int i = 0; i < numWts ; i++) {
            total += v1[i]*v2[i];
        }
        return total;
    }

    private void findFormulaTrainCnts() {
        for (int i = 0; i < domain_cnt; i++) {
            GroundMLN gm = inferences.get(i).state.groundMLN;
            Evidence truth = inferences.get(i).truth;
            for(GroundFormula gf : gm.groundFormulas)
            {
                List<Integer> parentFormulaId = gf.parentFormulaId;
                List<Integer> numCopies = gf.numCopies;
                if(parentFormulaId.get(0) == -1)
                    continue;
                boolean isFormulaSatisfied = true;
                for(GroundClause gc : gf.groundClauses)
                {
                    boolean isClauseSatisfied = false;
                    for(int gpId : gc.groundPredIndices)
                    {
                        int trueVal = 0;
                        if(truth.predIdVal.containsKey(gpId)) {
                            trueVal = truth.predIdVal.get(gpId);
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
                    formulaNumGroundings[i][parentFormulaId.get(j)] += numCopies.get(j);
                }


            }
        }
    }

    private void findFormulaTrainCountsLambda() {
        Arrays.fill(lambdaTrainCnts,0.0);
        for (int i = 0; i < domain_cnt; i++) {
            GroundMLN gm = inferences.get(i).state.groundMLN;
            Evidence truth = inferences.get(i).truth;
            for(int gfId : inferences.get(i).state.groundedGfIndicesList)
            {
                GroundFormula gf = gm.groundFormulas.get(gfId);
                boolean isFormulaSatisfied = true;
                for(GroundClause gc : gf.groundClauses)
                {
                    boolean isClauseSatisfied = false;
                    for(int gpId : gc.groundPredIndices)
                    {
                        int trueVal = 0;
                        if(truth.predIdVal.containsKey(gpId))
                            trueVal = truth.predIdVal.get(gpId);
                        BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                        isClauseSatisfied |= b.get(trueVal);
                        if(isClauseSatisfied)
                            break;
                    }
                    isFormulaSatisfied &= isClauseSatisfied;
                    if(!isFormulaSatisfied)
                        break;
                }
                if(isFormulaSatisfied) {
                    lambdaTrainCnts[i] += gf.originalWeight.getValue();
                }
            }
        }
    }

    private double findGradientForLambda(double[] prior) {
        double gradient = 0.0;
        System.out.println("SeLambda: "+weights[numFormulas]);
        for (int i = 0; i < domain_cnt; i++) {
            gradient -= findGradientForDomainForLambda(i);
        }
        gradient += (weights[numFormulas] - priorMeans[numFormulas]) * prior[numFormulas];
        return gradient/* / domain_cnt*/;
    }

    private double findGradientForDomainForLambda(int domainIndex) {
        double lambdaInferredCnts = inferences.get(domainIndex).numLambdaTrueCnts;
        lambdaInferredCnts /= inferences.get(domainIndex).numIter;
        if(withEM){
            double lambdaInferredCntsEM = 0;
            lambdaInferredCntsEM = inferencesEM.get(domainIndex).numLambdaTrueCnts;
            lambdaInferredCntsEM /= inferencesEM.get(domainIndex).numIter;
            return lambdaInferredCntsEM - lambdaInferredCnts;
        } else {
            System.out.println("FormulaNum\tActual Count\tInferred Count");
//            System.out.println("lambda\t\t" + lambdaTrainCnts[domainIndex] + "\t" + lambdaInferredCnts);
            double t1 = Math.round(lambdaTrainCnts[domainIndex] * 1000) / 1000.0;
            double t2 = Math.round(lambdaInferredCnts * 1000) / 1000.0;
            System.out.println("lambda\t\t" + t1 + "\t" + t2);

            return t1 - t2;
        }
    }

    private void findGradient(double[] prior) {
        Arrays.fill(gradient,0.0);
        for (int i = 0; i < domain_cnt ; i++) {
            if(dldebug)
                System.out.println("Finding gradient for domain " + i);
            getGradientForDomain(gradient, i);
        }

        //Adding gradient of Prior
        for (int i = 0; i < numFormulas; i++) {
            gradient[i] += (weights[i] - priorMeans[i]) * prior[i];
        }
        if(priorSoftEvidence){
            gradient[numFormulas] = findGradientForLambda(prior);
        }
    }

    private void getGradientForDomain(double []gradient, int domainIndex) {
        double []formulaInferredCnts = inferences.get(domainIndex).numFormulaTrueCnts;
        double []formulaInferredCntsEM = null;
        if(withEM)
            formulaInferredCntsEM = inferencesEM.get(domainIndex).numFormulaTrueCnts;

        if(dldebug) {
            if(withEM)
                System.out.println("FormulaNum\tEM Count\tInferred Count");
            else
                System.out.println("FormNum\tTtlGrnd\tActlCnt\tInfrdCnt");
        }

        for (int j = 0; j < numFormulas; j++) {
            double inferredCount = formulaInferredCnts[j]/inferences.get(domainIndex).numIter;
            double inferredCountEM = 0;
            if(withEM)
                inferredCountEM = formulaInferredCntsEM[j]/inferencesEM.get(domainIndex).numIter;
            if(dldebug)
            {
                if(withEM)
                    System.out.println(j + "\t" + inferredCountEM + "\t" + inferredCount);
                else
                    System.out.println(j + "\t" + formulaNumGroundings[domainIndex][j] + "\t" + formulaTrainCnts[domainIndex][j] + "\t" + inferredCount);
            }
            if(withEM)
                gradient[j] -= (inferredCountEM - inferredCount);
            else
                gradient[j] -= (formulaTrainCnts[domainIndex][j] - inferredCount);
        }
    }

    private void infer(boolean burningIn, boolean isInit) {
        InferPerMLN ipmln;
        Thread t[] = new Thread[domain_cnt];
        MLN mln = inferences.get(0).mln;
        Thread.currentThread().setPriority(10);
        for (int i = 0; i < domain_cnt; i++) {
            State state = inferences.get(i).state;
            if(priorSoftEvidence) {
                state.groundMLN.setGroundFormulaWtsToParentWtsSoftEvidence(mln, weights[numFormulas]);
            } else {
                state.groundMLN.setGroundFormulaWtsToParentWts(mln);
            }
            inferences.get(i).updateWtsForNextGndPred(0);
            ipmln = new InferPerMLN(inferences.get(i), burningIn, isInit);
            t[i] = new Thread(ipmln);
            t[i].setPriority(1);
        }

        try {
            for (int i = 0; i < domain_cnt; ++i){
                System.out.println("Doing inference for domain " + i);
                t[i].start();
            }
            System.out.println("Waiting for threads to finish.");
            Thread.currentThread().setPriority(1);
            for (int i = 0; i < domain_cnt; ++i){
                t[i].join();
            }
            Thread.currentThread().setPriority(10);
        } catch (InterruptedException e) {
            System.out.println("Main thread Interrupted");
        }

//        MLN mln = inferences.get(0).mln;
//        for (int i = 0; i < domain_cnt; i++) {
//            State state = inferences.get(i).state;
//            state.setGroundFormulaWtsToParentWts(mln);
//            inferences.get(i).updateWtsForNextGndPred(0);
//            System.out.println("Doing inference for domain " + i);
//            inferences.get(i).infer(burningIn, isInit);
//        }
    }

    private void inferEM(boolean burningIn, boolean isInit) {
        InferPerMLN ipmln;
        Thread t[] = new Thread[domain_cnt];
        MLN mln = inferencesEM.get(0).mln;
        Thread.currentThread().setPriority(10);
        for (int i = 0; i < domain_cnt; i++) {
            State state = inferencesEM.get(i).state;
            if(priorSoftEvidence) {
                state.groundMLN.setGroundFormulaWtsToParentWtsSoftEvidence(mln, weights[numFormulas]);
            } else {
                state.groundMLN.setGroundFormulaWtsToParentWts(mln);
            }
            inferencesEM.get(i).updateWtsForNextGndPred(0);
            ipmln = new InferPerMLN(inferencesEM.get(i), burningIn, isInit);
            t[i] = new Thread(ipmln);
            t[i].setPriority(1);
        }

        try {
            for (int i = 0; i < domain_cnt; ++i){
                System.out.println("Doing inference for EM domain " + i);
                t[i].start();
            }
            System.out.println("Waiting for threads to finish.");
            Thread.currentThread().setPriority(1);
            for (int i = 0; i < domain_cnt; ++i){
                t[i].join();
            }
            Thread.currentThread().setPriority(10);
        } catch (InterruptedException e) {
            System.out.println("Main thread Interrupted");
        }

//        MLN mln = inferencesEM.get(0).mln;
//        for (int i = 0; i < domain_cnt; i++) {
//            State state = inferencesEM.get(i).state;
//            state.setGroundFormulaWtsToParentWts(mln);
//            inferencesEM.get(i).updateWtsForNextGndPred(0);
//            System.out.println("Doing inference in EM for domain " + i);
//            inferencesEM.get(i).infer(burningIn, isInit);
//        }
    }

    private void setMLNWeights() {
        MLN mln = inferences.get(0).mln;
        MLN mlnEM = null;
        if(withEM)
            mlnEM = inferencesEM.get(0).mln;
        for (int i = 0; i < numFormulas ; i++)
        {
            mln.formulas.get(i).weight = new LogDouble(weights[i], true);
            if(withEM)
                mlnEM.formulas.get(i).weight = new LogDouble(weights[i], true);
        }
    }


}
