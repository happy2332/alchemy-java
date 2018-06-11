package org.utd.cs.mln.learning;

import org.utd.cs.mln.alchemy.util.ArrayUtils;
import org.utd.cs.mln.alchemy.util.ConditionalLogLikelihood;
import org.utd.cs.mln.inference.GibbsSampler_v3;
import org.utd.cs.mln.inference.Inference;
import org.utd.cs.mln.inference.MCMC;

import java.util.Arrays;
import java.util.List;

/**
 * Created by Happy on 4/24/18.
 */
public class ConjugateGradient {
    public final double cg_max_lambda;
    public final boolean preConditionCG;
    public final int maxBackTracks;
    public final double minllchange;
    public final int numIter;
    public boolean backtracked;
    public double[] delta_pred, oldWeights, old_gradient, d, oldd, averageWeights;
    public double cg_lambda;
    public int backtrackCount;
    public double alpha;
    public boolean cgdebug = false;

    public ConjugateGradient(CGParams cgParam, int numWts) {
        delta_pred = null;
        cg_lambda = cgParam.cg_lambda;
        cg_max_lambda = cgParam.cg_max_lambda;
        preConditionCG = cgParam.preConditionCG;
        maxBackTracks = cgParam.maxBackTracks;
        oldWeights = new double[numWts];
        old_gradient = new double[numWts];
        d = new double[numWts];
        oldd = new double[numWts];
        averageWeights = new double[numWts];
        minllchange = cgParam.min_ll_change;
        numIter = cgParam.numIter;
    }

    public int updateWts(double[] weights, double[] gradient, List<Inference> inferences, List<Inference> inferencesEM, int iter, boolean withEM, ConditionalLogLikelihood cllLoss) {
        double realdist = 1.0;
        double preddist = 1.0;
        int numWts = weights.length;
        int domain_cnt = inferences.size();
        if (iter > 1 && delta_pred != null && !backtracked) {
            double []dist = new double[numWts];
            for (int i = 0; i < numWts; i++) {
                dist[i] = weights[i] - oldWeights[i];
            }
            double []avgPred = new double[numWts];
            for (int i = 0; i < numWts; i++) {
                avgPred[i] = old_gradient[i] + delta_pred[i] / 2.0;
            }
            preddist = ArrayUtils.dotprod(avgPred, dist);

            // Real change is lower bound on actual change
            realdist = ArrayUtils.dotprod(gradient, dist);
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
            if (delta < 0.0 && backtrackCount < maxBackTracks) {
                System.out.println("Backtracking...");
                for (int i = 0; i < numWts; i++)
                    weights[i] = oldWeights[i];

                for (int i = 0; i < domain_cnt; i++)
                {
                    ((GibbsSampler_v3)inferences.get(i)).restoreCnts();
                    if(withEM)
                        ((GibbsSampler_v3)inferencesEM.get(i)).restoreCnts();
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
                ((GibbsSampler_v3)inferences.get(i)).saveToOldCnts();
                if(withEM)
                    ((GibbsSampler_v3)inferencesEM.get(i)).saveToOldCnts();
            }

            double preCond[] = new double[numWts];
            Arrays.fill(preCond, 1.0);
            if(preConditionCG)
            {
                double variance[] = getVariance(inferences, numWts);
                if(withEM)
                {
                    double varianceEM[] = getVariance(inferencesEM, numWts);
                    for (int i = 0; i < variance.length; i++) {
                        variance[i] -= varianceEM[i];
                    }
                }

                if(cgdebug)
                    System.out.println("variance : " + Arrays.toString(variance));
                for (int formulaNum = 0; formulaNum < numWts; formulaNum++)
                    preCond[formulaNum] = 1.0/(variance[formulaNum] + 0.00001);
            }
            double beta;

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

            if(cgdebug)
                System.out.println("beta = " + beta);

            // Compute new direction
            for (int w = 0; w < numWts; w++) {
                d[w] = -preCond[w] * gradient[w] + beta * oldd[w];
            }
            System.out.println("Direction D = " + Arrays.toString(d));

            double Hd[] = getHessianVectorProduct(inferences, d, cllLoss);
            if(withEM)
            {
                double []HdEM = getHessianVectorProduct(inferencesEM,d, cllLoss);
                for (int i = 0; i < Hd.length; i++) {
                    Hd[i] -= HdEM[i];
                }
            }
            alpha = computeQuadraticStepLength(Hd, gradient);
            if (alpha < 0.0)
            {
                for (int w = 0; w < numWts; w++)
                    d[w] = -preCond[w]*gradient[w];
                Hd = getHessianVectorProduct(inferences, d, cllLoss);
                if(withEM)
                {
                    double []HdEM = getHessianVectorProduct(inferencesEM, d, cllLoss);
                    for (int i = 0; i < Hd.length; i++) {
                        Hd[i] -= HdEM[i];
                    }
                }
                alpha = computeQuadraticStepLength(Hd, gradient);
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
            double maxchange = -ArrayUtils.dotprod(gradient, wchange);
            System.out.println("Maximum Estimated Improvement = " + maxchange);
            if (maxchange < minllchange) {
                System.out.println("Upper bound is less than " + minllchange + ", halting learning.");
                return -1;
            }

            // Save weights, gradient, and direction and adjust the weights
            for (int w = 0; w < numWts; w++) {
                oldWeights[w] = weights[w];
                oldd[w] = d[w];
                old_gradient[w] = gradient[w];
                weights[w] += wchange[w];

                averageWeights[w] = ((iter - 1) * averageWeights[w] + weights[w]) / iter;
            }
            delta_pred = getHessianVectorProduct(inferences, wchange, cllLoss);
            //TODO : delta pred for EM ?

            System.out.println("weights = " + Arrays.toString(weights));
            System.out.println("Avg weights = " + Arrays.toString(averageWeights));
        }
        return 0;
    }

    private double[] getVariance(List<Inference> inferencesParam, int numWts) {
        int domain_cnt = inferencesParam.size();
        double variance[] = new double[numWts];
        for (int domainId = 0; domainId < domain_cnt; domainId++) {
            double []trueCnts = inferencesParam.get(domainId).formulaTrueCnts;
            double []trueSqCnts = inferencesParam.get(domainId).formulaTrueSqCnts;
            for (int formulaId = 0; formulaId < trueCnts.length; formulaId++) {
                variance[formulaId] += trueSqCnts[formulaId] - trueCnts[formulaId]*trueCnts[formulaId];
            }
        }
        return variance;
    }

    private double[] getHessianVectorProduct(List<Inference> inferencesParam, double[] d, ConditionalLogLikelihood cllLoss) {
        int domain_cnt = inferencesParam.size();
        double Hd[] = ((MCMC)inferencesParam.get(0)).getHessianVectorProduct(d);
        for (int i = 1; i < domain_cnt; i++) {
            double Hd_i[] = ((MCMC)inferencesParam.get(i)).getHessianVectorProduct(d);
            for (int j = 0; j < Hd.length; j++) {
                Hd[j] += Hd_i[j];
            }
        }

        for (int i = 0; i < Hd.length; i++) {
            Hd[i] += d[i] * cllLoss.priorLambda[i]/(cllLoss.priorStdDevs[i]*cllLoss.priorStdDevs[i]);
        }
        return Hd;
    }

    private double computeQuadraticStepLength(double[] hd, double gradient[]) {
        int numWeights = d.length;
        // Compute step length using trust region approach
        double dHd = ArrayUtils.dotprod(d, hd);
        double dd = ArrayUtils.dotprod(d, d);
        double dg = ArrayUtils.dotprod(gradient, d);
        double alpha_denom = (dHd + cg_lambda * dd);
        if(alpha_denom == 0.0)
            return 0.0;
        double alpha = -dg/alpha_denom;

        if(cgdebug)
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

}
