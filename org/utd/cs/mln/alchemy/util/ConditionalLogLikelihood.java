package org.utd.cs.mln.alchemy.util;

import org.utd.cs.mln.inference.Inference;
import org.utd.cs.mln.learning.LearnTest;

import java.util.Arrays;
import java.util.List;

/**
 * Created by Happy on 4/24/18.
 */
public class ConditionalLogLikelihood extends Loss {
    private double []gradient;
    private List<Inference> inferences, inferencesEM;
    public boolean cllDebug = false;
    private int numWts;
    private int domain_cnt;
    public double[][] formulaTrainCnts;
    public boolean agg = false;

    public ConditionalLogLikelihood(List<Inference> inferencesParam, List<Inference> inferencesEMParam, boolean withEMParam, boolean usePriorParam, boolean aggParam) {
        super();
        init(inferencesParam, inferencesEMParam, withEMParam);
        if(usePriorParam)
        {
            Arrays.fill(priorLambda,1.0);
        }
        agg = aggParam;
    }

    private void init(List<Inference> inferencesParam, List<Inference> inferencesEMParam, boolean withEMParam) {
        this.inferences = inferencesParam;
        if(withEMParam)
            this.inferencesEM = inferencesEMParam;
        numWts = inferences.get(0).state.mln.formulas.size();
        if(inferences.get(0).state.mln.priorSoftEvidence)
        {
            numWts++;
        }
        super.init(numWts);
        gradient = new double[numWts];
        domain_cnt = inferences.size();
    }

    @Override
    public double getLossValue() {
        return 0;
    }

    @Override
    public double[] getGradient() {
        return new double[0];
    }

    /**
     * This method computes gradient of negative log likelihood of data and fills in with gradient field.
     * g[i] = (E[n_i] - n_i) + reg[i], where n_i is number of true groundings of ith formula in training data,
     * E[n_i] is expected number of true groundings of ith formula, where expectation is over P(y|x),
     * and reg[i] comes from regularization, which is priorLambda[i]*(weights[i] - priorMeans[i])/priorStdDev[i]^2.
     * n_i is stored in formulaTrainCnts[i]
     * To find E[n_i], run inference, and store results in formulaInferedCnts i.e. E[n_i] is formulaInferedCnts[i].
     */
    public double[] getGradient(double []weights) {
        // For finding gradient, we need to do inference first so that we can calculate expected true counts.
        infer(inferences);
        Arrays.fill(gradient, 0.0);
        if(cllDebug)
            System.out.println("formula\tActual Count\tInfered Count");
        for (int w = 0; w < numWts; w++) {
            for (int domainId = 0; domainId < domain_cnt; domainId++) {
                gradient[w] += (inferences.get(domainId).formulaTrueCnts[w] - formulaTrainCnts[domainId][w]);
                double wt = weights[w];
                int numFormulas = inferences.get(domainId).state.mln.formulas.size();
                if(agg)
                {
                    if(w < numFormulas)
                        wt = wt/ LearnTest.firstOrderNumConnections.get(w);
                }
                double regularization = priorLambda[w]*(wt - priorMeans[w])/(priorStdDevs[w]*priorStdDevs[w]);
                gradient[w] += regularization;
                if(cllDebug)
                    System.out.println(w + "\t" + formulaTrainCnts[domainId][w] + "\t" + inferences.get(domainId).formulaTrueCnts[w]);
            }
            //double regularization = priorLambda[w]*(weights[w] - priorMeans[w])/(priorStdDevs[w]*priorStdDevs[w]);
            //gradient[w] += regularization;

        }
        return gradient;
    }

    /**
     * This method runs inference for all domains in parallel.
     */
    public void infer(List<Inference> inferences) {
        ConditionalLogLikelihood.InferPerDomain inferPerDomain;
        // Whether to run inference in parallel or not.
        boolean withThreads = true;
        // Run inferences in parallel
        if(withThreads)
        {
            // Create separate thread for each domain
            Thread t[] = new Thread[domain_cnt];
            Thread.currentThread().setPriority(10);
            for (int i = 0; i < domain_cnt; i++) {
                inferPerDomain = new ConditionalLogLikelihood.InferPerDomain(i, inferences);
                t[i] = new Thread(inferPerDomain);
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
        }
        else
        {
            for (int i = 0; i <domain_cnt ; i++) {
                inferPerDomain = new InferPerDomain(i, inferences);
                inferPerDomain.run();
            }
        }
    }

    private class InferPerDomain implements Runnable{
        private final int domainId;
        private final List<Inference> inferences;

        public InferPerDomain(int domainId, List<Inference> inferences)
        {
            this.domainId = domainId;
            this.inferences = inferences;
        }

        @Override
        public void run() {
            inferences.get(domainId).infer();
        }
    }
}
