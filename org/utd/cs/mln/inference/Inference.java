package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.core.State;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

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

    // Whether Save all counts for all samples or not i.e. whether fill allFormulaTrueCnts and oldAllFormulaTrueCnts
    boolean saveAllCounts;

    public State state;

    int seed;

    // Indicates if need to store true counts for each first-order formula
    boolean trackFormulaTrueCnts;

    // sum of true counts and true squared counts of first order formulas over all samples.
    double [] formulaTrueCnts, formulaTrueSqCnts;

    List<Double> allLambdaTrueCnts, oldAllLambdaTrueCnts;

    double lambdaTrueCnts, lambdaTrueSqCnts;

    // allFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample
    // oldAllFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample in the previous iter of learning
    List<List<Double>> allFormulaTrueCnts, oldAllFormulaTrueCnts;

    // Number of samples taken of the true counts so far. This will be the size of
    int numSamples;
    Random rand;

    //Indicates whether softEvidence is given or not
    boolean priorSoftEvidence;

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
    public Inference(State state, int seed, boolean trackFormulaTrueCnts, boolean priorSoftEvidence)
    {
        this.state = state;
        this.seed = seed;
        if(seed == -1)
            this.seed = DEFAULT_SEED;
        this.trackFormulaTrueCnts = trackFormulaTrueCnts;
        this.priorSoftEvidence = priorSoftEvidence;
        this.saveAllCounts = false;
        this.numSamples = 0;
        this.rand = new Random();

        if(trackFormulaTrueCnts && state != null)
        {
            int numFormulas = state.mln.formulas.size();
            formulaTrueCnts = new double[numFormulas];
            formulaTrueSqCnts = new double[numFormulas];
        }
    }

    /**
     * Initializes the inference algorithm.
     */
    abstract void init();

    /**
     * Performs the inference algorithm.
     */
    abstract void infer();

    /**
     * Prints the probabilities of each predicate.
     */
    abstract void writeProbs(PrintWriter writer);

    /**
     * Computes Hessian vector product. Hessian is of the matrix allFormulaTrueCnts.
     * @param v vector to which Hessian is to be multiplied
     * @return resultant vector
     */
    public double[] getHessianVectorProduct(double[] v) {
        int numFormulas = state.mln.formulas.size();
        int numWts = numFormulas;
        if (priorSoftEvidence)
            numWts = numFormulas + 1;

        // For minimizing the negative log likelihood,
        // the ith element of H v is:
        //   E[n_i * vn] - E[n_i] E[vn]
        // where n is the vector of all clause counts
        // and vn is the dot product of v and n.

        double sumVN = 0;
        double []sumN = new double[numWts];
        double []sumNiVN = new double[numWts];

        // Get sufficient statistics from each sample,
        // so we can compute expectations
        int numSamples = allFormulaTrueCnts.size();
        for (int s = 0; s < numSamples; s++)
        {
            List<Double> n1 = allFormulaTrueCnts.get(s);
            double n2 = 0.0;
            if(priorSoftEvidence)
                n2 = allLambdaTrueCnts.get(s);

            // Compute v * n
            double vn = 0;

            for (int i = 0; i < numFormulas; i++)
                vn += v[i] * n1.get(i);

            if(priorSoftEvidence)
                vn += v[numFormulas] * n2;

            // Tally v*n, n_i, and n_i v*n
            sumVN += vn;
            for (int i = 0; i < numFormulas; i++)
            {
                sumN[i]    += n1.get(i);
                sumNiVN[i] += n1.get(i) * vn;
            }
            if(priorSoftEvidence){
                sumN[numFormulas] += n2;
                sumNiVN[numFormulas] += n2 * vn;
            }
        }

        // Compute actual product from the sufficient stats
        double []product = new double[numWts];
        for (int formulano = 0; formulano < numFormulas; formulano++)
        {
            double E_vn = sumVN/numSamples;
            double E_ni = sumN[formulano]/numSamples;
            double E_nivn = sumNiVN[formulano]/numSamples;
            product[formulano] = E_nivn - E_ni * E_vn;
        }
        if(priorSoftEvidence){
            double E_vn = sumVN/numSamples;
            double E_ni = sumN[numFormulas]/numSamples;
            double E_nivn = sumNiVN[numFormulas]/numSamples;
            product[numFormulas] = E_nivn - E_ni * E_vn;
        }

        return product;
    }

    /**
     * Resets all counts
     */
    public void resetCnts() {
        if(trackFormulaTrueCnts){
            if(priorSoftEvidence)
            {
                lambdaTrueCnts = 0.0;
                lambdaTrueSqCnts = 0.0;
            }

            Arrays.fill(formulaTrueCnts,0.0);
            Arrays.fill(formulaTrueSqCnts,0.0);
            numSamples = 0;
            if(saveAllCounts)
            {
                allFormulaTrueCnts = new ArrayList<List<Double>>();
                if(priorSoftEvidence)
                    allLambdaTrueCnts = new ArrayList<Double>();
            }

        }
    }

    /**
     * If we are storing all counts, then restore all counts from old counts
     */
    public void restoreCnts() {
        if (!saveAllCounts)
            return;

        resetCnts();
        System.out.println("Restoring counts...");
        if(trackFormulaTrueCnts){
            for (int i = 0; i < oldAllFormulaTrueCnts.size(); i++)
            {
                allFormulaTrueCnts.add(new ArrayList<Double>());
                int numcounts = oldAllFormulaTrueCnts.get(i).size();
                for (int j = 0; j < numcounts; j++)
                {
                    allFormulaTrueCnts.get(i).add(oldAllFormulaTrueCnts.get(i).get(j));
                    formulaTrueCnts[j] += allFormulaTrueCnts.get(i).get(j);
                    formulaTrueSqCnts[j] += allFormulaTrueCnts.get(i).get(j) * allFormulaTrueCnts.get(i).get(j);
                }
                if(priorSoftEvidence)
                {
                    allLambdaTrueCnts.add(oldAllLambdaTrueCnts.get(i));
                    lambdaTrueCnts += allLambdaTrueCnts.get(i);
                    lambdaTrueSqCnts += allLambdaTrueCnts.get(i) * allLambdaTrueCnts.get(i);
                }

                numSamples++;
            }
        }
    }

    /**
     * If we are storing all counts, save current all counts to old all counts
     */
    public void saveToOldCnts() {
        if (!saveAllCounts)
            return;

        if(trackFormulaTrueCnts){
            for (int i = 0; i < allFormulaTrueCnts.size(); i++)
            {
                int numcounts = allFormulaTrueCnts.get(i).size();
                for (int j = 0; j < numcounts; j++)
                {
                    oldAllFormulaTrueCnts.get(i).set(j, allFormulaTrueCnts.get(i).get(j));
                }
                if(priorSoftEvidence)
                    oldAllLambdaTrueCnts.set(i,allLambdaTrueCnts.get(i));
            }
        }
    }

    public void writeNetwork(PrintWriter writer) {
        writer.write(state.groundMLN.writableString());
    }
}
