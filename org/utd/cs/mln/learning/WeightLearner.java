package org.utd.cs.mln.learning;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.mln.alchemy.core.MLN;
import org.utd.cs.mln.alchemy.util.Parser;

import java.io.*;
import java.util.*;

/**
 * This is an abstract base class (can't instantiate it) for weight learning classes such as {@link DiscLearner} and {@link GenLearner}.
 * It stores data structures required for every weight learning algorithm.
 *
 * @author Happy
 * @since 03/28/18
 */
public abstract class WeightLearner {
    List<MLN> mlns, mlnsEM;

    // List of weights. Its size is number of formulas in MLN (+1 if priorSoftEvidence). If priorSoftEvidence, last
    // element is lambda learned.
    double [] weights, prior, priorMeans, priorStdDevs;
    public int domain_cnt, numWts;
    public boolean priorSoftEvidence;
    public int numFormulas;
    public boolean withEM;
    /**
     * Constructor : Initializes the fields
     * @param mlnsParam List of MLNs of size number of databases. All MLNs are same except for domains for predicates.
     */
    public WeightLearner(List<MLN> mlnsParam, LearnArgs lArgs) {
        mlns = new ArrayList<>();
        for(MLN mln : mlnsParam)
        {
            this.mlns.add(mln);
        }
        numFormulas = mlns.get(0).formulas.size();
        domain_cnt = mlns.size();
        numWts = numFormulas;
        if(lArgs.priorSoftEvidence)
        {
            priorSoftEvidence = true;
            numWts++;
        }
        weights = new double[numWts];
        priorMeans = new double[numWts];
        priorStdDevs = new double[numWts];
        prior = new double[numWts];
        if(priorSoftEvidence)
            weights[numWts-1] = lArgs.seLambda;
        if(lArgs.usePrior)
        {
            for (int i = 0; i < numFormulas; i++) {
                weights[i] = mlns.get(0).formulas.get(i).weight.getValue();
            }
        }
    }

    /**
     * Subclass must implement this method, which fills in list of weights.
     */
    public abstract void learnWeights();

    /**
     * Set MLN formula weights to learned weights.
     */
    void setMLNWeights()
    {
        // Since all mlns are same upto domain difference, we update only first MLN's weights
        MLN mln = mlns.get(0);
        MLN mlnEM = null;
        if(withEM)
             mlnEM = mlnsEM.get(0);
        for (int i = 0; i < numFormulas ; i++)
        {
            mln.formulas.get(i).weight = new LogDouble(weights[i], true);
//            if(withEM)
//                mlnEM.formulas.get(i).weight = new LogDouble(weights[i], true);
        }
    }

//    public double[] setPrior() {
//        double[] prior = new double[numWts];
//        for (int i = 0; i < numWts; ++i){
//            int temp = 0;
//            for (int j = 0; j < domain_cnt; ++j){
//                if(i == numFormulas){
////                    temp += inferences.get(j).state.groundMLN.sePredCount;
//                } else {
//                    temp += formulaTrainCnts[j][i];
////                    temp += formulaNumGroundings[j][i];
//                }
//            }
////            temp /= domain_cnt;
////            temp = 1;
//            prior[i] = temp / (priorStdDevs[i] * priorStdDevs[i]);
//        }
//        return prior;
//    }

    /**
     * Write output MLN file with learned weights.
     * @param mln_file Input MLN file, required to copy formulas
     * @param out_file Output file to which it will write
     * @throws FileNotFoundException
     */
    public void writeWeights(String mln_file, String out_file) throws FileNotFoundException {
        Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(mln_file))));
        PrintWriter pw = new PrintWriter(out_file);
        boolean isformula = false;
        int formulaNum = 0;
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine().replaceAll("\\s", "");

            if(line.isEmpty() || line.contains(Parser.COMMENT)) {
                pw.write(line+"\n");
                continue;
            }
            if (line.contains("#formulas")) {
                pw.write("#formulas\n");
                isformula = true;
                continue;
            }
            if (isformula == false) {
                pw.write(line + "\n");
                continue;
            } else {
                String[] formulaArr = line.split(Parser.WEIGHTSEPARATOR);
                pw.printf(formulaArr[0] + Parser.WEIGHTSEPARATOR + "%.3f\n",weights[formulaNum]);
                formulaNum++;
            }
        }
        if(priorSoftEvidence)
            pw.write("#lambda::" + weights[numWts-1]);
        scanner.close();
        pw.close();
    }

}
