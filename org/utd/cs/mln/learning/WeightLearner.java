package org.utd.cs.mln.learning;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.MyAssert;
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
    // List of states. Each state is corresponding to a domain.
    List<State> states, statesEM;

    // List of weights. Its size is number of formulas in MLN (+1 if priorSoftEvidence). If priorSoftEvidence, last
    // element is lambda learned.
    double [] weights;
    // For each domain d, for each ground predicate gp in that d, stores softevidence weight (without multiplication with lambda) for each possible value val
    // of gp
    // softEvidencePerPredPerVal[d][gp][val]
    List<Map<Integer, List<Double>>> softEvidencePerPredPerVal;

    // priorLambda is lambda for weight regularization. For each weight, there can be different lambda. If an entry of prior
    // is zero, it means there is no regularization on that corresponding weight.
    // priorMeans and priorStdDevs are the mean and standard deviation vectors for gaussian prior. Note that standard
    // deviation is a vector here, meaning its a diagonal matrix.
    // By default, priorMeans = 0, and priorStdDevs = 1, in that case, gaussian prior is just l2 regularization.
    double [] priorLambda, priorMeans, priorStdDevs;
    public int domain_cnt, numWts;
    public boolean priorSoftEvidence;
    public int numFormulas;
    public boolean withEM;
    /**
     * Constructor : Initializes the fields
     * @param mlnsParam List of MLNs of size number of databases. All MLNs are same except for domains for predicates.
     */
    public WeightLearner(List<MLN> mlnsParam, List<GroundMLN> groundMlnsParam, List<GroundMLN> groundMlnsEMParam,
                         List<Evidence> truthsParam, List<Evidence> truthsEMParam, LearnArgs lArgs) throws FileNotFoundException {

        withEM = lArgs.withEM;
        //create states
        states = new ArrayList<>();
        if(withEM)
            statesEM = new ArrayList<>();
        domain_cnt = mlnsParam.size();
        for (int i = 0; i < domain_cnt; i++) {
            states.add(new State(mlnsParam.get(i), groundMlnsParam.get(i), truthsParam.get(i)));
            if(withEM)
                statesEM.add(new State(mlnsParam.get(i), groundMlnsEMParam.get(i), truthsEMParam.get(i)));
        }

        numFormulas = states.get(0).mln.formulas.size();
        numWts = numFormulas;
        if(lArgs.priorSoftEvidence)
        {
            priorSoftEvidence = true;
            numWts++;
        }
        weights = new double[numWts];
        priorMeans = new double[numWts];
        priorStdDevs = new double[numWts];
        Arrays.fill(priorStdDevs,2.0);
        priorLambda = new double[numWts];
        Arrays.fill(priorLambda,1.0);

        if(priorSoftEvidence){
            weights[numWts-1] = lArgs.seLambda;
            softEvidencePerPredPerVal = new ArrayList<>();
            createSoftEvidencePerPredPerVal(lArgs);
        }

        // If useMlnWts is true, then take initial weights as specified in mln file.
        if(lArgs.useMlnWts)
        {
            for (int i = 0; i < numFormulas; i++) {
                weights[i] = states.get(0).mln.formulas.get(i).weight.getValue();
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
        for (int domainId = 0; domainId < domain_cnt; domainId++) {
            MLN mln = states.get(domainId).mln;

            for (int i = 0; i < numFormulas ; i++)
            {
                mln.formulas.get(i).weight = new LogDouble(weights[i], true);
            }
            if(priorSoftEvidence)
                mln.softEvidenceLambda = weights[numWts-1];
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

    public void setPriorLambda()
    {
        Arrays.fill(priorLambda, 0.0);
    }

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

    private void createSoftEvidencePerPredPerVal(LearnArgs lArgs) throws FileNotFoundException {
        for (int domainId = 0; domainId < domain_cnt; domainId++) {
            State state = states.get(domainId);
            String softEvidenceFile = lArgs.softEvidenceFiles.get(domainId);
            String sePredName = lArgs.sePred;
            Map<Integer, List<Double>> softEvidencePerPredPerValPerDomain = new HashMap<>();
            Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(softEvidenceFile))));
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine().replaceAll("\\s", "");
                if (line.isEmpty()) {
                    continue;
                }
                String words[] = line.split(",");
                int constant = Integer.parseInt(words[0]);
                int value = Integer.parseInt(words[1]);
                double weight = Double.parseDouble(words[2]);
                GroundPredicate gp = new GroundPredicate();
                for(PredicateSymbol ps : states.get(0).mln.symbols)
                {
                    if(ps.symbol.equals(sePredName))
                    {
                        gp.symbol = new GroundPredicateSymbol(ps.id, ps.symbol, ps.values, ps.variable_types);
                        break;
                    }
                }
                gp.numPossibleValues = gp.symbol.values.values.size();
                gp.terms.add(constant);
                MyAssert.assume(state.groundMLN.groundPredToIntegerMap.containsKey(gp));
                int gpIndex = state.groundMLN.groundPredToIntegerMap.get(gp);
                if(!softEvidencePerPredPerValPerDomain.containsKey(gpIndex))
                    softEvidencePerPredPerValPerDomain.put(gpIndex, new ArrayList<Double>(Collections.nCopies(gp.numPossibleValues, 0.0)));
                softEvidencePerPredPerValPerDomain.get(gpIndex).set(value, weight);
            }
            softEvidencePerPredPerVal.add(softEvidencePerPredPerValPerDomain);
        }
    }

}
