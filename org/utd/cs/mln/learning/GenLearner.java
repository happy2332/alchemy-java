package org.utd.cs.mln.learning;
import org.utd.cs.gm.utility.*;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
import org.utd.cs.mln.alchemy.util.OtherUtils;
import util.LBFGS;
import util.U;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * This class implements generative learning of MLN by maximizing pseudologlikelihood of data.
 * The algorithm used is LBFGS, for which we use @see <a href="https://github.com/brendano/myutil">this library</a>.
 * @author Happy
 * @since 27/03/18
 */
public class GenLearner extends WeightLearner{

    // List of states. Each state is corresponding to a domain.
    private List<State> states;

    // For each domain d, for each ground predicate gp in that d, for each first order formula f in which gp appears,
    // stores number of satisfied groundings of f for each possible value val of gp.
    // satCounts[d][gp][f][val]
    private List<List<Map<Integer, List<Integer>>>> satCounts;

    // For each domain d, for each ground predicate gp in that d, stores softevidence weight for each possible value val
    // of gp
    // softEvidencePerPredPerVal[d][gp][val]
    private List<Map<Integer, List<Double>>> softEvidencePerPredPerVal;

    // For each domain d, stores number of ground atoms for each first order predicate p.
    // predToNumGndingsMap[d][p]
    private List<Map<String, Integer>> predToNumGndingsMap;
    private boolean genLearnDebug = false;

    private double []gradient;
    private long time;


    public GenLearner(List<MLN> mlns, List<GroundMLN> groundMlns, List<Evidence> truths, LearnArgs lArgs) throws FileNotFoundException {
        super(mlns, lArgs);
        init(groundMlns, truths, lArgs);
    }

    void init(List<GroundMLN> groundMlns, List<Evidence> truths, LearnArgs lArgs) throws FileNotFoundException {
        states = new ArrayList<>();
        satCounts = new ArrayList<>();
        predToNumGndingsMap = new ArrayList<>();
        for (int i = 0; i < domain_cnt; i++) {
            this.states.add(new State(mlns.get(i), groundMlns.get(i), truths.get(i)));
            predToNumGndingsMap.add(new HashMap<String,Integer>());
        }
        gradient = new double[numWts];
        createPredToNumGndingsMap();
        long time = System.currentTimeMillis();
        countSatVals();
        System.out.println("Counting took " + org.utd.cs.gm.utility.Timer.time((System.currentTimeMillis() - time)/1000.0));
        if(priorSoftEvidence){
            softEvidencePerPredPerVal = new ArrayList<>();
            createSoftEvidencePerPredPerVal(lArgs);
        }
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
                for(PredicateSymbol ps : mlns.get(0).symbols)
                {
                    if(ps.symbol.equals(sePredName))
                    {
                        gp.symbol = new GroundPredicateSymbol(ps.id, ps.symbol, ps.values, ps.variable_types);
                        break;
                    }
                }
                gp.numPossibleValues = gp.symbol.values.values.size();
                gp.terms.add(constant);
                assert state.groundMLN.groundPredToIntegerMap.containsKey(gp);
                int gpIndex = state.groundMLN.groundPredToIntegerMap.get(gp);
                if(!softEvidencePerPredPerValPerDomain.containsKey(gpIndex))
                    softEvidencePerPredPerValPerDomain.put(gpIndex, new ArrayList<Double>(Collections.nCopies(gp.numPossibleValues, 0.0)));
                softEvidencePerPredPerValPerDomain.get(gpIndex).set(value, weight);
            }
            softEvidencePerPredPerVal.add(softEvidencePerPredPerValPerDomain);
        }
    }

    private void createPredToNumGndingsMap() {
        for(int domainId = 0 ; domainId < domain_cnt ; domainId++)
        {
            Map<String,Integer> predToNumGndingsMapPerDomain = predToNumGndingsMap.get(domainId);
            GroundMLN groundMln = states.get(domainId).groundMLN;
            List<GroundPredicate> groundPreds = groundMln.groundPredicates;
            for(GroundPredicate gp : groundPreds)
            {
                String predName = gp.symbol.symbol;
                if(!predToNumGndingsMapPerDomain.containsKey(predName))
                {
                    predToNumGndingsMapPerDomain.put(predName,0);
                }
                int oldVal = predToNumGndingsMapPerDomain.get(predName);
                predToNumGndingsMapPerDomain.put(predName,oldVal+1);
            }
        }
    }


    /**
     * For each domain, for each groundPred p, create Map of satCounts. Key : First order formulaId in which p occurs,
     * value : List of true groundings of this formulaId for each possible value of p
     */

    void countSatVals() {
        // for each domain
        for (int domainId = 0 ; domainId < domain_cnt ; domainId++) {
            System.out.println("Finding counts for domain "+(domainId+1));
            State state = states.get(domainId);
            List<Map<Integer, List<Integer>>> satCountsPerDomain = new ArrayList<>();
            for (int i = 0; i < state.groundMLN.groundPredicates.size(); i++) {
                satCountsPerDomain.add(new HashMap<Integer, List<Integer>>());
            }
            MLN mln = mlns.get(domainId);
            // for each formula, compute sat counts
            for (int formulaId = 0; formulaId < mln.formulas.size(); formulaId++) {
                Formula formula = mln.formulas.get(formulaId);
                Set<Term> formulaWiseTermToGround = new HashSet<Term>();
                for (WClause clause : formula.clauses) {
                    for (Atom atom : clause.atoms) {
                        for (int j = 0; j < atom.terms.size(); j++) {
                            Term term = atom.terms.get(j);
                            formulaWiseTermToGround.add(term);
                        }
                    }
                }
                if(genLearnDebug)
                    System.out.println("Finding counts for formula : " + formulaId);
                countSatValsPerFormula(domainId, formulaId, satCountsPerDomain, new ArrayList<Term>(formulaWiseTermToGround));
            }
            satCounts.add(satCountsPerDomain);
        }
    }


    private void countSatValsPerFormula(int domainId, int formulaId, List<Map<Integer, List<Integer>>> satCountsPerDomain, ArrayList<Term> terms) {
        State state = states.get(domainId);
        MLN mln = mlns.get(domainId);
        int[][] permutations = FullyGrindingMill.permute(terms);
        if(genLearnDebug)
            System.out.println("Number of permutations : " + permutations.length);
        List<GroundPredicate> groundPredicates = state.groundMLN.groundPredicates;
        Map<GroundPredicate, Integer> groundPredToIntegerMap = state.groundMLN.groundPredToIntegerMap;
        Formula formula = mln.formulas.get(formulaId);
        for (int i = 0; i < permutations.length; i++) {
            if(genLearnDebug)
            {
                if (i % 100000 == 0)
                    System.out.println("i : " + i);
            }
            GroundFormula newFormula = new GroundFormula();
            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
            Set<Integer> formulaGpIndices = new HashSet<>();
            for (int c = 0; c < formula.clauses.size(); c++) {
                WClause clause = formula.clauses.get(c);
                GroundClause newGroundClause = new GroundClause();
                Map<Integer, BitSet> gpIndexToSatVals = new HashMap<>();
                // Iterate over each first order atom, and create ground atom for it.
                for (int j = 0; j < clause.atoms.size(); j++) {
                    boolean sign = clause.sign.get(j);
                    Atom oldAtom = clause.atoms.get(j); // first order atom
                    int valTrue = clause.valTrue.get(j);
                    GroundPredicate gp = new GroundPredicate(); // GroundPredicate to create
                    gp.symbol = new GroundPredicateSymbol(oldAtom.symbol.id, oldAtom.symbol.symbol, oldAtom.symbol.values, oldAtom.symbol.variable_types);
                    // Fill in the terms with constants
                    for (Term term : oldAtom.terms) {
                        int termIndex = terms.indexOf(term);
                        gp.terms.add(permutations[i][termIndex]);
                    }

                    int numPossibleValues = oldAtom.symbol.values.values.size();
                    gp.numPossibleValues = numPossibleValues;
                    int gpIndex = groundPredToIntegerMap.get(gp);
                    formulaGpIndices.add(gpIndex);

                    // Check if this groundPredicate occurs first time in this ground clause. then update
                    // groundClause's data structures about this groundPredicate.
                    int gpIndexInClause = newGroundClause.groundPredIndices.indexOf(gpIndex);
                    //GroundAtom newGroundAtom = new GroundAtom(gpIndex, gpIndexInClause, valTrue, sign);
                    if (gpIndexInClause == -1) {
                        newGroundClause.groundPredIndices.add(gpIndex);
                        gpIndexInClause = newGroundClause.groundPredIndices.size() - 1;
                        newGroundClause.globalToLocalPredIndex.put(gpIndex, gpIndexInClause);
                        gpIndexToSatVals.put(gpIndexInClause, new BitSet(gp.numPossibleValues));
                    }

                    // Now once we have added new ground Atom, we need to check if ground clause gets satisfied or not.
                    BitSet gpBitSet = new BitSet(gp.numPossibleValues);
                    gpBitSet.set(valTrue);
                    if (sign == true)
                        gpBitSet.flip(0, gp.numPossibleValues);
                    gpBitSet.or(gpIndexToSatVals.get(gpIndexInClause));

                    // If all bits are set for this groundPred, then this clause will always be satisfied and hence,
                    // shouldn't be added into groundformula. Note that, although at this point, we know that
                    // this clause shouldn't be added, but still we shouldn't just break out of this loop, as we
                    // need to add groundPredicates, but we shouldn't add any clauseInfo into groundPredicates appearing
                    // in this clause.

                    gpIndexToSatVals.put(gpIndexInClause, gpBitSet);
                }
                for (int gpId = 0; gpId < newGroundClause.groundPredIndices.size(); gpId++) {
                    BitSet b = gpIndexToSatVals.get(gpId);
                    newGroundClause.grounPredBitSet.add(b);
                }
                newGroundClauseList.add(newGroundClause);
            }
            newFormula.groundClauses.addAll(newGroundClauseList);

            for (Integer gpIndex : formulaGpIndices) {
                int numPossibleValues = groundPredicates.get(gpIndex).numPossibleValues;
                int tempTrueValue = state.truthVals.get(gpIndex);
                int ftcPerValue[] = new int[numPossibleValues];

                for (int predValue = 0; predValue < numPossibleValues; ++predValue) {
                    state.truthVals.set(gpIndex, predValue);
                    boolean isFormulaSatisfied = true;
                    for (GroundClause gc : newFormula.groundClauses) {
                        boolean isClauseSatisfied = false;
                        for (int gpId : gc.groundPredIndices) {
                            BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                            isClauseSatisfied = b.get(state.truthVals.get(gpId));
                            if (isClauseSatisfied) {
                                break;
                            }
                        }
                        isFormulaSatisfied = isClauseSatisfied;
                        if (!isFormulaSatisfied) {
                            break;
                        }
                    }
                    if (isFormulaSatisfied) {

                        ++ftcPerValue[predValue];

                    }
                }
                Map<Integer, List<Integer>> satCountsMap = satCountsPerDomain.get(gpIndex);
                if (!satCountsMap.containsKey(formulaId)) {
                    satCountsMap.put(formulaId, new ArrayList<>(Collections.nCopies(numPossibleValues, 0)));
                }
                for (int j = 0; j < numPossibleValues; j++) {
                    satCountsMap.get(formulaId).set(j, satCountsMap.get(formulaId).get(j) + ftcPerValue[j]);
                }
                state.truthVals.set(gpIndex, tempTrueValue);
            }
        }
    }

    public void learnWeights() {
        time = System.currentTimeMillis();
        setMLNWeights();
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
        LBFGS.Function f = new LBFGS.Function() {
            @Override
            public double evaluate(double[] beta, double[] g, int n, double step) {
                double logLoss = 0;
                Arrays.fill(g, 0);
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
                    logLoss -= getPseudoLogLikelihood(prior);
                    double []tempGrad = getPseudoGradient(prior);
                    for (int i = 0; i < g.length; i++) {
                        g[i] -= tempGrad[i];
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
            p.epsilon = 1.0E-4D;
//            p.past = 5;
//            p.delta = 1.0E-7D;

        LBFGS.ProgressCallback cb = new LBFGS.ProgressCallback() {
            @Override
            public int apply(double[] x, double[] g, double fx, double xnorm,
                             double gnorm, double step, int n, int k, LBFGS.Status ls) {
                if(k % 10 == 0) {
                    U.pf("ITER %d obj=%g\n", k, fx);
                    System.out.println("Weights: "+Arrays.toString(x));
                    System.out.println("Gradient: "+Arrays.toString(g));
                    System.out.println("Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
                }
                return 0;
            }
        };
//            weights = new double[0];
        LBFGS.Result r = LBFGS.lbfgs(weights, f, cb, p);
        U.p(r);
    }

    private double[] getPseudoGradient(double[] prior) {
        double[] gradientPerDomain;
        Arrays.fill(gradient,0.0);
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            gradientPerDomain = getPseudoGradient(domainID);
            for (int formulaID = 0; formulaID < numWts; ++formulaID){
                gradient[formulaID] += gradientPerDomain[formulaID];
            }
        }

        for (int formulaID = 0; formulaID < numWts; ++formulaID){
            gradient[formulaID] -= (weights[formulaID] - priorMeans[formulaID]) * prior[formulaID];
        }
        return gradient;
    }

    private double[] getPseudoGradient(int domainID) {
        State state = states.get(domainID);
        List<Map<Integer, List<Integer>>> satCountsPerDomain = satCounts.get(domainID);
        double[] gradientLocal = new double[numWts];
        int numGroundPreds = state.groundMLN.groundPredicates.size();
        for (int predId = 0; predId < numGroundPreds; predId++) {
            int trueVal = state.truthVals.get(predId);
            int numPossibleVals = state.groundMLN.groundPredicates.get(predId).numPossibleValues;
            List<Double> probDistribution = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(predId));
            // For all the formulas in which this pred occurs, update their gradients
            Map<Integer, List<Integer>> satCountsPerDomainPerPred = satCountsPerDomain.get(predId);
            double[] gradientPerGp = new double[numWts];
            String predName = states.get(domainID).groundMLN.groundPredicates.get(predId).symbol.symbol;
            for(int formulaId : satCountsPerDomainPerPred.keySet())
            {
                // for ith formula, add n_i@(predId = trueVal)
                gradientPerGp[formulaId] += satCountsPerDomainPerPred.get(formulaId).get(trueVal);

                // Now subtract expected value of n_i(predId).
                // To calculate expected value, for each possible value val of predId, sum n_i@val * P(predId = val)
                double E_N_i = 0;
                for(int val = 0 ; val < numPossibleVals ; val++)
                {
                    E_N_i += satCountsPerDomainPerPred.get(formulaId).get(val) * probDistribution.get(val);
                }
                gradientPerGp[formulaId] -= E_N_i;
                gradientPerGp[formulaId] /= predToNumGndingsMap.get(domainID).get(predName);
            }
            for (int formulaId = 0; formulaId < numFormulas; formulaId++) {
                gradientLocal[formulaId] += gradientPerGp[formulaId];
            }
            if(priorSoftEvidence)
            {
                if(softEvidencePerPredPerVal.get(domainID).containsKey(predId))
                {
                    gradientLocal[numWts-1] += softEvidencePerPredPerVal.get(domainID).get(predId).get(trueVal);
                    double expectedVal = 0.0;
                    for (int val = 0; val < numPossibleVals; val++) {
                        expectedVal += softEvidencePerPredPerVal.get(domainID).get(predId).get(val);
                    }
                    gradientLocal[numWts-1] -= expectedVal;
                    gradientLocal[numWts-1] /= predToNumGndingsMap.get(domainID).get(predName);
                }

            }
        }

//        if(priorSoftEvidence){
//            for (Integer predID : lambdaPerPredPerValue.keySet()){
//                int numPossibleValues = predList.get(predID).numPossibleValues;
//                gradient[numFormulas] += lambdaPerPredPerValue.get(predID).get(state.truthVals.get(predID)).getValue();
//                List<Double> probDistribution = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(predID));
//                for (int value = 0; value < numPossibleValues; ++value){
//                    gradient[numFormulas] -= probDistribution.get(value) * lambdaPerPredPerValue.get(predID).get(value).getValue();
//                }
//            }
//        }
        return gradientLocal;
    }

    private double getPseudoLogLikelihood(double[] prior) {
        double likelihood = 0;
        for (int domainIndex = 0; domainIndex < domain_cnt; ++domainIndex){
            likelihood += getPseudoLogLikelihood(domainIndex);
        }

        for (int formulaID = 0; formulaID < numWts; ++formulaID){
            double temp = weights[formulaID] - priorMeans[formulaID];
            likelihood -= temp * temp * prior[formulaID] / 2;
        }
        return likelihood;
    }

    private double getPseudoLogLikelihood(int domainIndex) {
        if(genLearnDebug) {
            System.out.println("Calculating Pseudo Log Likelihood...");
        }
        resetWtsPerPredPerVal(domainIndex);
        setWtsPerPredPerVal(domainIndex);
        State state = states.get(domainIndex);
        int numGps = state.groundMLN.groundPredicates.size();
        double pll = 0;
        for(int gpId = 0 ; gpId < numGps ; gpId++) {
            double d = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(gpId)).get(state.truthVals.get(gpId));
            String predName = state.groundMLN.groundPredicates.get(gpId).symbol.symbol;
            pll += (Math.log(d)/predToNumGndingsMap.get(domainIndex).get(predName));
        }

        return pll;
    }

    private void setWtsPerPredPerVal(int domainIndex) {
        List<Map<Integer, List<Integer>>> satCountsPerDomain = satCounts.get(domainIndex);
        int numGroundPreds = satCountsPerDomain.size();
        MLN mln = mlns.get(0);
        State state = states.get(domainIndex);
        for (int gpId = 0; gpId < numGroundPreds; gpId++) {
            Map<Integer, List<Integer>> satCountsPerDomainPerPred = satCountsPerDomain.get(gpId);
            int numPossibleVals = state.groundMLN.groundPredicates.get(gpId).numPossibleValues;
            for(int formulaId : satCountsPerDomainPerPred.keySet())
            {
                double formulaWt = weights[formulaId];
                List<Integer> satCountsPerDomainPerPredPerFormula = satCountsPerDomainPerPred.get(formulaId);
                for (int val = 0; val < numPossibleVals; val++) {
                    int satCount = satCountsPerDomainPerPredPerFormula.get(val);
                    double oldVal = state.wtsPerPredPerVal.get(gpId).get(val);
                    state.wtsPerPredPerVal.get(gpId).set(val, oldVal + satCount*formulaWt);
                }
            }
            if(priorSoftEvidence)
            {
                if(softEvidencePerPredPerVal.get(domainIndex).containsKey(gpId))
                {
                    double lambda = weights[numWts-1];
                    for (int val = 0; val < numPossibleVals; val++) {
                        double oldVal = state.wtsPerPredPerVal.get(gpId).get(val);
                        state.wtsPerPredPerVal.get(gpId).set(val, oldVal + lambda * softEvidencePerPredPerVal.get(domainIndex).get(gpId).get(val));
                    }
                }
            }

        }
    }

    private void resetWtsPerPredPerVal(int domainIndex) {
        State state = states.get(domainIndex);
        int numGroundPreds = state.groundMLN.groundPredicates.size();
        for (int gpId = 0; gpId < numGroundPreds; gpId++) {
            int numPossibleVals = state.wtsPerPredPerVal.get(gpId).size();
            for(int val = 0 ; val < numPossibleVals ; val++)
            {
                state.wtsPerPredPerVal.get(gpId).set(val,0.0);
            }
        }
    }
}
