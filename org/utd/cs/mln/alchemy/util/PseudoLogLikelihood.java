package org.utd.cs.mln.alchemy.util;

import org.utd.cs.mln.alchemy.core.*;

import java.util.*;

/**
 * Created by Happy on 4/4/18.
 */
public class PseudoLogLikelihood extends Loss{
    // For each domain d, for each ground predicate gp in that d, for each first order formula f in which gp appears,
    // stores number of satisfied groundings of f for each possible value val of gp.
    // satCounts[d][gp][f][val]
    private List<Map<Integer, Map<Integer, List<Integer>>>> satCounts;

    // For each domain d, stores number of ground atoms for each first order predicate p.
    // predToNumGndingsMap[d][p]
    private List<Map<String, Integer>> predToNumGndingsMap;
    private double []gradient;
    List<State> states;
    boolean pllDebug = false;

    // priorLambda is lambda for weight regularization. For each weight, there can be different lambda. If an entry of prior
    // is zero, it means there is no regularization on that corresponding weight.
    // priorMeans and priorStdDevs are the mean and standard deviation vectors for gaussian prior. Note that standard
    // deviation is a vector here, meaning its a diagonal matrix.
    // By default, priorMeans = 0, and priorStdDevs = 1, in that case, gaussian prior is just l2 regularization.
    double [] priorLambda, priorMeans, priorStdDevs;

    public PseudoLogLikelihood(List<State> statesParam){
        init(statesParam);
    }
    //TODO : Implement c'tor for EM also

    /**
     * Initializes fields :
     * <ol>
     *     <li>Create States according to GroundMLNs and truths</li>
     *     <li>Fills in predToNumGndingsMap</li>
     *     <li>Fills in satCounts</li>
     * </ol>
     */
    public void init(List<State> statesParam){
        this.states = statesParam;
        satCounts = new ArrayList<>();
        predToNumGndingsMap = new ArrayList<>();
        int domain_cnt = states.size();
        for (int i = 0; i < domain_cnt; i++) {
            predToNumGndingsMap.add(new HashMap<String,Integer>());
        }
        int numWts = states.get(0).mln.formulas.size();
        if(states.get(0).mln.priorSoftEvidence)
        {
            numWts++;
        }
        gradient = new double[numWts];
        priorMeans = new double[numWts];
        priorStdDevs = new double[numWts];
        Arrays.fill(priorStdDevs,2.0);
        priorLambda = new double[numWts];
        Arrays.fill(priorLambda,0.0);
        createPredToNumGndingsMap();
        long time = System.currentTimeMillis();
        countSatVals();
        System.out.println("Counting took " + org.utd.cs.gm.utility.Timer.time((System.currentTimeMillis() - time)/1000.0));
    }

    /**
     * Fills in {@link PseudoLogLikelihood#predToNumGndingsMap}
     */
    private void createPredToNumGndingsMap() {
        int domain_cnt = states.size();
        for(int domainId = 0 ; domainId < domain_cnt ; domainId++)
        {
            Map<String,Integer> predToNumGndingsMapPerDomain = predToNumGndingsMap.get(domainId);
            GroundMLN groundMln = states.get(domainId).groundMLN;
            for(GroundPredicate gp : groundMln.groundPredToIntegerMap.keySet())
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
        int domain_cnt = states.size();
        // for each domain
        for (int domainId = 0 ; domainId < domain_cnt ; domainId++) {
            System.out.println("Finding counts for domain "+(domainId+1));
            State state = states.get(domainId);
            Map<Integer, Map<Integer, List<Integer>>> satCountsPerDomain = new HashMap<>();
            for (Integer i : state.groundMLN.indexToGroundPredMap.keySet()) {
                satCountsPerDomain.put(i, new HashMap<Integer, List<Integer>>());
            }
            MLN mln = states.get(domainId).mln;
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
                if(pllDebug)
                    System.out.println("Finding counts for formula : " + formulaId);
                countSatValsPerFormula(domainId, formulaId, satCountsPerDomain, new ArrayList<Term>(formulaWiseTermToGround));
            }
            satCounts.add(satCountsPerDomain);
        }
    }


    private void countSatValsPerFormula(int domainId, int formulaId, Map<Integer, Map<Integer, List<Integer>>> satCountsPerDomain, ArrayList<Term> terms) {
        State state = states.get(domainId);
        MLN mln = states.get(domainId).mln;
        int[][] permutations = FullyGrindingMill.permute(terms);
//        if(genLearnDebug)
//            System.out.println("Number of permutations : " + permutations.length);
        Map<GroundPredicate, Integer> groundPredToIntegerMap = state.groundMLN.groundPredToIntegerMap;
        Formula formula = mln.formulas.get(formulaId);
        for (int i = 0; i < permutations.length; i++) {
//            if(genLearnDebug)
//            {
//                if (i % 100000 == 0)
//                    System.out.println("i : " + i);
//            }
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
                int numPossibleValues = state.groundMLN.indexToGroundPredMap.get(gpIndex).numPossibleValues;
                int tempTrueValue = state.truthVals.get(gpIndex);
                int ftcPerValue[] = new int[numPossibleValues];

                for (int predValue = 0; predValue < numPossibleValues; ++predValue) {
                    state.truthVals.put(gpIndex, predValue);
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
                state.truthVals.put(gpIndex, tempTrueValue);
            }
        }
    }

    @Override
    public double getLossValue() {
        double likelihood = 0;
        int domain_cnt = states.size();
        for (int domainIndex = 0; domainIndex < domain_cnt; ++domainIndex){
            likelihood += getPseudoLogLikelihood(domainIndex);
        }
        int numWts = gradient.length;

        double weights[] = new double[numWts];
        MLN firstMLN = states.get(0).mln;
        for (int formulaId = 0; formulaId < firstMLN.formulas.size(); formulaId++) {
            weights[formulaId] = firstMLN.formulas.get(formulaId).weight.getValue();
        }
        if(firstMLN.priorSoftEvidence)
        {
            weights[numWts-1] = firstMLN.softEvidenceLambda;
        }
        for (int wtId = 0; wtId < numWts; ++wtId){
            double temp = (weights[wtId] - priorMeans[wtId])/priorStdDevs[wtId];
            likelihood -= temp * temp * priorLambda[wtId] / 2;
        }
        return -likelihood;
    }

    private double getPseudoLogLikelihood(int domainIndex) {
        if(pllDebug) {
            System.out.println("Calculating Pseudo Log Likelihood...");
        }
        setWtsPerPredPerVal(domainIndex);
        State state = states.get(domainIndex);
        double pll = 0;
        for(Integer gpId : state.groundMLN.indexToGroundPredMap.keySet()) {
            double d = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(gpId)).get(state.truthVals.get(gpId));
            String predName = state.groundMLN.indexToGroundPredMap.get(gpId).symbol.symbol;
            pll += (Math.log(d)/predToNumGndingsMap.get(domainIndex).get(predName));
        }

        return pll;
    }

    private void setWtsPerPredPerVal(int domainIndex) {
        resetWtsPerPredPerVal(domainIndex);
        Map<Integer, Map<Integer, List<Integer>>> satCountsPerDomain = satCounts.get(domainIndex);
        State state = states.get(domainIndex);
        MLN mln = state.mln;
        GroundMLN gMln = state.groundMLN;
        for (Integer gpId : satCountsPerDomain.keySet()) {
            Map<Integer, List<Integer>> satCountsPerDomainPerPred = satCountsPerDomain.get(gpId);
            int numPossibleVals = state.groundMLN.indexToGroundPredMap.get(gpId).numPossibleValues;
            for(int formulaId : satCountsPerDomainPerPred.keySet())
            {
                double formulaWt = mln.formulas.get(formulaId).weight.getValue();
                List<Integer> satCountsPerDomainPerPredPerFormula = satCountsPerDomainPerPred.get(formulaId);
                for (int val = 0; val < numPossibleVals; val++) {
                    int satCount = satCountsPerDomainPerPredPerFormula.get(val);
                    double oldVal = state.wtsPerPredPerVal.get(gpId).get(val);
                    state.wtsPerPredPerVal.get(gpId).set(val, oldVal + satCount*formulaWt);
                }
            }
            if(mln.priorSoftEvidence)
            {
                if(gMln.softEvidencePerPredPerVal.containsKey(gpId))
                {
                    double lambda = mln.softEvidenceLambda;
                    for (int val = 0; val < numPossibleVals; val++) {
                        double oldVal = state.wtsPerPredPerVal.get(gpId).get(val);
                        state.wtsPerPredPerVal.get(gpId).set(val, oldVal + lambda * gMln.softEvidencePerPredPerVal.get(gpId).get(val));
                    }
                }
            }

        }
    }

    public double[] getGradient() {
        double[] gradientPerDomain;
        Arrays.fill(gradient,0.0);
        int domain_cnt = states.size();
        int numWts = gradient.length;
        for (int domainID = 0; domainID < domain_cnt; ++domainID){
            gradientPerDomain = getPseudoGradient(domainID);
            for (int formulaID = 0; formulaID < numWts; ++formulaID){
                gradient[formulaID] += gradientPerDomain[formulaID];
            }
        }

        double weights[] = new double[numWts];
        MLN firstMLN = states.get(0).mln;
        for (int formulaId = 0; formulaId < firstMLN.formulas.size(); formulaId++) {
            weights[formulaId] = firstMLN.formulas.get(formulaId).weight.getValue();
        }
        if(firstMLN.priorSoftEvidence)
        {
            weights[numWts-1] = firstMLN.softEvidenceLambda;
        }
        for (int formulaID = 0; formulaID < numWts; ++formulaID){
            gradient[formulaID] -= (weights[formulaID] - priorMeans[formulaID]) * priorLambda[formulaID] / (priorStdDevs[formulaID]*priorStdDevs[formulaID]);
        }
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] *= -1;
        }
        return gradient;
    }

    private double[] getPseudoGradient(int domainID) {
        State state = states.get(domainID);
        Map<Integer, Map<Integer, List<Integer>>> satCountsPerDomain = satCounts.get(domainID);
        int numWts = gradient.length;
        double[] gradientLocal = new double[numWts];
        for (Integer predId : state.groundMLN.indexToGroundPredMap.keySet()) {
            int trueVal = state.truthVals.get(predId);
            int numPossibleVals = state.groundMLN.indexToGroundPredMap.get(predId).numPossibleValues;
            List<Double> probDistribution = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(predId));
            // For all the formulas in which this pred occurs, update their gradients
            Map<Integer, List<Integer>> satCountsPerDomainPerPred = satCountsPerDomain.get(predId);
            double[] gradientPerGp = new double[numWts];
            String predName = states.get(domainID).groundMLN.indexToGroundPredMap.get(predId).symbol.symbol;
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
            for (int formulaId = 0; formulaId < state.mln.formulas.size(); formulaId++) {
                gradientLocal[formulaId] += gradientPerGp[formulaId];
            }
            if(state.mln.priorSoftEvidence)
            {
                if(state.groundMLN.softEvidencePerPredPerVal.containsKey(predId))
                {
                    gradientLocal[numWts-1] += state.groundMLN.softEvidencePerPredPerVal.get(predId).get(trueVal);
                    double expectedVal = 0.0;
                    for (int val = 0; val < numPossibleVals; val++) {
                        expectedVal += state.groundMLN.softEvidencePerPredPerVal.get(predId).get(val);
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


    private void resetWtsPerPredPerVal(int domainIndex) {
        State state = states.get(domainIndex);
        for (Integer gpId : state.groundMLN.indexToGroundPredMap.keySet()) {
            int numPossibleVals = state.wtsPerPredPerVal.get(gpId).size();
            for(int val = 0 ; val < numPossibleVals ; val++)
            {
                state.wtsPerPredPerVal.get(gpId).set(val,0.0);
            }
        }
    }

}
