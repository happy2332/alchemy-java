//package org.utd.cs.mln.inference;
//
//import org.utd.cs.gm.core.LogDouble;
//import org.utd.cs.gm.utility.Timer;
//import org.utd.cs.mln.alchemy.core.*;
//import org.utd.cs.mln.alchemy.util.OtherUtils;
//
//import java.io.PrintWriter;
//import java.util.*;
//
///**
// * Created by Happy on 3/2/17.
// */
//public class GibbsSampler_v2 {
//
//    public MLN mln;
//    public Evidence truth;
//    public State state;
//    private List<List<Integer>> countNumAssignments = new ArrayList<>(); // For each groundPred in state.mln.groundPreds, stores how many times this groundpred gets assigned to a  particular value. Used for calculating marginal prob
//    private List<List<Double>> marginals = new ArrayList<>();
//    private double[][] allFormulaTrueCnts; // allFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample
//    private double[][] oldAllFormulaTrueCnts; // oldAllFormulaTrueCnts[i][j] is the number of true groundings of jth formula in ith sample in the previous iter of learning
//    public double[] numFormulaTrueCnts, numFormulaTrueSqCnts, allLambdaFormulaTrueCnts, oldAllLambdaFormulaTrueCnts; // sum of true counts of formulas over all samples.
//    public double numLambdaTrueCnts, numLambdaTrueSqCnts;
//    public List<Integer> mlnStatePerIter[]; // For each groundPredicate in mln.groundPredicates, stores its truthval
//
//    public boolean priorSoftEvidence;
//
//    private boolean trackState = false, trackFormulaCounts = false, saveAllCounts = false, gsdebug=true, calculate_marginal=true;
//    private Random rand;
//    public int numBurnSteps, numIter;
//
//    ArrayList<Integer> allGndPredIndices;
//    int[][][] ftcPerPredPerValue;
//    HashMap<Integer, HashMap<Integer, LogDouble>> lambdaPerPredPerValue = new HashMap<>();
//
//    public GibbsSampler_v2(MLN mln, GroundMLN groundMLN, Evidence truth, int numBurnSteps, int numIter, boolean trackFormulaCounts, boolean calculate_marginal, boolean priorSoftEvidence, boolean trackState) {
//        this.mln = mln;
//        this.truth = truth;
//        state = new State(mln, groundMLN);
//        this.numBurnSteps = numBurnSteps;
//        this.numIter = numIter;
//        this.trackFormulaCounts = trackFormulaCounts;
//        this.trackState = trackState;
//        this.calculate_marginal = calculate_marginal;
//        int numGroundPreds = groundMLN.indexToGroundPredMap.size();
//
//        // Initialize all counts of assignments to 0.
//        for(GroundPredicate gp : groundMLN.indexToGroundPredMap.values())
//        {
//            int numValues = gp.numPossibleValues;
//            countNumAssignments.add(new ArrayList<>(Collections.nCopies(numValues,0)));
//        }
//
//        // Intialize all counts of assignments to 0.
//        for(int i = 0 ; i < numGroundPreds ; i++)
//        {
//            int numValues = state.groundMLN.groundPredicates.get(i).numPossibleValues;
//            marginals.add(new ArrayList<>(Collections.nCopies(numValues,0.0)));
//        }
//
//        allGndPredIndices = new ArrayList<>();
//        ftcPerPredPerValue = new int[state.groundMLN.groundPredicates.size()][][];
//
//        if(trackFormulaCounts){
//            numFormulaTrueCnts = new double[mln.formulas.size()];
//            numFormulaTrueSqCnts = new double[mln.formulas.size()];
//            allFormulaTrueCnts = new double[numIter][mln.formulas.size()];
//            allLambdaFormulaTrueCnts = new double[numIter];
//        }
//        if(trackState){
//            mlnStatePerIter = new ArrayList[numIter];
//        }
//        this.priorSoftEvidence = priorSoftEvidence;
////        rand = new Random(LearnTest.getSeed());
//        rand = new Random(123456);
//    }
//
//    private void init() {
//        if(gsdebug)
//            System.out.println("Inference initialization");
//        doInitialRandomAssignment();
//        initializeNumSatValues(); // Initialize numFalseClauseIds and numSatLiterals in the clauses
//        ArrayList<Integer> allGndPredIndices = new ArrayList<>();
//        int numGndPreds = state.groundMLN.groundPredicates.size();
//
//        for(int i = 0 ; i < numGndPreds ; i++) {
//            allGndPredIndices.add(i);
//        }
//
//        updateWtsForGndPreds(allGndPredIndices);
//
//        if(trackFormulaCounts){
//            for (int i = 0; i < numFormulaTrueCnts.length; i++) {
//                numFormulaTrueCnts[i] = 0.0;
//                numFormulaTrueSqCnts[i] = 0.0;
//            }
//        }
//        numLambdaTrueCnts = 0.0;
//        numLambdaTrueSqCnts = 0;
//    }
//
//    //without EM
//    public void updateMLN(boolean doInit) {
//        if(doInit){
//            doInitialTrueAssignment();
//            initializeNumSatValues(); // Initialize numFalseClauseIds and numSatLiterals in the clauses
//            int numGndPreds = state.groundMLN.groundPredicates.size();
//
//            for(int i = 0 ; i < numGndPreds ; i++) {
//                allGndPredIndices.add(i);
//            }
//        }
//
//        updateWtsForGndPreds(allGndPredIndices);
//    }
//
//    //with EM
//    public void updateMLN(boolean init, List<Integer> stTruthVals, Map<Integer, Integer> subtypeMap) {
//        doInitialTrueAssignment(init, stTruthVals, subtypeMap);
//        initializeNumSatValues(); // Initialize numFalseClauseIds and numSatLiterals in the clauses
//        if(init){
//            int numGndPreds = state.groundMLN.groundPredicates.size();
//
//            allGndPredIndices.clear();
//            for(int i = 0 ; i < numGndPreds ; i++) {
//                allGndPredIndices.add(i);
//            }
//        }
//
//        updateWtsForGndPreds(allGndPredIndices);
//    }
//
//    public double getPseudoLogLikelihood(){
//        if(gsdebug) {
//            System.out.println("Calculating Pseudo Log Likelihood...");
//        }
//
//        int numGps = state.groundMLN.groundPredicates.size();
//        double pll = 0;
//        for(int i = 0 ; i < numGps ; i++) {
//            double d = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(i)).get(state.truthVals.get(i));
//            pll += Math.log(d);
//        }
//
//        return pll;
//    }
//
//    public double[] getPseudoGradient(boolean priorSoftEvidence) {
//        int numFormulas = mln.formulas.size();
//        double[] gradient;
//        if(priorSoftEvidence){
//            gradient = new double[numFormulas + 1];
//        } else {
//            gradient = new double[numFormulas];
//        }
////        int x = 0;
//        List<GroundPredicate> predList = state.groundMLN.groundPredicates;
//        int numPred = predList.size();
////        System.out.println("Test Gradient: ");
//        for (int predID = 0; predID < numPred; ++predID) {
////            System.out.println("PredID: "+predID);
////            System.out.println("Truth Value: "+state.truthVals.get(predID));
//            int numPossibleValues = predList.get(predID).numPossibleValues;
//            List<Double> probDistribution = OtherUtils.getProbabilityDistribution(state.wtsPerPredPerVal.get(predID));
////            System.out.println("probDistribution: "+probDistribution);
////            for (int formulaID = x; formulaID < x + 1; ++formulaID){
//            for (int formulaID = 0; formulaID < numFormulas; ++formulaID){
////                double deltaG;
////                System.out.println("formulaID: "+formulaID);
////                System.out.println("ftcPerValue: "+Arrays.toString(ftcPerPredPerValue[predID][formulaID]));
//                gradient[formulaID] += ftcPerPredPerValue[predID][formulaID][state.truthVals.get(predID)];
//                for (int value = 0; value < numPossibleValues; ++value){
//                    double temp = probDistribution.get(value) * ftcPerPredPerValue[predID][formulaID][value];
//                    gradient[formulaID] -= temp;
//                }
//            }
//        }
////        System.out.println("Gradient["+x+"]: " + gradient[x]);
//
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
//        return gradient;
//    }
//
//    public void findFtcPerPredicate(){
//        int numFormulas = mln.formulas.size();
//        List<GroundPredicate> predList = state.groundMLN.groundPredicates;
//        int numPred = predList.size();
//        for (int predID = 0; predID < numPred; ++predID){
//            int numPossibleValues = predList.get(predID).numPossibleValues;
//            int[][] ftcPerValue = new int[numFormulas][numPossibleValues];
//            Map<Integer, Set<Integer>> formulaHavingPred = predList.get(predID).groundFormulaIds;
//
//            for(Integer formulaID : formulaHavingPred.keySet()){
//                GroundFormula gf = state.groundMLN.groundFormulas.get(formulaID);
//                List<Integer> parentFormulaId = gf.parentFormulaId;
//                List<Integer> numCopies = gf.numCopies;
//                if(gf.parentFormulaId.isEmpty()) {
//                    BitSet b = gf.groundClauses.get(0).grounPredBitSet.get(0);
//                    Integer predValue = 0;
//                    for(int i = 0; i < b.length(); ++i) {
//                        if (b.get(i)) {
//                            predValue = i;
//                        }
//                    }
//                    HashMap<Integer, LogDouble> subtypeValueToWeight;
//                    if(lambdaPerPredPerValue.containsKey(predID)){
//                        subtypeValueToWeight = lambdaPerPredPerValue.get(predID);
//                    } else {
//                        subtypeValueToWeight = new HashMap<>();
//                    }
////                    HashMap<Integer, LogDouble> subtypeValueToWeight = new HashMap<>();
//                    subtypeValueToWeight.put(predValue, gf.originalWeight);
//                    lambdaPerPredPerValue.put(predID, subtypeValueToWeight);
//                } else {
//                    int tempTrueValue = state.truthVals.get(predID);
//                    for (int predValue = 0; predValue < numPossibleValues; ++predValue){
//                        state.truthVals.set(predID, predValue);
//                        boolean isFormulaSatisfied = true;
//                        for (GroundClause gc : gf.groundClauses) {
//                            boolean isClauseSatisfied = false;
//                            for (int gpId : gc.groundPredIndices) {
//                                BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
//                                isClauseSatisfied = b.get(state.truthVals.get(gpId));
////                                System.out.println(gc);
//                                if (isClauseSatisfied) {
//                                    break;
//                                }
//                            }
//                            isFormulaSatisfied = isClauseSatisfied;
//                            if (!isFormulaSatisfied) {
//                                break;
//                            }
//                        }
//                        if (isFormulaSatisfied) {
//                            for (int i = 0; i < parentFormulaId.size(); i++) {
//                                ftcPerValue[parentFormulaId.get(i)][predValue] += numCopies.get(i) ;
//                            }
//
//                        }
////                        System.out.println(gf);
//                    }
//                    state.truthVals.set(predID, tempTrueValue);
//                }
//            }
////            int assignment = 0;
////            if(truth.predIdVal.containsKey(predID))
////                assignment = truth.predIdVal.get(predID);
////            else(){
////
////            }
////            state.truthVals.set(predID, assignment);
//            ftcPerPredPerValue[predID] = ftcPerValue;
//        }
//    }
//
////    without EM
//    private void doInitialTrueAssignment() {
//        GroundMLN gmln = state.groundMLN;
//        List<GroundPredicate> gpList = gmln.groundPredicates;
//        int numGps = gpList.size();
//
//        for(int i = 0 ; i < numGps ; i++) {
//            int assignment = 0;
//            if(truth.predIdVal.containsKey(i))
//                assignment = truth.predIdVal.get(i);
//            state.truthVals.set(i, assignment);
//        }
//    }
//
////    with EM
//    private void doInitialTrueAssignment(boolean init, List<Integer> stTruthVals, Map<Integer, Integer> subtypeMap) {
//        GroundMLN gmln = state.groundMLN;
//        List<GroundPredicate> gpList = gmln.groundPredicates;
//        int numGps = gpList.size();
//        int assignment;
//        if(init){
//            for(int i = 0 ; i < numGps ; i++) {
//                assignment = 0;
//                if(truth.predIdVal.containsKey(i))
//                    assignment = truth.predIdVal.get(i);
//                state.truthVals.set(i, assignment);
//            }
//        }
//        for (Integer predID : subtypeMap.keySet()){
//            if(subtypeMap.get(predID) == -1){
//                assignment = getUniformAssignment(gpList.get(predID).numPossibleValues);
//            } else {
//                assignment = stTruthVals.get(subtypeMap.get(predID));
//            }
//            state.truthVals.set(predID, assignment);
//        }
//
//    }
//
//    private void doInitialRandomAssignment() {
//        GroundMLN gmln = state.groundMLN;
//        List<GroundPredicate> gpList = gmln.groundPredicates;
//        int numGps = gpList.size();
//
//        for(int i = 0 ; i < numGps ; i++)
//        {
//            int numPossibleVals = gpList.get(i).numPossibleValues;
//            int assignment = getUniformAssignment(numPossibleVals);
//            state.truthVals.set(i,assignment);
//        }
//    }
//
//    // According to present state, initialize state.numFalseClausesSet and state.numTrueLiterals
//    private void initializeNumSatValues() {
//        GroundMLN gm = state.groundMLN;
//        int numFormulas = state.groundMLN.groundFormulas.size();
//        for(int i = 0 ; i < numFormulas ; i++) {
//            GroundFormula gf = gm.groundFormulas.get(i);
//            Set<Integer> falseClauseIds = new HashSet<>();
//            int numClauses = gf.groundClauses.size();
//
//            for(int j = 0 ; j < numClauses ; j++) {
//                GroundClause gc = gf.groundClauses.get(j);
//                int numPreds = gc.groundPredIndices.size();
//                //int numAtoms = gc.groundAtoms.size();
//                int numSatLiterals = 0;
//                for(int k = 0 ; k < numPreds ; k++) {
//                    int globalPredIndex = gc.groundPredIndices.get(k);
//                    int currentAssignment = state.truthVals.get(globalPredIndex);
//                    if (gc.grounPredBitSet.get(k).get(currentAssignment)) {
//                        numSatLiterals++;
//                    }
//                }
//                state.numTrueLiterals.get(i).set(j, numSatLiterals);
//
//                if (numSatLiterals == 0) {
//                    falseClauseIds.add(j);
//                }
//            } // end of clause
//            state.falseClausesSet.set(i,falseClauseIds);
//        }// end of formula
//    }
//
//    public void infer(boolean burningIn, boolean isInit) {
//        if(isInit)
//            init();
//
//        int numGndPreds = state.groundMLN.groundPredicates.size();
//
//        if(burningIn)
//        {
//            // Burning in
//
//            System.out.println("Burning in started...");
//            long time = System.currentTimeMillis();
//
//            for(int i =1 ; i <= numBurnSteps; i++)
//            {
//                for(int gpId =0; gpId < numGndPreds; gpId++){
//                    performGibbsStep(gpId, (gpId+1)%numGndPreds);
//                }
//                if(i%100 == 0) {
//                    System.out.println("iter : " + i + ", Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
//                }
//            }
//
//            System.out.println("Burning completed in : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
//
//        }
//
//        System.out.println("Gibbs sampling started...");
//        long time = System.currentTimeMillis();
//        for(int i = 1 ; i <= numIter; i++) {
//            for(int gpId = 0; gpId < numGndPreds; gpId++){
//                int assignment = performGibbsStep(gpId, (gpId+1)%numGndPreds);
//                if(calculate_marginal)
//                    countNumAssignments.get(gpId).set(assignment, countNumAssignments.get(gpId).get(assignment)+1);
//            }
//
//            if(trackFormulaCounts) {
//                int numWts = mln.formulas.size();
//                double []numTrueGndings = state.getNumTrueGndings(numWts);
//                for (int j = 0; j < numWts; j++) {
//                    numFormulaTrueCnts[j] += numTrueGndings[j];
//                    allFormulaTrueCnts[i-1][j] += numTrueGndings[j]; //i-1 because iterations start from 1
//                    numFormulaTrueSqCnts[j] += numTrueGndings[j]*numTrueGndings[j];
//                }
//                if(priorSoftEvidence){
//                    numLambdaTrueCnts += numTrueGndings[numWts];
//                    numLambdaTrueSqCnts += numTrueGndings[numWts] * numTrueGndings[numWts];
//                    allLambdaFormulaTrueCnts[i - 1] = numTrueGndings[numWts];
//                }
//            }
//            if(trackState){
//                mlnStatePerIter[i - 1] = new ArrayList<>(state.truthVals);
//            }
//            if(i%100 == 0) {
//                System.out.println("iter : " + i + ", Elapsed Time : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
//            }
//        }
//
//        if(calculate_marginal)
//        {
//            for(int i = 0 ; i < numGndPreds ; i++)
//            {
//                double sum = 0.0;
//                int numValues = countNumAssignments.get(i).size();
//                for(int j = 0 ; j < numValues ; j++)
//                {
//                    sum += countNumAssignments.get(i).get(j);
//                }
//                for(int j = 0 ; j < numValues ; j++)
//                {
//                    double marg = countNumAssignments.get(i).get(j)/sum;
//                    marginals.get(i).set(j,marg);
//                }
//            }
//        }
//
//        System.out.println("Gibbs Sampling completed in : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
//    }
//
//    public void writeMarginal(PrintWriter writer) {
//
//        for(int i = 0 ; i < marginals.size() ; i++)
//        {
//            for(int j = 0 ; j < marginals.get(i).size() ; j++)
//            {
//                double marginal = marginals.get(i).get(j);
//                writer.println(state.groundMLN.groundPredicates.get(i) + "=" + j + " " + marginal);
//            }
//        }
//
//    }
//
//    private int performGibbsStep(int gpId, int nextGpId) {
//        int assignment = get_probabilistic_assignment(state.wtsPerPredPerVal.get(gpId));
//        int prev_assignment = state.truthVals.get(gpId);
//        state.truthVals.set(gpId, assignment);
//        //List<Integer> affectedGndPredIndices = new ArrayList<>();
//        if(assignment != prev_assignment)
//        {
//             // Markov Blanket for current flipped atom. When flipping an atom, if the value changes then we need to update satWeights for all these M.B predicates.
//            //findMarkovBlanket(gpId, assignment, prev_assignment, affectedGndPredIndices);
//            updateSatCounts(gpId, assignment, prev_assignment);
//            //updateWtsForGndPreds(affectedGndPredIndices);
//        }
//
////        affectedGndPredIndices.clear();
////        affectedGndPredIndices.add(nextGpId);
////        updateWtsForGndPreds(affectedGndPredIndices);
//        updateWtsForNextGndPred(nextGpId);
//
//        return assignment;
//    }
//
//    private void updateSatCounts(int gpId, int assignment, int prev_assignment) {
//        GroundPredicate gp = state.groundMLN.groundPredicates.get(gpId);
//        for(int formulaId : gp.groundFormulaIds.keySet())
//        {
//            for(int cid : gp.groundFormulaIds.get(formulaId))
//            {
//                GroundClause gc = state.groundMLN.groundFormulas.get(formulaId).groundClauses.get(cid);
//                int localPredindex = gc.globalToLocalPredIndex.get(gpId);
//                BitSet b = gc.grounPredBitSet.get(localPredindex);
//                int cur_val = b.get(assignment) ? 1 : 0;
//                int prev_val = b.get(prev_assignment) ? 1 : 0;
//                int satDifference = cur_val - prev_val; // numsatLiterals according to new assignment - old assignment
//                int curSatLiterals = state.numTrueLiterals.get(formulaId).get(cid);
//                curSatLiterals += satDifference;
//                state.numTrueLiterals.get(formulaId).set(cid, curSatLiterals);
//
//                if(curSatLiterals > 0)
//                {
//                    state.falseClausesSet.get(formulaId).remove(cid);
//                }
//                else
//                {
//                    state.falseClausesSet.get(formulaId).add(cid);
//                }
//            }
//        }
//    }
//
//    private void findMarkovBlanket(int gpId, int assignment, int prev_assignment, List<Integer> markov_blanket) {
//        //Set<Integer> mbSet = new HashSet<>();
//        GroundPredicate gp = state.groundMLN.groundPredicates.get(gpId);
//        for(int formulaId : gp.groundFormulaIds.keySet())
//        {
//            for(int cid : gp.groundFormulaIds.get(formulaId))
//            {
//                GroundClause gc = state.groundMLN.groundFormulas.get(formulaId).groundClauses.get(cid);
//                int localPredindex = gc.globalToLocalPredIndex.get(gpId);
//                BitSet b = gc.grounPredBitSet.get(localPredindex);
//                int cur_val = b.get(assignment) ? 1 : 0;
//                int prev_val = b.get(prev_assignment) ? 1 : 0;
//                int satDifference = cur_val - prev_val; // numsatLiterals according to new assignment - old assignment
//                int curSatLiterals = state.numTrueLiterals.get(formulaId).get(cid);
//                curSatLiterals += satDifference;
//                state.numTrueLiterals.get(formulaId).set(cid, curSatLiterals);
//
//                /*
//                for(int otherGpId : gc.groundPredIndices)
//                {
//                    if(otherGpId != gpId)
//                    {
//                        mbSet.add(otherGpId);
//                    }
//                }*/
//                if(curSatLiterals > 0)
//                {
//                    state.falseClausesSet.get(formulaId).remove(cid);
//                }
//                else
//                {
//                    state.falseClausesSet.get(formulaId).add(cid);
//                }
//            }
//
////            GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
////            mbSet.addAll(gf.groundPredIndices);
////            mbSet.remove(gpId);
//        }
////        markov_blanket.addAll(mbSet);
//    }
//
//    public void updateWtsForNextGndPred(int i) {
//
//        List<GroundPredicate> groundPredicates = state.groundMLN.groundPredicates;
//        GroundPredicate gp = groundPredicates.get(i);
//        int numPossibleVals =gp.numPossibleValues;
//        double wtsPerVal[] = new double[numPossibleVals];
//        Map<Integer, Set<Integer>> formulaIds = gp.groundFormulaIds;
//
//        for(Integer formulaId : formulaIds.keySet())
//        {
//            GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
//            double wt = gf.weight.getValue();
//            Set<Integer> tempSet = new HashSet<Integer>();
//            tempSet.addAll(state.falseClausesSet.get(formulaId));
//            tempSet.removeAll(formulaIds.get(formulaId));
//            BitSet formulaBitSet = new BitSet(numPossibleVals);
//            formulaBitSet.flip(0,numPossibleVals);
//
//            // If there is a clause which is false, and doesn't contain gp, then this formula is always false
//            if(tempSet.size() == 0)
//            {
//
//                for(Integer cid : formulaIds.get(formulaId))
//                {
//                    BitSet clauseBitSet = new BitSet(numPossibleVals);
//                    GroundClause gc = gf.groundClauses.get(cid);
//                    if(!gc.globalToLocalPredIndex.containsKey(i))
//                    {
//                        System.out.println("problem....");
//                    }
//                    int localPredIndex = gc.globalToLocalPredIndex.get(i);
//                    int numSatLiterals = state.numTrueLiterals.get(formulaId).get(cid);
//                    if(numSatLiterals > 1)
//                        clauseBitSet.flip(0,numPossibleVals);
//                    else if(numSatLiterals == 1)
//                    {
//                        BitSet b = gc.grounPredBitSet.get(localPredIndex);
//                        if(b.get(state.truthVals.get(i)))
//                        {
//                            clauseBitSet = (BitSet) b.clone();
//                        }
//                        else
//                        {
//                            clauseBitSet.flip(0,numPossibleVals);
//                        }
//                    }
//                    else {
//                        BitSet b = gc.grounPredBitSet.get(localPredIndex);
//                        clauseBitSet = (BitSet) b.clone();
//                    }
//                    formulaBitSet.and(clauseBitSet);
//                }// end clauses loop
//
//                int startIndex = 0;
//                while(startIndex < numPossibleVals)
//                {
//                    int index = formulaBitSet.nextSetBit(startIndex);
//                    if(index == -1)
//                        break;
//                    wtsPerVal[index] += wt;
//                    startIndex = index+1;
//                }
//            }// end if condition
//
//        } //end formulas loops
//
//        List<Double> tempWts = new ArrayList<>();
//        for(double d : wtsPerVal)
//        {
//            tempWts.add(d);
//        }
//        state.wtsPerPredPerVal.set(i, tempWts);
//    }
//
//    public void updateWtsForGndPreds(List<Integer> affectedGndPredIndices) {
//        List<GroundPredicate> groundPredicates = state.groundMLN.groundPredicates;
//        for(Integer i : affectedGndPredIndices)
//        {
//            GroundPredicate gp = groundPredicates.get(i);
//            int numPossibleVals =gp.numPossibleValues;
//            double wtsPerVal[] = new double[numPossibleVals];
//            Map<Integer, Set<Integer>> formulaIds = gp.groundFormulaIds;
//
//            for(Integer formulaId : formulaIds.keySet())
//            {
//                GroundFormula gf = state.groundMLN.groundFormulas.get(formulaId);
//                double wt = gf.weight.getValue();
//                Set<Integer> tempSet = new HashSet<Integer>();
//                tempSet.addAll(state.falseClausesSet.get(formulaId));
//                tempSet.removeAll(formulaIds.get(formulaId));
//                BitSet formulaBitSet = new BitSet(numPossibleVals);
//                formulaBitSet.flip(0,numPossibleVals);
//
//                // If there is a clause which is false, and not doesn't contain gp, then this formula is always false
//                if(tempSet.size() == 0)
//                {
//                    for(Integer cid : formulaIds.get(formulaId))
//                    {
//                        BitSet clauseBitSet = new BitSet(numPossibleVals);
//                        GroundClause gc = gf.groundClauses.get(cid);
//                        int localPredIndex = gc.globalToLocalPredIndex.get(i);
//                        int numSatLiterals = state.numTrueLiterals.get(formulaId).get(cid);
//                        if(numSatLiterals > 1)
//                            clauseBitSet.flip(0,numPossibleVals);
//                        else if(numSatLiterals == 1)
//                        {
//                            BitSet b = gc.grounPredBitSet.get(localPredIndex);
////                            try {
//                                if (b.get(state.truthVals.get(i))) {
//                                    clauseBitSet = (BitSet) b.clone();
//                                } else {
//                                    clauseBitSet.flip(0, numPossibleVals);
//                                }
////                            } catch (NullPointerException npe){
////                                System.out.println("Global PredID: "+i);
////                                System.out.println("Local PredID: "+localPredIndex);
////                                System.out.println("BitSet: "+b);
////                                System.out.println("Symbol: "+groundPredicates.get(i).symbol.symbol+", Value: "+groundPredicates.get(i).terms);
////                                System.out.println("State: "+state);
//////                                for (Integer key : state.truthVals.)
////                                System.out.println("State.truthVals: "+state.truthVals);
////                                npe.printStackTrace();
////                                System.exit(0);
////                            }
//                        }
//                        else {
//                            BitSet b = gc.grounPredBitSet.get(localPredIndex);
//                            clauseBitSet = (BitSet) b.clone();
//                        }
//                        formulaBitSet.and(clauseBitSet);
//                    }// end clauses loop
//
//                    int startIndex = 0;
//                    while(startIndex < numPossibleVals)
//                    {
//                        int index = formulaBitSet.nextSetBit(startIndex);
//                        if(index == -1)
//                            break;
//                        wtsPerVal[index] += wt;
//                        startIndex = index+1;
//                    }
//                }// end if condition
//
//            } //end formulas loops
//
//            List<Double> tempWts = new ArrayList<>();
//            for(double d : wtsPerVal)
//            {
//                tempWts.add(d);
//            }
//            state.wtsPerPredPerVal.set(i, tempWts);
//        }
//    }
//
//
//
//    private int get_probabilistic_assignment(List<Double> satWeights) {
//        List<Double> probabilities = OtherUtils.getProbabilityDistribution(satWeights);
//        return getRandomAssignment(probabilities);
//    }
//
//    private int getRandomAssignment(List<Double> probabilities) {
////        double p = Math.random();
////        Random rand = new Random(LearnTest.getSeed());
////        Random rand = new Random();
//        double p = rand.nextDouble();
//        double cumulativeProbability = 0.0;
//        for (int i = 0 ; i < probabilities.size() ; i++) {
//            cumulativeProbability += probabilities.get(i);
//            if (p <= cumulativeProbability) {
//                return i;
//            }
//        }
//        return probabilities.size()-1;
//    }
//
//    private int getUniformAssignment(int numPossibleVals) {
////        Random rand = new Random(LearnTest.getSeed());
////        Random rand = new Random();
//        return rand.nextInt(numPossibleVals);
//    }
//
//    public double[] getHessianVectorProduct(double[] v) {
//        int numFormulas = mln.formulas.size();
//        int numWts = numFormulas;
//        if (priorSoftEvidence)
//            numWts = numFormulas + 1;
//        int numSamples = numIter;
//
//        // For minimizing the negative log likelihood,
//        // the ith element of H v is:
//        //   E[n_i * vn] - E[n_i] E[vn]
//        // where n is the vector of all clause counts
//        // and vn is the dot product of v and n.
//
//        double sumVN = 0;
//        double []sumN = new double[numWts];
//        double []sumNiVN = new double[numWts];
//
//        // Get sufficient statistics from each sample,
//        // so we can compute expectations
//        for (int s = 0; s < numSamples; s++)
//        {
//            double[] n1 = allFormulaTrueCnts[s];
//            double n2 = allLambdaFormulaTrueCnts[s];
//
//            // Compute v * n
//            double vn = 0;
//
//            for (int i = 0; i < numFormulas; i++)
//                vn += v[i] * n1[i];
//
//            if(priorSoftEvidence)
//                vn += v[numFormulas] * n2;
//
//            // Tally v*n, n_i, and n_i v*n
//            sumVN += vn;
//            for (int i = 0; i < numFormulas; i++)
//            {
//                sumN[i]    += n1[i];
//                sumNiVN[i] += n1[i] * vn;
//            }
//            if(priorSoftEvidence){
//                sumN[numFormulas] += n2;
//                sumNiVN[numFormulas] += n2 * vn;
//            }
//        }
//
//        // Compute actual product from the sufficient stats
//        double []product = new double[numWts];
//        for (int formulano = 0; formulano < numFormulas; formulano++)
//        {
//            double E_vn = sumVN/numSamples;
//            double E_ni = sumN[formulano]/numSamples;
//            double E_nivn = sumNiVN[formulano]/numSamples;
//            product[formulano] = E_nivn - E_ni * E_vn;
//        }
//        if(priorSoftEvidence){
//            double E_vn = sumVN/numSamples;
//            double E_ni = sumN[numFormulas]/numSamples;
//            double E_nivn = sumNiVN[numFormulas]/numSamples;
//            product[numFormulas] = E_nivn - E_ni * E_vn;
//        }
//
//        return product;
//    }
//
//    public void resetCnts() {
//        numLambdaTrueCnts = 0.0;
//        numLambdaTrueSqCnts = 0.0;
//        if(trackFormulaCounts){
//            Arrays.fill(numFormulaTrueCnts,0.0);
//            Arrays.fill(numFormulaTrueSqCnts,0.0);
//            allFormulaTrueCnts = new double[numIter][mln.formulas.size()];
//            allLambdaFormulaTrueCnts = new double[numIter];
//        }
//    }
//
//    public void saveAllCounts(boolean saveCounts) {
//        if(saveAllCounts == saveCounts)
//            return;
//
//        saveAllCounts = saveCounts;
//        if(trackFormulaCounts){
//            allFormulaTrueCnts = new double[numIter][mln.formulas.size()];
//            oldAllFormulaTrueCnts = new double[numIter][mln.formulas.size()];
//            allLambdaFormulaTrueCnts = new double[numIter];
//            oldAllLambdaFormulaTrueCnts = new double[numIter];
//        }
//    }
//
//    public void restoreCnts() {
//        if (!saveAllCounts)
//            return;
//
//        resetCnts();
//        System.out.println("Restoring counts...");
//        if(trackFormulaCounts){
//            for (int i = 0; i < oldAllFormulaTrueCnts.length; i++)
//            {
//                int numcounts = oldAllFormulaTrueCnts[i].length;
//                for (int j = 0; j < numcounts; j++)
//                {
//                    allFormulaTrueCnts[i][j] = oldAllFormulaTrueCnts[i][j];
//                    numFormulaTrueCnts[j] += allFormulaTrueCnts[i][j];
//                    numFormulaTrueSqCnts[j] += allFormulaTrueCnts[i][j] * allFormulaTrueCnts[i][j];
//                }
//                allLambdaFormulaTrueCnts[i] = oldAllLambdaFormulaTrueCnts[i];
//                numLambdaTrueCnts += allLambdaFormulaTrueCnts[i];
//                numLambdaTrueSqCnts += allLambdaFormulaTrueCnts[i] * allLambdaFormulaTrueCnts[i];
//            }
//        }
//    }
//
//    public void saveCnts() {
//        if (!saveAllCounts)
//            return;
//
//        if(trackFormulaCounts){
//            for (int i = 0; i < allFormulaTrueCnts.length; i++)
//            {
//                int numcounts = allFormulaTrueCnts[i].length;
//                for (int j = 0; j < numcounts; j++)
//                {
//                    oldAllFormulaTrueCnts[i][j] = allFormulaTrueCnts[i][j];
//                }
//                oldAllLambdaFormulaTrueCnts[i] = allLambdaFormulaTrueCnts[i];
//            }
//        }
//    }
//
//    public void debugFTC(int formulaID) {
//        for (int predID = 0; predID < ftcPerPredPerValue.length; ++predID){
//            for (int predVal = 0; predVal < ftcPerPredPerValue[predID][formulaID].length; ++predVal){
//                if(ftcPerPredPerValue[predID][formulaID][predVal] != 0){
//                    System.out.println(state.groundMLN.groundPredicates.get(predID) + "=" + predVal);
//                    if(true){
//
//                    }
//                }
//            }
//        }
//    }
//}