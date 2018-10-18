package org.utd.cs.mln.learning;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
import org.utd.cs.mln.alchemy.util.MyAssert;
import org.utd.cs.mln.alchemy.util.Parser;
import org.utd.cs.mln.inference.GibbsSampler_v2;
import org.utd.cs.mln.lmap.MlnToHyperCube;

import java.io.FileNotFoundException;
import java.text.CollationElementIterator;
import java.util.*;

/**
 * @author Happy
 * @since 3/11/17
 *
 * This class implements weight learning. It parses all arguments, creates required MLNs and groundMLNs,
 * and then calls appropriate weight learner.
 * @see WeightLearner
 */
public class LearnTest {

    public static List<Double> firstOrderNumConnections = new ArrayList<>();

    public static void main(String []args) throws FileNotFoundException, CloneNotSupportedException, PredicateNotFound {
        long totaltime = System.currentTimeMillis();
        System.out.println("Parsing command line arguments...");
        LearnArgs lArgs = parseArgs(args);
        int numDb = lArgs.truthFiles.size();

        // Create required data structures

        // List of MLNs of size numDb. All MLNs are same except for domains for predicates.
        List<MLN> mlns = new ArrayList<>();

        // List of truth evidence of size numDb.
        List<Evidence> truths = new ArrayList<>();
        // List of truth evidence of size numDb.
        List<Evidence> truthsEStep = new ArrayList<>();
        // List of truth evidence of size numDb.
        List<Evidence> truthsMStep = new ArrayList<>();

        // List of groundMlns of size numDb.
        List<GroundMLN> groundMlns = new ArrayList<>();
        // List of groundMlns of size numDb.
        List<GroundMLN> groundMlnsEStep = new ArrayList<>();
        // List of groundMlns of size numDb.
        List<GroundMLN> groundMlnsMStep = new ArrayList<>();

        Set<String> evidPreds = null, queryPreds = null, hiddenPreds = null, queryPredsEStep = null, evidPredsEStep = null, queryPredsMStep = null, evidPredsMStep = null;
        List<Map<Integer,Integer>> groundHiddenPredMapMToEStep = new ArrayList<>();
        WeightLearner wl;
        for(int i = 0 ; i < numDb ; i++) {
            FullyGrindingMill fgm = new FullyGrindingMill();
            MLN mln = new MLN(lArgs.priorSoftEvidence);
            mlns.add(mln);
            Parser parser = new Parser(mln);
            System.out.println("Reading DB file "+(i+1));
            List<String> files = new ArrayList<>();
            files.add(lArgs.truthFiles.get(i));
            if(lArgs.evidenceFiles != null)
                files.add(lArgs.evidenceFiles.get(i));

            // Set parser's truthEvidFiles. This is used for setting varTypeToDomain structure.
            parser.setTruthEvidFiles(files);
            parser.parseInputMLNFile(lArgs.mlnFile);
            // Need to set query and evidence predicates only once, they will be same for all DBs
            if(i == 0)
            {
                evidPreds = parser.evidPreds;
                queryPreds = parser.queryPreds;
                hiddenPreds = parser.hiddenPreds;
                if(hiddenPreds.size() > 0 && !lArgs.withEM)
                {
                    throw new ParameterException("hidden Preds found in MLN : " + hiddenPreds + ", provide withEM flag!!!");
                }
                else if(lArgs.withEM && hiddenPreds.isEmpty())
                {
                    throw new ParameterException("Must provide hidden preds withEM!!!");
                }
                firstOrderNumConnections = new ArrayList<>(Collections.nCopies(mln.formulas.size(),0.0));
            }

            List<FirstEvidence> evidList = parser.parseInputEvidenceFile(lArgs.evidenceFiles.get(i));
            List<FirstEvidence> truthList = parser.parseInputEvidenceFile(lArgs.truthFiles.get(i));
            MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
            int origNumFormulas = mln.formulas.size();
            if(lArgs.withEM)
            {
                FullyGrindingMill fgmEStep = new FullyGrindingMill();
                FullyGrindingMill fgmMStep = new FullyGrindingMill();
                queryPredsEStep = hiddenPreds;
                evidPredsEStep = new HashSet<>();
                evidPredsEStep.addAll(queryPreds);
                evidPredsEStep.addAll(evidPreds);
                List<FirstEvidence> evidListEStep = new ArrayList<>();
                evidListEStep.addAll(evidList);
                evidListEStep.addAll(truthList);
                List<FirstEvidence> truthListEStep = new ArrayList<>();
                HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMapEStep = mlnToHyperCube.createPredsHyperCube(evidListEStep, evidPredsEStep, truthListEStep, queryPredsEStep, mln, lArgs.queryEvidence);
                // Create a copy of this mln without formulas. In that copy, we will add hypercube formulas.
                MLN hyperCubeMLNEStep = mln.copyMLNWithoutformulas();
                for(int formulaId = 0 ; formulaId < origNumFormulas ; formulaId++){
                    hyperCubeMLNEStep.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMapEStep, false));
                    if((formulaId+1)%10 == 0)
                        System.out.println(formulaId+1 + " formulas done...");
                }
                // This is a reverse mapping from groundPredicates to their indices in global groundPredicates list. This
                // mapping will be filled by following function call : fgm.createGroundPredicatesFromPredsHashMap
                Map<GroundPredicate, Integer> groundPredicateIntegerMapEStep = new HashMap<>();
                List<GroundPredicate> groundPredicatesEStep = fgmEStep.createGroundPredicatesFromPredsHashMap(predsHyperCubeHashMapEStep, groundPredicateIntegerMapEStep);

                queryPredsMStep = new HashSet<>();
                queryPredsMStep.addAll(queryPreds);
                queryPredsMStep.addAll(hiddenPreds);
                evidPredsMStep = evidPreds;
                List<FirstEvidence> truthListMStep = truthList;
                List<FirstEvidence> evidListMStep = evidList;
                HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMapMStep = mlnToHyperCube.createPredsHyperCube(evidListMStep, evidPredsMStep, truthListMStep, queryPredsMStep, mln, lArgs.queryEvidence);
                // Create a copy of this mln without formulas. In that copy, we will add hypercube formulas.
                MLN hyperCubeMLNMStep = mln.copyMLNWithoutformulas();
                for(int formulaId = 0 ; formulaId < origNumFormulas ; formulaId++){
                    hyperCubeMLNMStep.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMapMStep, false));
                    if((formulaId+1)%10 == 0)
                        System.out.println(formulaId+1 + " formulas done...");
                }
                // This is a reverse mapping from groundPredicates to their indices in global groundPredicates list. This
                // mapping will be filled by following function call : fgm.createGroundPredicatesFromPredsHashMap
                Map<GroundPredicate, Integer> groundPredicateIntegerMapMStep = new HashMap<>();
                List<GroundPredicate> groundPredicatesMStep = fgmMStep.createGroundPredicatesFromPredsHashMap(predsHyperCubeHashMapMStep, groundPredicateIntegerMapMStep);
                groundHiddenPredMapMToEStep.add(mapPredIdMStepToEStep(groundPredicateIntegerMapMStep, groundPredicateIntegerMapEStep, hiddenPreds));
                if(lArgs.pll)
                {
                    // This groundMLN will contain only list of groundPredicates, not groundformulas since we are not grounding MLN.
                    GroundMLN groundMLNEStep = new GroundMLN();

                    for (int gpId = 0; gpId< groundPredicatesEStep.size(); gpId++) {
                        groundMLNEStep.indexToGroundPredMap.put(gpId, groundPredicatesEStep.get(gpId));
                    }
                    groundMLNEStep.groundPredToIntegerMap = groundPredicateIntegerMapEStep;
                    groundMlnsEStep.add(groundMLNEStep);
                    Evidence truthEStep = new Evidence();
                    truthsEStep.add(truthEStep);

                    // This groundMLN will contain only list of groundPredicates, not groundformulas since we are not grounding MLN.
                    GroundMLN groundMLNMStep = new GroundMLN();

                    for (int gpId = 0; gpId< groundPredicatesMStep.size(); gpId++) {
                        groundMLNMStep.indexToGroundPredMap.put(gpId, groundPredicatesMStep.get(gpId));
                    }
                    groundMLNMStep.groundPredToIntegerMap = groundPredicateIntegerMapMStep;
                    groundMlnsMStep.add(groundMLNMStep);
                    Evidence truthMStep = parser.parseEvidence(groundMLNMStep,lArgs.truthFiles.get(i), queryPredsMStep);
                    truthsMStep.add(truthMStep);
//                if(evidPreds.size()!=0)
//                {
//                    Evidence evidence = parser.parseEvidence(groundMLN, lArgs.evidenceFiles.get(i), queryClusterPreds);
//                    reduceDomainByEvidence(mln, groundMLN, evidence);
//                    fgm.removeEvidenceGroundPreds(groundMLN, evidence);
//                }
                }

                else
                {
                    long time = System.currentTimeMillis();

                    System.out.println("Creating MRF...");
                    GroundMLN groundMlnEStep = fgmEStep.groundWithHyperCubes(hyperCubeMLNEStep, groundPredicatesEStep, groundPredicateIntegerMapEStep);
                    if(lArgs.agg)
                    {
                        groundMlnEStep.setNumConnections();
                        int count = 0;
                        for(GroundFormula gf : groundMlnEStep.groundFormulas)
                        {
                            for(int pfId : gf.parentFormulaId)
                            {
                                if(firstOrderNumConnections.get(pfId) == 0.0)
                                {
                                    firstOrderNumConnections.set(pfId, Collections.max(gf.numConnections.get(pfId)));
                                    count++;
                                }
                            }
                            if(count >= mln.formulas.size())
                                break;
                        }
                        groundMlnEStep.setEffWts(mln);
                    }
                    else
                    {
                        groundMlnEStep.setGroundFormulaWtsToSumOfParentWts(mln);
                    }
//
                    Evidence truthEStep = new Evidence();
                    truthsEStep.add(truthEStep);
                    if(lArgs.priorSoftEvidence) {
                        groundMlnEStep = fgmEStep.addSoftEvidence(groundMlnEStep, lArgs.softEvidenceFiles.get(i), lArgs.seLambda, lArgs.sePred);
                    }
                    groundMlnsEStep.add(groundMlnEStep);

                    GroundMLN groundMlnMStep = fgmMStep.groundWithHyperCubes(hyperCubeMLNMStep, groundPredicatesMStep, groundPredicateIntegerMapMStep);
                    if(lArgs.agg)
                    {
                        groundMlnMStep.setNumConnections();
                        for(Formula formula : mln.formulas)
                            firstOrderNumConnections.add(0.0);
                        int count = 0;
                        for(GroundFormula gf : groundMlnMStep.groundFormulas)
                        {
                            for(int pfId : gf.parentFormulaId)
                            {
                                if(firstOrderNumConnections.get(pfId) == 0.0)
                                {
                                    firstOrderNumConnections.set(pfId, Collections.max(gf.numConnections.get(pfId)));
                                    count++;
                                }
                            }
                            if(count >= mln.formulas.size())
                                break;
                        }
                        groundMlnMStep.setEffWts(mln);
                    }
                    else
                    {
                        groundMlnMStep.setGroundFormulaWtsToSumOfParentWts(mln);
                    }
//
                    Evidence truthMStep = parser.parseEvidence(groundMlnMStep,lArgs.truthFiles.get(i), queryPredsMStep);
                    truthsMStep.add(truthMStep);
                    if(lArgs.priorSoftEvidence) {
                        groundMlnMStep = fgmMStep.addSoftEvidence(groundMlnMStep, lArgs.softEvidenceFiles.get(i), lArgs.seLambda, lArgs.sePred);
                    }
                    groundMlnsMStep.add(groundMlnMStep);
                    System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
                    System.out.println("Total number of ground formulas in E step : " + groundMlnEStep.groundFormulas.size());
                    System.out.println("Total number of ground formulas in M step : " + groundMlnMStep.groundFormulas.size());
                    System.out.println("Total number of ground preds in E Step : " + groundMlnEStep.indexToGroundPredMap.size());
                    System.out.println("Total number of ground preds in M Step : " + groundMlnMStep.indexToGroundPredMap.size());

                }
                System.out.println("ssffdf");
            }
            else
            {
                HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList, evidPreds, truthList, queryPreds, mln, lArgs.queryEvidence);

                System.out.println("Creating Hypercubes for " + origNumFormulas + " formulas...");
                boolean isNormal = false;

                // Create a copy of this mln without formulas. In that copy, we will add hypercube formulas.
                MLN hyperCubeMLN = mln.copyMLNWithoutformulas();
                for(int formulaId = 0 ; formulaId < origNumFormulas ; formulaId++){
                    hyperCubeMLN.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMap, isNormal));
                    if((formulaId+1)%10 == 0)
                        System.out.println(formulaId+1 + " formulas done...");
                }

                // This is a reverse mapping from groundPredicates to their indices in global groundPredicates list. This
                // mapping will be filled by following function call : fgm.createGroundPredicatesFromPredsHashMap
                Map<GroundPredicate, Integer> groundPredicateIntegerMap = new HashMap<>();
                List<GroundPredicate> groundPredicates = fgm.createGroundPredicatesFromPredsHashMap(predsHyperCubeHashMap, groundPredicateIntegerMap);

                if(lArgs.pll)
                {
                    // This groundMLN will contain only list of groundPredicates, not groundformulas since we are not grounding MLN.
                    GroundMLN groundMLN = new GroundMLN();

                    for (int gpId = 0; gpId< groundPredicates.size(); gpId++) {
                        groundMLN.indexToGroundPredMap.put(gpId, groundPredicates.get(gpId));
                    }
                    groundMLN.groundPredToIntegerMap = groundPredicateIntegerMap;
                    groundMlns.add(groundMLN);
                    Evidence truth = parser.parseEvidence(groundMLN,lArgs.truthFiles.get(i), queryPreds);
                    truths.add(truth);
                    if(lArgs.agg)
                    {
                        for (int j = 0; j < mln.formulas.size(); j++) {
                            Formula f = mln.formulas.get(j);
                            Set<Term> fTerms = new HashSet<Term>(f.terms);
                            double maxConnections = 0.0;
                            for(WClause clause : f.clauses)
                            {
                                for(Atom atom : clause.atoms)
                                {
                                    Set<Term> atomTerms = new HashSet<Term>(atom.terms);
                                    Set<Term> tempFTerms = new HashSet<Term>(fTerms);
                                    tempFTerms.removeAll(atomTerms);
                                    int num = 1;
                                    for(Term t : tempFTerms)
                                    {
                                        num *= t.domain.size();
                                    }
                                    if(num > maxConnections)
                                    {
                                        maxConnections = num;
                                    }
                                }
                            }
                            firstOrderNumConnections.set(j, firstOrderNumConnections.get(j)+maxConnections);
                        }
                    }

//                if(evidPreds.size()!=0)
//                {
//                    Evidence evidence = parser.parseEvidence(groundMLN, lArgs.evidenceFiles.get(i), queryClusterPreds);
//                    reduceDomainByEvidence(mln, groundMLN, evidence);
//                    fgm.removeEvidenceGroundPreds(groundMLN, evidence);
//                }
                }
                else
                {
                    long time = System.currentTimeMillis();

                    System.out.println("Creating MRF...");
                    GroundMLN groundMln = fgm.groundWithHyperCubes(hyperCubeMLN, groundPredicates, groundPredicateIntegerMap);
                    if(lArgs.agg)
                    {
//                        groundMln.setNumConnections();
//                        int count = 0;
//                        for(GroundFormula gf : groundMln.groundFormulas)
//                        {
//                            for(int pfId : gf.parentFormulaId)
//                            {
//                                if(firstOrderNumConnections.get(pfId) == 0.0)
//                                {
//                                    firstOrderNumConnections.set(pfId, Collections.max(gf.numConnections.get(pfId)));
//                                    count++;
//                                }
//                            }
//                            if(count >= mln.formulas.size())
//                                break;
//                        }
                        if(lArgs.agg)
                        {
                            for (int j = 0; j < mln.formulas.size(); j++) {
                                Formula f = mln.formulas.get(j);
                                Set<Term> fTerms = new HashSet<Term>(f.terms);
                                double maxConnections = 0.0;
                                for(WClause clause : f.clauses)
                                {
                                    for(Atom atom : clause.atoms)
                                    {
                                        Set<Term> atomTerms = new HashSet<Term>(atom.terms);
                                        Set<Term> tempFTerms = new HashSet<Term>(fTerms);
                                        tempFTerms.removeAll(atomTerms);
                                        int num = 1;
                                        for(Term t : tempFTerms)
                                        {
                                            num *= t.domain.size();
                                        }
                                        if(num > maxConnections)
                                        {
                                            maxConnections = num;
                                        }
                                    }
                                }
                                firstOrderNumConnections.set(j, firstOrderNumConnections.get(j)+maxConnections);
                            }
                        }

                        groundMln.setEffWts(mln);
                    }
                    else
                    {
                        groundMln.setGroundFormulaWtsToSumOfParentWts(mln);
                    }

                    groundMlns.add(groundMln);
                    Evidence truth = parser.parseEvidence(groundMln,lArgs.truthFiles.get(i), queryPreds);
                    truths.add(truth);
                    if(lArgs.priorSoftEvidence) {
                        groundMln = fgm.addSoftEvidence(groundMln, lArgs.softEvidenceFiles.get(i), lArgs.seLambda, lArgs.sePred);
                    }
                    System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
                    System.out.println("Total number of ground formulas : " + groundMln.groundFormulas.size());
                    System.out.println("Total number of ground preds : " + groundMln.indexToGroundPredMap.size());

                }
            }

        }

        if(lArgs.withEM)
        {
            if(lArgs.pll)
                wl = new GenLearner(mlns, groundMlnsMStep, groundMlnsEStep, truthsMStep, truthsEStep, lArgs, groundHiddenPredMapMToEStep);
            else
                wl = new DiscLearn(mlns, groundMlnsMStep, groundMlnsEStep, truthsMStep, truthsEStep, lArgs, groundHiddenPredMapMToEStep);
        }
        else
        {
            if(lArgs.pll)
                wl = new GenLearner(mlns, groundMlns, null, truths, null, lArgs, null);
            else
                wl = new DiscLearn(mlns, groundMlns, null, truths, null, lArgs, null);
        }


        wl.learnWeights();
        wl.writeWeights(lArgs.mlnFile, lArgs.outFile);
        System.out.println("Learning done !!!");
        System.out.println("Final weights are : " + Arrays.toString(wl.weights));
        System.out.println("Total Time taken : " + Timer.time((System.currentTimeMillis() - totaltime)/1000.0));
    }

    private static void reduceDomainByEvidence(MLN mln, GroundMLN groundMln, Evidence evidence) {
        // Key : clusterNumber of constant, value : set of constants in that key cluster
        Map<Integer, Set<Integer>> valToConstantsMap = new HashMap<>();

        for(Integer gpId : evidence.predIdVal.keySet())
        {
            GroundPredicate gp = groundMln.indexToGroundPredMap.get(gpId);
            int trueVal = evidence.predIdVal.get(gpId);
            MyAssert.assume(gp.terms.size() == 1);
            int constant = gp.terms.get(0);
            if(!valToConstantsMap.containsKey(trueVal))
            {
                valToConstantsMap.put(trueVal, new HashSet<Integer>());
            }
            valToConstantsMap.get(trueVal).add(constant);
        }

        for(Formula f : mln.formulas)
        {
            Map<Term, Integer> stTermsToClusterMap = new HashMap<>();
            for (int clauseId = f.clauses.size()-1; clauseId >=0 ; clauseId--) {
                WClause clause = f.clauses.get(clauseId);
                if(clause.atoms.size() == 1 && clause.atoms.get(0).symbol.symbol.equals("st"))
                {
                    stTermsToClusterMap.put(clause.atoms.get(0).terms.get(0), clause.valTrue.get(0));
                    f.clauses.remove(clauseId);
                }
            }

            for(WClause clause : f.clauses)
            {
                for(Atom atom : clause.atoms)
                {
                    for(Term term : atom.terms)
                    {
                        if(stTermsToClusterMap.containsKey(term))
                        {
                            int clusterNum = stTermsToClusterMap.get(term);
                            term.domain = new ArrayList<>(valToConstantsMap.get(clusterNum));
                        }
                    }
                }
            }
        }

    }

    private static Set<String> getHiddenPreds(Parser parser, MLN mln, LearnArgs lArgs, Set<String> closedWorldPreds) {
        Set<String> allPreds = new HashSet<String>();
        for(PredicateSymbol ps : mln.symbols)
        {
            allPreds.add(ps.symbol);
        }
        allPreds.removeAll(parser.evidPreds);
        allPreds.removeAll(parser.queryPreds);
        if(allPreds.size() > 0 && !lArgs.withEM)
        {
            throw new ParameterException("hidden Preds found in MLN : " + allPreds + ", provide withEM flag!!!");
        }
        else if(lArgs.withEM && allPreds.isEmpty())
        {
            throw new ParameterException("Must provide hidden preds withEM!!!");
        }
        closedWorldPreds.addAll(parser.evidPreds);
        //closedWorldPreds.addAll(parser.queryPreds);
        return allPreds;
    }

    private static Map<Integer,Integer> mapPredIdMStepToEStep(Map<GroundPredicate, Integer> gpToIntegerMapMStep, Map<GroundPredicate, Integer> gpToIntegerMapEStep, Set<String> hiddenPreds) {
        Map<Integer, Integer> subtypeMap = new HashMap<>();
        for(GroundPredicate gp : gpToIntegerMapMStep.keySet()){
            if(hiddenPreds.contains(gp.symbol.symbol)){
                if(gpToIntegerMapEStep.containsKey(gp))
                    subtypeMap.put(gpToIntegerMapMStep.get(gp), gpToIntegerMapEStep.get(gp));
                else
                    subtypeMap.put(gpToIntegerMapMStep.get(gp), -1);
            }
        }
        return subtypeMap;
    }

    private static LearnArgs parseArgs(String[] args){
        LearnArgs la = new LearnArgs();
        JCommander jc = JCommander.newBuilder().addObject(la).build();
        jc.setColumnSize(200);
        jc.setProgramName("LearnTest");
        try
        {
            jc.parse(args);
            la.validate();
            System.out.println(la);
        }
        catch(ParameterException p)
        {
            System.out.println(p.getMessage());
            jc.usage();
            LearnArgs.GibbsParamConverter.jc.setProgramName("-infer");
            LearnArgs.GibbsParamConverter.jc.usage();
            LearnArgs.CGParamsConverter.jc.setProgramName("-cg");
            LearnArgs.CGParamsConverter.jc.usage();
            System.exit(-1);
        }
        return la;
    }

    public enum Method {
        CG, LBFGS, VP

    }
}
