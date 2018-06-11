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
        //String filename = "/Users/Happy/phd/experiments/without/data/Imdb/mln/imdb_mln.txt";
        //String evidence_files[] = "/Users/Happy/phd/experiments/without/data/Imdb/db/empty_file.txt".split(",");
        //String train_files[] = "/Users/Happy/phd/experiments/without/data/Imdb/db/imdb.1_train.txt".split(",");
        //String out_file = "/Users/Happy/phd/experiments/without/data/Imdb/results/imdb_results_out.txt";
        //List<String> evidence_preds = Arrays.asList("".split(","));
        //List<String> query_preds = Arrays.asList("actor,director,movie,workedUnder".split(","));;

        // Create required data structures

        // List of MLNs of size numDb. All MLNs are same except for domains for predicates.
        List<MLN> mlns = new ArrayList<>();

        // List of truth evidence of size numDb.
        List<Evidence> truths = new ArrayList<>();

        // List of groundMlns of size numDb.
        List<GroundMLN> groundMlns = new ArrayList<>();

        List<GibbsSampler_v2> inferences = new ArrayList<>();
        List<GibbsSampler_v2> inferencesEM = new ArrayList<GibbsSampler_v2>();
        Map<Integer, Integer>[] subtypeMaps = new Map[numDb];
        Set<String> evidPreds = null, queryPreds = null, hiddenPreds = null;
        Set<String> closedWorldPreds = new HashSet<>();

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
            //TODO : for now, we send only truth file to collectDomain, otherwise send both truth and evidence.
            // We assume that there is no new constant in evidence.
            //Map<String, Set<Integer>> varTypeToDomain = parser.collectDomain(files);
            parser.setTruthEvidFiles(files);
            parser.parseInputMLNFile(lArgs.mlnFile);
            // Need to set query and evidence predicates only once, they will be same for all DBs
            if(i == 0)
            {
                evidPreds = parser.evidPreds;
                queryPreds = parser.queryPreds;
                closedWorldPreds.addAll(evidPreds);
                closedWorldPreds.removeAll(queryPreds);
                hiddenPreds = getHiddenPreds(parser, mln, lArgs, closedWorldPreds);
            }

            Map<GroundPredicate, Integer> groundPredicateIntegerMap = new HashMap<>();
            Set<String> queryClusterPreds = new HashSet<String>(queryPreds);
            queryClusterPreds.add("st");
            List<GroundPredicate> groundPredicates = fgm.createGroundPredicates(mln, groundPredicateIntegerMap, queryClusterPreds);
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
                if(evidPreds.size()!=0)
                {
                    Evidence evidence = parser.parseEvidence(groundMLN, lArgs.evidenceFiles.get(i), queryClusterPreds);
                    reduceDomainByEvidence(mln, groundMLN, evidence);
                    fgm.removeEvidenceGroundPreds(groundMLN, evidence);
                }
            }
            else
            {
                long time = System.currentTimeMillis();
                List<FirstEvidence> evidList = parser.parseInputEvidenceFile(lArgs.evidenceFiles.get(i));
                MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
                HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList,mln,closedWorldPreds);
                System.out.println("Creating Hypercubes...");
                int origNumFormulas = mln.formulas.size();
                boolean isNormal = false;

                // Create a copy of this mln without formulas. In that copy, we will add hypercube formulas.
                MLN hyperCubeMLN = mln.copyMLNWithoutformulas();
                for(int formulaId = 0 ; formulaId < origNumFormulas ; formulaId++){
                    hyperCubeMLN.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMap, isNormal));
                }
                System.out.println("Creating MRF...");
                GroundMLN groundMln = fgm.groundWithHyperCubes(hyperCubeMLN, groundPredicates, groundPredicateIntegerMap);
                if(lArgs.agg)
                {
                    groundMln.setNumConnections();
                    for(Formula formula : mln.formulas)
                        firstOrderNumConnections.add(0.0);
                    int count = 0;
                    for(GroundFormula gf : groundMln.groundFormulas)
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
                    groundMln.setEffWts(mln);
                }
                else
                {
                    groundMln.setGroundFormulaWtsToSumOfParentWts(mln);
                }

                groundMlns.add(groundMln);
                if(evidPreds.size()!=0)
                {
                    Evidence evidence = parser.parseEvidence(groundMln, lArgs.evidenceFiles.get(i), queryClusterPreds);
                    reduceDomainByEvidence(mln, groundMln, evidence);
                    fgm.removeEvidenceGroundPreds(groundMln, evidence);
                }
                Evidence truth = parser.parseEvidence(groundMln,lArgs.truthFiles.get(i), queryPreds);
                truths.add(truth);
                if(lArgs.priorSoftEvidence) {
                    groundMln = fgm.addSoftEvidence(groundMln, lArgs.softEvidenceFiles.get(i), lArgs.seLambda, lArgs.sePred);
                }
                boolean withEM = false;
//                if(lArgs.withEM){
//                    withEM = true;
//                    Set<String> evidEmPreds = new HashSet<String>();
//                    evidEmPreds.addAll(evidPreds);
//                    evidEmPreds.addAll(queryPreds);
////                Map<GroundPredicate,Integer> gpToIntegerMap = new HashMap<>();
//                    GroundMLN EMNewGroundMln = fgm.handleEvidence(groundMln, Evidence.mergeEvidence(evidence,truth), null, evidEmPreds, null, hiddenPreds, withEM, lArgs.queryEvidence);
////                newGroundMln = fgm.handleEvidence(groundMln, evidence, truth, evidPreds, queryPreds, hiddenPreds, false, gpToIntegerMap);
//                    subtypeMaps[i] = mapPredIdEmToNormal(newGroundMln.groundPredToIntegerMap, EMNewGroundMln.groundPredToIntegerMap, hiddenPreds);
//
//                    //GibbsSampler_v2 gsEM = new GibbsSampler_v2(mln, EMNewGroundMln, truth, 20, lArgs.numEMSamples, false, false, lArgs.priorSoftEvidence, true);
//                    //inferencesEM.add(gsEM);
//                }

                //GibbsSampler_v2 gs = new GibbsSampler_v2(mln, newGroundMln, truth, 100, lArgs.gibbsParam.samplesPerTest, true, false, lArgs.priorSoftEvidence, false);
                //inferences.add(gs);

                System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
                System.out.println("Total number of ground formulas : " + groundMln.groundFormulas.size());
                System.out.println("Total number of ground preds : " + groundMln.indexToGroundPredMap.size());

            }
        }

        WeightLearner wl;
        if(lArgs.pll)
            wl = new GenLearner(mlns, groundMlns, null, truths, null, lArgs);
        else
            wl = new DiscLearn(mlns, groundMlns, null, truths, null, lArgs);

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

    private static Map<Integer,Integer> mapPredIdEmToNormal(Map<GroundPredicate, Integer> gpToIntegerMap, Map<GroundPredicate, Integer> gpToIntegerMapEM, Set<String> hiddenPreds) {
        Map<Integer, Integer> subtypeMap = new HashMap<>();
        for(GroundPredicate gp : gpToIntegerMap.keySet()){
            if(hiddenPreds.contains(gp.symbol.symbol)){
                if(gpToIntegerMapEM.containsKey(gp))
                    subtypeMap.put(gpToIntegerMap.get(gp), gpToIntegerMapEM.get(gp));
                else
                    subtypeMap.put(gpToIntegerMap.get(gp), -1);
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
