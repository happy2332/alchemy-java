package org.utd.cs.mln.learning;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
import org.utd.cs.mln.alchemy.util.Parser;
import org.utd.cs.mln.inference.GibbsSampler_v2;

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
        FullyGrindingMill fgm = new FullyGrindingMill();
        Map<Integer, Integer>[] subtypeMaps = new Map[numDb];
        List<String> closedWorldPreds = new ArrayList<>();

        for(int i = 0 ; i < numDb ; i++) {
            MLN mln = new MLN();
            mlns.add(mln);
            Parser parser = new Parser(mln);
            parser.parseInputMLNFile(lArgs.mlnFile);
            // Need to set query and evidence predicates only once, they will be same for all DBs
            if(i == 0)
                setQueryEvidPreds(mln, lArgs, closedWorldPreds);
            System.out.println("Reading DB file "+(i+1));
            String files[] = new String[1];
            //files[0] = evidenceFiles[i];
            files[0] = lArgs.truthFiles.get(i);
            //TODO : for now, we send only truth file to collectDomain, otherwise send both truth and evidence.
            // We assume that there is no new constant in evidence.
            Map<String, Set<Integer>> varTypeToDomain = parser.collectDomain(files);
            mln.overWriteDomain(varTypeToDomain);

            if(lArgs.genLearn)
            {
                // This groundMLN will contain only list of groundPredicates, not groundformulas since we are not grounding MLN.
                GroundMLN groundMLN = new GroundMLN();
                List<GroundPredicate> groundPredicates = fgm.createGroundPredicates(mln, varTypeToDomain);
                groundMLN.groundPredicates = groundPredicates;
                groundMLN.setGroundPredToIntegerMap();
                groundMlns.add(groundMLN);
                Evidence truth = parser.parseEvidence(groundMLN,lArgs.truthFiles.get(i));
                truths.add(truth);
                // Create genlearner object

                //gl.addState(mln, groundMLN, truth);
                continue;
            }
            System.out.println("hello");
            System.out.println("Creating MRF...");
            long time = System.currentTimeMillis();
            GroundMLN groundMln = fgm.ground(mln);
            System.out.println("Total number of ground formulas before handling evidence : " + groundMln.groundFormulas.size());
            Evidence evidence = parser.parseEvidence(groundMln,lArgs.evidenceFiles.get(i));
            Evidence truth = parser.parseEvidence(groundMln,lArgs.truthFiles.get(i));
            GroundMLN newGroundMln = null;

            Map<GroundPredicate,Integer> gpToIntegerMap = new HashMap<>();
            if(lArgs.priorSoftEvidence) {
                groundMln = fgm.addSoftEvidence(groundMln, lArgs.softEvidenceFiles.get(i), lArgs.seLambda, lArgs.sePred);
            }
            boolean withEM = false;
            newGroundMln = fgm.handleEvidence(groundMln, evidence, truth, lArgs.evidPreds, lArgs.queryPreds, lArgs.hiddenPreds, lArgs.withEM, gpToIntegerMap, lArgs.queryEvidence);
//            newGroundMln = groundMln;

            if(lArgs.withEM){
                withEM = true;
                List<String> evidEmPreds = new ArrayList<String>();
                evidEmPreds.addAll(lArgs.evidPreds);
                evidEmPreds.addAll(lArgs.queryPreds);
//                Map<GroundPredicate,Integer> gpToIntegerMap = new HashMap<>();
                Map<GroundPredicate,Integer> gpToIntegerMapEM = new HashMap<>();
                GroundMLN EMNewGroundMln = fgm.handleEvidence(groundMln, Evidence.mergeEvidence(evidence,truth), null, evidEmPreds, null, lArgs.hiddenPreds, withEM, gpToIntegerMapEM, lArgs.queryEvidence);
//                newGroundMln = fgm.handleEvidence(groundMln, evidence, truth, evidPreds, queryPreds, hiddenPreds, false, gpToIntegerMap);
                subtypeMaps[i] = mapPredIdEmToNormal(gpToIntegerMap, gpToIntegerMapEM, lArgs.hiddenPreds);

                GibbsSampler_v2 gsEM = new GibbsSampler_v2(mln, EMNewGroundMln, truth, 20, lArgs.numEMSamples, false, false, lArgs.priorSoftEvidence, true);
                inferencesEM.add(gsEM);
            }

            GibbsSampler_v2 gs = new GibbsSampler_v2(mln, newGroundMln, truth, 100, lArgs.gibbsParam.samplesPerTest, true, false, lArgs.priorSoftEvidence, false);
            inferences.add(gs);

            System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time)/1000.0));
            System.out.println("Total number of ground formulas : " + newGroundMln.groundFormulas.size());
            System.out.println("Total number of ground preds : " + newGroundMln.groundPredicates.size());
        }

        WeightLearner wl = null;
        if(lArgs.genLearn)
        {
            wl = new GenLearner(mlns, groundMlns, truths, lArgs);
            wl.learnWeights();
        }
        else
        {
            // Start learning
            DiscLearner dl = new DiscLearner(inferences, inferencesEM, subtypeMaps, lArgs.numIter, 100.0, lArgs.minllChange, Double.MAX_VALUE, lArgs.withEM, true, lArgs.usePrior, lArgs.priorSoftEvidence, lArgs.seLambda, Method.CG);

            dl.learnWeights();
        }
        wl.writeWeights(lArgs.mlnFile, lArgs.outFile);
        System.out.println("Learning done !!!");
        System.out.println("Final weights are : " + Arrays.toString(wl.weights));
        System.out.println("Total Time taken : " + Timer.time((System.currentTimeMillis() - totaltime)/1000.0));
    }

    private static void setQueryEvidPreds(MLN mln, LearnArgs lArgs, List<String> closedWorldPreds) {
        Set<String> allPreds = new HashSet<String>();
        for(PredicateSymbol ps : mln.symbols)
        {
            allPreds.add(ps.symbol);
        }

        // set query, evidence, hidden preds as well as closedworldpreds
        // If discriminative learning, then all other predicates except query and hidden are evidence predicates
        if(lArgs.discLearn)
        {

            Set<String> allbutEvidPreds = new HashSet<String>(lArgs.queryPreds);
            if(lArgs.withEM)
            {
                allbutEvidPreds.addAll(lArgs.hiddenPreds);
            }
            allPreds.removeAll(allbutEvidPreds);
            lArgs.evidPreds = new ArrayList<>(allPreds);
            closedWorldPreds.addAll((lArgs.evidPreds));
        }
        // If generative learning, then all except hidden preds are query preds and there is no evidence predicates
        else
        {
            if(lArgs.withEM)
            {
                allPreds.removeAll(lArgs.hiddenPreds);
            }
            lArgs.queryPreds = new ArrayList<>(allPreds);
            closedWorldPreds.addAll((lArgs.queryPreds));
        }

    }

    private static Map<Integer,Integer> mapPredIdEmToNormal(Map<GroundPredicate, Integer> gpToIntegerMap, Map<GroundPredicate, Integer> gpToIntegerMapEM, List<String> hiddenPreds) {
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
            LearnArgs.GibbsParamConverter.jc.usage();
            System.exit(-1);
        }
        return la;
    }

    public enum Method {
        CG, LBFGS, VP

    }
}
