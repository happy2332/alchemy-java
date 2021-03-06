package org.utd.cs.mln.inference;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
import org.utd.cs.mln.alchemy.util.MyAssert;
import org.utd.cs.mln.alchemy.util.Parser;
import org.utd.cs.mln.learning.LearnTest;
import org.utd.cs.mln.lmap.MlnToHyperCube;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by Happy on 2/28/17.
 */
public class InferTest {

    public static void main(String[] args) throws FileNotFoundException, CloneNotSupportedException, PredicateNotFound {
        System.out.println("Parsing command line arguments...");
        InferArgs iArgs = parseArgs(args);
        FullyGrindingMill fgm = new FullyGrindingMill();
        MLN mln = new MLN(iArgs.priorSoftEvidence);
        Parser parser = new Parser(mln);
        List<String> files = new ArrayList<>();
        files.add(iArgs.goldFile);
        files.add(iArgs.evidFile);
        parser.setTruthEvidFiles(files);
        parser.parseInputMLNFile(iArgs.mlnFile);
        LearnTest.firstOrderNumConnections = new ArrayList<>(Collections.nCopies(mln.formulas.size(),0.0));

        Set<String> evidPreds = parser.evidPreds, queryPreds = parser.queryPreds;
        validateEvidQueryPreds(mln, evidPreds, queryPreds);

        boolean isgroundwithhypercube = true;
        MLN hyperCubeMLN = null;
        HashMap<PredicateSymbol, ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = null;
        if(isgroundwithhypercube) {
            List<FirstEvidence> evidList = parser.parseInputEvidenceFile(iArgs.evidFile);
            List<FirstEvidence> truthList = parser.parseInputEvidenceFile(iArgs.goldFile);
            MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
            predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList, evidPreds, truthList, queryPreds, mln, iArgs.queryEvidence);
            System.out.println("Creating Hypercubes...");
            int origNumFormulas = mln.formulas.size();
            boolean isNormal = false;

            // Create a copy of this mln without formulas. In that copy, we will add hypercube formulas.
            hyperCubeMLN = mln.copyMLNWithoutformulas();
            for (int formulaId = 0; formulaId < origNumFormulas; formulaId++) {
                hyperCubeMLN.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMap, isNormal));
            }
        }
        GroundMLN groundMln;
        Evidence gold;
        long time = System.currentTimeMillis();
        Map<GroundPredicate, Integer> groundPredicateIntegerMap = new HashMap<>();

        if(isgroundwithhypercube){
            List<GroundPredicate> groundPredicates = fgm.createGroundPredicatesFromPredsHashMap(predsHyperCubeHashMap, groundPredicateIntegerMap);
            System.out.println("Creating MRF...");
            groundMln = fgm.groundWithHyperCubes(hyperCubeMLN, groundPredicates, groundPredicateIntegerMap);
            if(iArgs.agg)
            {
                groundMln.setNumConnections();
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
                    LearnTest.firstOrderNumConnections.set(j, LearnTest.firstOrderNumConnections.get(j)+maxConnections);
                }
                groundMln.setEffWts(mln);
            }
            else
            {
                groundMln.setGroundFormulaWtsToSumOfParentWts(mln);
            }

//            Evidence evidence = parser.parseEvidence(groundMln, iArgs.evidFile, queryClusterPreds);
//            fgm.removeEvidenceGroundPreds(groundMln, evidence);
            gold = parser.parseEvidence(groundMln, iArgs.goldFile, queryPreds);
        }

        else
        {
            List<GroundPredicate> groundPredicates = fgm.createGroundPredicates(mln, groundPredicateIntegerMap, queryPreds);
            groundMln = fgm.ground(mln, groundPredicates, groundPredicateIntegerMap);
            Evidence evidence = parser.parseEvidence(groundMln, iArgs.evidFile, evidPreds);
            Set<String> hiddenPreds = new HashSet<>();
//        Map<Integer, List<Integer>> featureVectors = fgm.getFeatureVectors(groundMln, mln.formulas.size(), evidence, "person", varTypeToDomain.get("person"), false);
//        writeFeatures(featureVectors,5,100);

            gold = parser.parseEvidence(groundMln, iArgs.goldFile, queryPreds);
//            System.out.println("Total number of ground formulas (before handling evidence): " + groundMln.groundFormulas.size());
            groundMln = fgm.handleEvidence(groundMln, evidence, gold, evidPreds, queryPreds, hiddenPreds, false, iArgs.queryEvidence);

//            System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));
//
//            System.out.println("Total number of ground formulas (after handling evidence): " + groundMln.groundFormulas.size());
            //fgm.removeEvidenceGroundPreds(groundMln, evidence);

        }


        if(iArgs.priorSoftEvidence) {
            groundMln = fgm.addSoftEvidence(groundMln, iArgs.softEvidenceFile, iArgs.seLambda, iArgs.sePred);
        }

//        GibbsSampler_v2 gs = new GibbsSampler_v2(mln, groundMln, gold, 100, 1000, false, true, false, false);

        System.out.println("Total number of ground formula : " + groundMln.groundFormulas.size());
        System.out.println("Total permutation count : " + FullyGrindingMill.permutationCount);
        System.out.println("Total number of ground predicates : " + groundMln.indexToGroundPredMap.size());
        System.out.println("Time taken to ground MLN : "+Timer.time((System.currentTimeMillis() - time) / 1000.0));
        State state = new State(mln, groundMln);
        state.setTruthVals(gold);
        GibbsParams gibbsparams = new GibbsParams(iArgs.gibbsParam);
        Inference inference = new GibbsSampler_v3(state,-1,false,false, false,gibbsparams, true);
        PrintWriter writer = null;
        writer = new PrintWriter(new FileOutputStream(iArgs.outFile));
//        gs.infer(true, true);
//        gs.writeMarginal(writer);
        //inference.writeNetwork(writer);
        //PseudoLogLikelihood pll = new PseudoLogLikelihood(state);
        //System.out.println("pll : "+pll.getPseudoLogLikelihood());
        inference.infer();
        inference.writeProbs(writer);
        writer.close();
    }

    private static void validateEvidQueryPreds(MLN mln, Set<String> evidPreds, Set<String> queryPreds) {
        Set<String> allPreds = new HashSet<>();
        for(PredicateSymbol ps : mln.symbols)
        {
            allPreds.add(ps.symbol);
        }
        allPreds.removeAll(evidPreds);
        allPreds.removeAll(queryPreds);
        MyAssert.assume(allPreds.size()==0);
    }

    private static void writeFeatures(Map<Integer, List<Integer>> featureVectors, int db, int evidPer) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter("data/temp/feature/features_infer."+db+"."+evidPer+".txt");
        for(int constant : featureVectors.keySet())
        {
            for(int val : featureVectors.get(constant))
            {
                pw.write(val+",");
            }
            pw.write(constant+"\n");
        }
        pw.close();
    }

    private static InferArgs parseArgs(String[] args){
        InferArgs ia = new InferArgs();
        JCommander jc = JCommander.newBuilder().addObject(ia).build();
        jc.setColumnSize(200);
        jc.setProgramName("InferTest");
        try
        {
            jc.parse(args);
            ia.validate();
            System.out.println(ia);
        }
        catch(ParameterException p)
        {
            System.out.println(p.getMessage());
            jc.usage();
            System.exit(-1);
        }
        return ia;
    }

}
