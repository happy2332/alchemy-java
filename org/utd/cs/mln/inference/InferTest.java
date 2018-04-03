package org.utd.cs.mln.inference;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import org.utd.cs.gm.utility.Timer;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.alchemy.util.FullyGrindingMill;
import org.utd.cs.mln.alchemy.util.Parser;
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
        MLN mln = new MLN();
        Parser parser = new Parser(mln);
        parser.parseInputMLNFile(iArgs.mlnFile);
        Set<String> allPreds = new HashSet<String>();
        List<String> closedWorldPreds = new ArrayList<>();
        for(PredicateSymbol ps : mln.symbols)
        {
            allPreds.add(ps.symbol);
        }
        if(iArgs.evidPreds == null)
        {
            iArgs.evidPreds = new ArrayList<>();
            iArgs.evidPreds.addAll(allPreds);
            iArgs.evidPreds.removeAll(iArgs.queryPreds);
        }
        closedWorldPreds.addAll(iArgs.evidPreds);
        closedWorldPreds.removeAll(iArgs.queryPreds);
        String files[] = new String[2];
        files[0] = iArgs.evidFile;
        files[1] = iArgs.goldFile;
        Map<String, Set<Integer>> varTypeToDomain = parser.collectDomain(files);
        mln.overWriteDomain(varTypeToDomain);
        boolean isgroundwithhypercube = true;
        GroundMLN groundMln = null;
        Evidence gold = null;
        if(isgroundwithhypercube)
        {
            List<FirstEvidence> evidList = parser.parseInputEvidenceFile(iArgs.evidFile);

            MlnToHyperCube mlnToHyperCube = new MlnToHyperCube();
            HashMap<PredicateSymbol,ArrayList<ArrayList<HyperCube>>> predsHyperCubeHashMap = mlnToHyperCube.createPredsHyperCube(evidList,mln,closedWorldPreds);
            System.out.println("Creating Hypercubes...");
            long time = System.currentTimeMillis();
            int origNumFormulas = mln.formulas.size();
            boolean isNormal = false;

            for(int formulaId = 0 ; formulaId < origNumFormulas ; formulaId++){
                mln.formulas.addAll(mlnToHyperCube.createFormulaHyperCube(mln.formulas.get(formulaId), predsHyperCubeHashMap, isNormal));
            }
            for(int formulaId = origNumFormulas-1 ; formulaId >= 0 ; formulaId--){
                mln.formulas.remove(formulaId);
            }
            System.out.println("Creating MRF...");
            groundMln = fgm.groundWithHyperCubes(mln);
            gold = parser.parseEvidence(groundMln, iArgs.goldFile);
            System.out.println("Total number of ground formula : " + groundMln.groundFormulas.size());
            System.out.println("Total permutation count : " + FullyGrindingMill.permutationCount);
            System.out.println("Time taken to ground MLN : "+Timer.time((System.currentTimeMillis() - time) / 1000.0));
        }

        else
        {
            groundMln = fgm.ground(mln);
            long time = System.currentTimeMillis();
            Evidence evidence = parser.parseEvidence(groundMln, iArgs.evidFile);
            List<String> hiddenPreds = new ArrayList<>();
//        Map<Integer, List<Integer>> featureVectors = fgm.getFeatureVectors(groundMln, mln.formulas.size(), evidence, "person", varTypeToDomain.get("person"), false);
//        writeFeatures(featureVectors,5,100);

            Map<GroundPredicate, Integer> groundPredicateIntegerMap = null;
            gold = parser.parseEvidence(groundMln, iArgs.goldFile);
            System.out.println("Total number of ground formulas (before handling evidence): " + groundMln.groundFormulas.size());
            groundMln = fgm.handleEvidence(groundMln, evidence, gold, iArgs.evidPreds, iArgs.queryPreds, hiddenPreds, false, groundPredicateIntegerMap, iArgs.queryEvidence);

            System.out.println("Time taken to create MRF : " + Timer.time((System.currentTimeMillis() - time) / 1000.0));

            System.out.println("Total number of ground formulas (after handling evidence): " + groundMln.groundFormulas.size());

        }


        if(iArgs.priorSoftEvidence) {
            groundMln = fgm.addSoftEvidence(groundMln, iArgs.softEvidenceFile, iArgs.seLambda, iArgs.sePred);
        }

        //GibbsSampler_v2 gs = new GibbsSampler_v2(mln, groundMln, gold, NumBurnIn, NumSamples, trackFormulaCounts, calculateMarginal, priorSoftEvidence, false);
        State state = new State(mln, groundMln);
        GibbsParams gibbsparams = new GibbsParams();
        Inference inference = new GibbsSampler_v3(state,-1,false,false,gibbsparams);
        PrintWriter writer = null;
        try{
            writer = new PrintWriter(new FileOutputStream(iArgs.outFile));
        }
        catch (IOException e) {
        }
//        gs.infer(true, true);
//        gs.writeMarginal(writer);
        //inference.writeNetwork(writer);
        inference.init();
        inference.infer();
        inference.writeProbs(writer);
        writer.close();
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
