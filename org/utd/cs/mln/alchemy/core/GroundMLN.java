package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.gm.utility.Pair;
import org.utd.cs.mln.alchemy.util.Aggregator;
import org.utd.cs.mln.alchemy.util.VecOperations;

import java.util.*;

/**
 * Created by Happy on 2/28/17.
 */
public class GroundMLN {
    public List<GroundPredicateSymbol> symbols = new ArrayList<>();
    public Map<Integer, GroundPredicate> indexToGroundPredMap = new HashMap<>();
    public Map<GroundPredicate,Integer> groundPredToIntegerMap = new HashMap<>();
    public List<GroundFormula> groundFormulas = new ArrayList<>();
    public Map<Integer, List<Double>> softEvidencePerPredPerVal = new HashMap<>();
    public boolean priorSoftEvidence = false;
    // TODO : Fill softEvidencPerPredPerVal in addSoftEvidence function

//    /**
//     * This function fills in {@link GroundMLN#groundPredToIntegerMap} which creates a mapping from groundpredicate to its index in {@link GroundMLN#groundPredicates} list.
//     * So if groundPredicates[i] = gp, then groundPredToIntegerMap[gp] = i.
//     * This is useful to search a groundpredicate gp in groundPredicates list.
//     */
//    public void setGroundPredToIntegerMap()
//    {
//        for (int i = 0; i < groundPredicates.size(); i++) {
//            groundPredToIntegerMap.put(groundPredicates.get(i),i);
//        }
//    }

    public int getNumGroundFormulas() {
        int count = 0;
        for(GroundFormula gf : groundFormulas)
        {
            for(int nc : gf.numCopies)
            {
                count += nc;
            }
        }
        return count;
    }

//    public void setGroundFormulaWtsToParentWts(MLN mln) {
//        for(GroundFormula gf : groundFormulas)
//        {
//            List<Integer> parentFormulaId = gf.parentFormulaId;
//            List<Integer> numCopies = gf.numCopies;
//            double weight = 0.0;
//            for (int i = 0; i < parentFormulaId.size(); i++) {
//                weight += numCopies.get(i) * mln.formulas.get(parentFormulaId.get(i)).weight.getValue();
//            }
//            gf.weight = new LogDouble(weight, true);
//        }
//    }

    public void setGroundFormulaWtsToSumOfParentWts(MLN mln) {
        for(GroundFormula gf : groundFormulas)
        {
            List<Integer> parentFormulaId = gf.parentFormulaId;
            List<Integer> numCopies = gf.numCopies;
            if(!parentFormulaId.isEmpty())
            {
                double weight = 0.0;
                for (int i = 0; i < parentFormulaId.size(); i++) {
                    weight += numCopies.get(i) * mln.formulas.get(parentFormulaId.get(i)).weight.getValue();
                }
                gf.weight = new LogDouble(weight, true);
            }
            else
                gf.weight = new LogDouble(mln.softEvidenceLambda * gf.originalWeight.getValue(),true);
        }
    }

    public String writableString(){
        String result = "";
        for(GroundFormula gf : groundFormulas)
        {
            result += gf.writableString(indexToGroundPredMap);
            result += " :: " + gf.weight.getValue() +"\n";
        }
        return result;
    }

    public void setNumConnections()
    {
        Map<Integer, Map<Pair, Double>> connections = new HashMap<>();
        long time = System.currentTimeMillis();
        for(GroundFormula gf : groundFormulas)
        {
            List<Integer> formulaIds = gf.parentFormulaId;
            for(int formulaId : formulaIds)
            {
                if(!connections.containsKey(formulaId))
                {
                    connections.put(formulaId, new HashMap<Pair, Double>());
                }

                for (int i = 0 ; i < gf.allGroundPredIndices.size() ; i++) {
                    int gpIndex = gf.allGroundPredIndices.get(i);
                    Pair key = new Pair(i,gpIndex);
                    if(!connections.get(formulaId).containsKey(key))
                    {
                        connections.get(formulaId).put(key,0.0);
                    }
                    double val = connections.get(formulaId).get(key);
                    connections.get(formulaId).put(key,val+1);
                }
            }
        }
        System.out.println("Time taken ; " + org.utd.cs.gm.utility.Timer.time((System.currentTimeMillis() - time)/1000.0));
        for(GroundFormula gf : groundFormulas)
        {
            List<Integer> formulaIds = gf.parentFormulaId;
            for(int formulaId : formulaIds)
            {
                gf.numConnections.put(formulaId, new ArrayList<Double>());
                for (int i = 0 ; i < gf.allGroundPredIndices.size() ; i++) {
                    int gpIndex = gf.allGroundPredIndices.get(i);
                    Pair key = new Pair(i,gpIndex);
                    gf.numConnections.get(formulaId).add(connections.get(formulaId).get(key));
                }
            }
            //System.out.println("parentId : "+gf.parentFormulaId);
            //System.out.println("numConnections : "+gf.numConnections);
        }
        System.gc();
    }

    public void setEffWts(MLN mln)
    {
        for(GroundFormula gf : groundFormulas) {
            List<Integer> parentFormulaId = gf.parentFormulaId;
            List<Integer> numCopies = gf.numCopies;
            double effWeight = 0.0;
            if (!parentFormulaId.isEmpty()) {
                for (int i = 0; i < parentFormulaId.size(); i++) {
                    int formulaId = parentFormulaId.get(i);
                    List<Double> doubWeight = new ArrayList<>();
                    for (int j = 0; j < gf.numConnections.get(formulaId).size(); j++) {
                        doubWeight.add(mln.formulas.get(formulaId).weight.getValue());
                    }
                    effWeight += numCopies.get(i) * VecOperations.dotprod(doubWeight, Aggregator.aggregator(gf.numConnections.get(formulaId)));
                }
                gf.weight = new LogDouble(effWeight, true);
            }
            else
                gf.weight = new LogDouble(mln.softEvidenceLambda * gf.originalWeight.getValue(), true);
        }
    }
}
