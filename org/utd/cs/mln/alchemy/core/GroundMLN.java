package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;

import java.util.*;

/**
 * Created by Happy on 2/28/17.
 */
public class GroundMLN {
    public List<GroundPredicateSymbol> symbols = new ArrayList<>();
    public List<GroundPredicate> groundPredicates = new ArrayList();
    public Map<GroundPredicate,Integer> groundPredToIntegerMap = new HashMap<>();
    public List<GroundFormula> groundFormulas = new ArrayList<>();

    /**
     * This function fills in {@link GroundMLN#groundPredToIntegerMap} which creates a mapping from groundpredicate to its index in {@link GroundMLN#groundPredicates} list.
     * So if groundPredicates[i] = gp, then groundPredToIntegerMap[gp] = i.
     * This is useful to search a groundpredicate gp in groundPredicates list.
     */
    public void setGroundPredToIntegerMap()
    {
        for (int i = 0; i < groundPredicates.size(); i++) {
            groundPredToIntegerMap.put(groundPredicates.get(i),i);
        }
    }

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

    public void setGroundFormulaWtsToParentWts(MLN mln) {
        for(GroundFormula gf : groundFormulas)
        {
            List<Integer> parentFormulaId = gf.parentFormulaId;
            List<Integer> numCopies = gf.numCopies;
            double weight = 0.0;
            for (int i = 0; i < parentFormulaId.size(); i++) {
                weight += numCopies.get(i) * mln.formulas.get(parentFormulaId.get(i)).weight.getValue();
            }
            gf.weight = new LogDouble(weight, true);
        }
    }

    public void setGroundFormulaWtsToParentWtsSoftEvidence(MLN mln, double lambda) {
        for(GroundFormula gf : groundFormulas)
        {
            List<Integer> parentFormulaId = gf.parentFormulaId;
            List<Integer> numCopies = gf.numCopies;
            if(parentFormulaId.get(0) != -1)
            {
                double weight = 0.0;
                for (int i = 0; i < parentFormulaId.size(); i++) {
                    weight += numCopies.get(i) * mln.formulas.get(parentFormulaId.get(i)).weight.getValue();
                }
            }
            else
                gf.weight = new LogDouble(lambda*gf.originalWeight.getValue(),true);
        }
    }

    public String writableString(){
        String result = "";
        for(GroundFormula gf : groundFormulas)
        {
            result += gf.writableString(groundPredicates);
            result += " :: " + gf.weight.getValue() +"\n";
        }
        return result;
    }
}
