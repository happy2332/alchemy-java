package org.utd.cs.mln.alchemy.core;

import org.utd.cs.mln.alchemy.util.MyAssert;

import java.util.*;

/**
 * Created by Happy on 8/16/18.
 */
public class SimpleGroundMln {
    // It is associated to particular MLN
    public MLN mln;
    public List<SimpleGroundPred> groundPreds; // List starts from index 1
    public List<SimpleGroundClause> groundClauses;
    public SimpleGroundMln(MLN mln_)
    {
        mln = mln_;
        groundClauses = new ArrayList<>();
        groundPreds = new ArrayList<>();
    }
    public void createFromGroundMLN(GroundMLN groundMLN)
    {
        // First copy and create all ground preds
        int numGroundPreds = groundMLN.groundPredToIntegerMap.size();
        // First create 1 empty groundPred
        groundPreds.add(null);

        // Now create a temporary map from old groundPredicates to new Continous indices
        Map<GroundPredicate, Integer> tempMap = new HashMap<>();
        int predIndex = 1;
        for(GroundPredicate gp : groundMLN.groundPredToIntegerMap.keySet())
        {
            tempMap.put(gp,predIndex++);
            SimpleGroundPred sgp = new SimpleGroundPred();
            sgp.symbolName = new String(gp.symbol.symbol);
            sgp.terms = new ArrayList<>(gp.terms);
            groundPreds.add(sgp);
        }

        // Now create simpleGroundClauses
        for(GroundFormula gf : groundMLN.groundFormulas)
        {
            GroundClause gc = gf.groundClauses.get(0);
            SimpleGroundClause sgc = new SimpleGroundClause();
            sgc.weight = gf.weight.getValue();

            // Fill in sgc.groundAtomIndices
            for(Integer gpId : gc.groundPredIndices)
            {
                BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                MyAssert.assume(b.size()==2 && b.get(0) != b.get(1));
                // Check if this atom is with sign or not
                int sign = 1;
                if(b.get(0) == true)
                {
                    sign = -1;
                }
                // Find new index
                int index = tempMap.get(groundMLN.indexToGroundPredMap.get(gpId));
                // new gpId is +1
                sgc.groundAtomIndices.add(sign*index);
            }

            // Fill in parentFormulaIds and numCopies
            sgc.parentFormulaIds = new ArrayList<>(gf.parentFormulaId);
            sgc.numCopies = new ArrayList<>(gf.numCopies);
        }
    }
}
