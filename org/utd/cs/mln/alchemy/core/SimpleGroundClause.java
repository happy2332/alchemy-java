package org.utd.cs.mln.alchemy.core;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Happy on 8/16/18.
 */
public class SimpleGroundClause {
    public double weight;
    public List<Integer> groundAtomIndices; // List of indices of groundAtoms. Ex : If groundMln.groundPreds = [NA, S(A), S(B), C(A)],
    // and if groundAtomIndices = [-1,2], then it means !S(A) | S(B)
    public List<Integer> parentFormulaIds;
    public List<Integer> numCopies;
    public SimpleGroundClause()
    {
        groundAtomIndices = new ArrayList<>();
        parentFormulaIds = new ArrayList<>();
        numCopies = new ArrayList<>();
    }
}
