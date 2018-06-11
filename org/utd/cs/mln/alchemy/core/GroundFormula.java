package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;

import java.util.*;

/**
 * Created by Happy on 2/21/17.
 */
public class GroundFormula {
    public List<GroundClause> groundClauses = new ArrayList<>();
    public Set<Integer> groundPredIndices = new HashSet<>();
    public int formulaId; // index of this formula in the MLN's list of ground formulas
    public List<Integer> parentFormulaId = new ArrayList<Integer>(); // id of first order formula from which this came
    public List<Integer> numCopies = new ArrayList<Integer>(); // for each first order formula from which this groundformula is coming, it stores number of copies of it.
    public LogDouble weight;
    public  LogDouble originalWeight; // weight without multiplication with softEvidenceLambda
    public Map<Integer, List<Double>> numConnections = new HashMap<>();
    public List<Integer> allGroundPredIndices = new ArrayList<>();

    @Override
    public String toString() {
        return "GroundFormula{" +
                "\nparentFormulaId=" + parentFormulaId +
                ", \nformulaId=" + formulaId +
                ", \ngroundClauses=" + groundClauses +
                ", \nweight=" + weight +
                '}';
    }

    public String writableString(Map<Integer, GroundPredicate> indexToGroundPredMap){
        String result = "";
        for (int gcId = 0; gcId < groundClauses.size(); gcId++) {
            GroundClause gc = groundClauses.get(gcId);
            result += gc.writableString(indexToGroundPredMap);
            if(gcId < groundClauses.size()-1)
                result += " ^ ";
        }
        return result;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GroundFormula that = (GroundFormula) o;

        return groundClauses.equals(that.groundClauses);

    }

    @Override
    public int hashCode() {
        return groundClauses.hashCode();
    }
}