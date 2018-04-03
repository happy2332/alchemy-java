package org.utd.cs.mln.alchemy.core;

import java.util.ArrayList;

public class FirstEvidence {
    public FirstEvidence(PredicateSymbol symbol, ArrayList<Integer> terms, int truthValue) {
        super();
        this.symbol = symbol;
        this.terms = terms;
        this.truthValue = truthValue;
    }
    public PredicateSymbol symbol;
    public ArrayList<Integer> terms;
    public int truthValue;
}


