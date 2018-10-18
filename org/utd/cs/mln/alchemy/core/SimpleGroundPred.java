package org.utd.cs.mln.alchemy.core;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Happy on 8/16/18.
 */
public class SimpleGroundPred {
    public String symbolName;
    public List<Integer> terms;
    public SimpleGroundPred()
    {
        symbolName = new String();
        terms = new ArrayList<>();
    }
}
