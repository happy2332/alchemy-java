package org.utd.cs.mln.alchemy.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.utd.cs.mln.alchemy.util.OtherUtils.getUniformAssignment;

/**
 * Created by Happy on 8/16/18.
 */
public class SimpleState {
    // This state is associated to a particular SimpleGroundMLN
    public SimpleGroundMln simpleGroundMln;
    public List<Integer> world;
    public SimpleState(SimpleGroundMln simpleGroundMln_)
    {
        simpleGroundMln = simpleGroundMln_;
        world = new ArrayList<>(Collections.nCopies(simpleGroundMln.groundPreds.size(),0));
    }

    public void randomize() {
        for (int i = 1; i < world.size(); i++) {
            int assignment = getUniformAssignment(2);
            world.set(i, assignment);
        }
    }
}
