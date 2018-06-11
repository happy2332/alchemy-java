package org.utd.cs.mln.learning;

/**
 * Created by Happy on 4/21/18.
 */
public class CGParams {
    public double cg_lambda;
    public double cg_max_lambda;
    public int maxBackTracks;
    public boolean preConditionCG;
    //max number of learning iterations
    public int numIter;
    public double min_ll_change;

    public CGParams(CGParams cgParam) {
        this.cg_lambda = cgParam.cg_lambda;
        this.cg_max_lambda = cgParam.cg_max_lambda;
        this.maxBackTracks = cgParam.maxBackTracks;

        this.preConditionCG = cgParam.preConditionCG;
        this.numIter = cgParam.numIter;
        this.min_ll_change = cgParam.min_ll_change;
    }

    public CGParams() {
        this.cg_lambda = 100.0;
        this.cg_max_lambda = Double.MAX_VALUE;
        this.maxBackTracks = 1000;

        this.preConditionCG = true;
        this.numIter = 100;
        this.min_ll_change = 0.001;
    }
}
