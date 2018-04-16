package org.utd.cs.mln.inference;

/**
 * Created by Happy on 3/29/18.
 * This class holds parameters needed to run Gibbs sampling.
 * @see MCMCParams
 */
public class GibbsParams extends MCMCParams{
    public double  gamma;
    public double  epsilonError;
    double  fracConverged;
    public boolean testConvergence;
    public int     samplesPerTest;

    @Override
    public String toString() {
        return "GibbsParams{" +
                "gamma=" + gamma +
                ", epsilonError=" + epsilonError +
                ", fracConverged=" + fracConverged +
                ", testConvergence=" + testConvergence +
                ", samplesPerTest=" + samplesPerTest +
                "} " + super.toString();
    }

    public GibbsParams() {
        numChains       = 5;
        burnMinSteps    = 100;
        burnMaxSteps    = 100;
        minSteps        = -1;
        maxSteps        = 2000;
        maxSeconds      = -1;
        gamma           = 1 - 0.05;;
        epsilonError    = 0.01;
        fracConverged   = 0.95;
        testConvergence = true;
        samplesPerTest  = 100;
    }

}
