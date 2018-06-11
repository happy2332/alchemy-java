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

    public GibbsParams(GibbsParams other) {
        super(other);
        this.gamma = other.gamma;
        this.epsilonError = other.epsilonError;
        this.fracConverged = other.fracConverged;
        this.testConvergence = other.testConvergence;
        this.samplesPerTest = other.samplesPerTest;
    }

    public GibbsParams() {
        super();
        gamma           = 1 - 0.05;
        epsilonError    = 0.01;
        fracConverged   = 0.95;
        testConvergence = true;
        samplesPerTest  = 10;
    }

}
