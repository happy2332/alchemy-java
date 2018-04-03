package org.utd.cs.mln.inference;

/**
 * This class holds parameters common to all MCMC inference algorithms.
 * @author Happy
 * @since 03/28/18
 * @see GibbsSampler_v3,GibbsParams
 */
public class MCMCParams {
    // No. of chains which MCMC will use
    public int numChains;
    // Min. no. of burn-in steps MCMC will take per chain
    public int burnMinSteps;

    @Override
    public String toString() {
        return "MCMCParams{" +
                "numChains=" + numChains +
                ", burnMinSteps=" + burnMinSteps +
                ", burnMaxSteps=" + burnMaxSteps +
                ", minSteps=" + minSteps +
                ", maxSteps=" + maxSteps +
                ", maxSeconds=" + maxSeconds +
                '}';
    }

    // Max. no. of burn-in steps MCMC will take per chain
    int burnMaxSteps;
    // Min. no. of sampling steps MCMC will take per chain
    int minSteps;
    // Max. no. of sampling steps MCMC will take per chain
    int maxSteps;
    // Max. no. of seconds MCMC should run
    int maxSeconds;
}
