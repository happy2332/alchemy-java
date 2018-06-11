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
    // Max. no. of burn-in steps MCMC will take per chain
    public int burnMaxSteps;
    // Min. no. of sampling steps MCMC will take per chain
    public int minSteps;
    // Max. no. of sampling steps MCMC will take per chain
    public int maxSteps;
    // Max. no. of seconds MCMC should run
    public int maxSeconds;

    public MCMCParams(MCMCParams other) {
        this.numChains = other.numChains;
        this.burnMinSteps = other.burnMinSteps;
        this.burnMaxSteps = other.burnMaxSteps;
        this.minSteps = other.minSteps;
        this.maxSteps = other.maxSteps;
        this.maxSeconds = other.maxSeconds;
    }

    public MCMCParams() {
        numChains       = 10;
        burnMinSteps    = 100;
        burnMaxSteps    = 100;
        minSteps        = -1;
        maxSteps        = 1000;
        maxSeconds      = -1;
    }

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
}
