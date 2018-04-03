package org.utd.cs.mln.inference;

import org.utd.cs.mln.alchemy.util.MeanVariance;

/**
 * Tests convrgence of MCMC algorithm by Gelman-Rubin algorithm.
 * <p>
 *     It compares between-chains and within-chains variances, and if difference is small, then algorithm has convereged.
 * </p>
 */
public class GelmanConvergenceTest {
    int numChains_;
    MeanVariance withinChainMeanVars_[];
    int numSamples_;
    MeanVariance betweenChainsMeanVar_;

    GelmanConvergenceTest(int numChains)
    {
        numChains_ = numChains;
        withinChainMeanVars_ = new MeanVariance[numChains_];
        for (int c = 0; c < numChains; c++) {
            withinChainMeanVars_[c] = new MeanVariance();
        }
        betweenChainsMeanVar_ = new MeanVariance();
        numSamples_ = 0;
    }

    // values is array of size numChains_
    void appendNewValues(double []values)
    {
        for (int i = 0; i < numChains_; i++)
            withinChainMeanVars_[i].appendValue(values[i]);
        numSamples_++;
    }

    double getConvergenceScore()
    {
        betweenChainsMeanVar_.reset();
        double totalW = 0;
        for (int i = 0; i < numChains_; i++)
        {
            betweenChainsMeanVar_.appendValue( withinChainMeanVars_[i].getMean() );
            totalW += withinChainMeanVars_[i].getVariance();
        }
        int numValues = withinChainMeanVars_[0].numValues_;

        double B = betweenChainsMeanVar_.getVariance() * numValues;
        double W = totalW / numChains_;

        // score as stated in "Probability and Statistics", DeGroot and Schervish
        double score = B/W;
        return score;
    }

    public static boolean checkConvergenceOfAll(GelmanConvergenceTest tests[], int numTests, boolean print)
    {
        // threshold as stated in "Probability and Statistics", DeGroot & Schervish
        double threshold = 1 + 0.44 * tests[0].numSamples_;
        int maxItem     = -1;
        double maxScore = -1;
        int numbad      = 0;

        for (int f = 0; f < numTests; f++)
        {
            double score = tests[f].getConvergenceScore();

            if (Double.isInfinite(score) || Double.isNaN(score)) { numbad++; continue; }

            if (score > threshold)
            {
                if (print)
                    System.out.println(" Item " + f + "'s score of " + score + " > threshold of "
                            + threshold);
                return false;
            }

            if (score > maxScore)
            {
                maxScore = score;
                maxItem = f;
            }
        }

        if (numbad == numTests)
        {
            if (print) System.out.println(" All scores were inf or Nan!");
            return false;
        }

        // at this point all scores are less than the threshold

        if (print)
            System.out.println(" max item is " + maxItem + " with score " + maxScore
                    + " < threshold of " + threshold);

        return true;
    }
}
