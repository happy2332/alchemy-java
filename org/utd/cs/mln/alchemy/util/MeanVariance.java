package org.utd.cs.mln.alchemy.util;

/**
 * Created by Happy on 3/29/18.
 */
public class MeanVariance {
    public double totalX_;
    public double totalXSquared_;
    public int numValues_;

    public MeanVariance()
    {
        totalX_ = 0;
        totalXSquared_ = 0;
        numValues_ = 0;
    }

    public void reset()
    {
        totalX_ = 0;
        totalXSquared_ = 0;
        numValues_ = 0;
    }

    public void appendValue(double x)
    {
        totalX_ += x;
        totalXSquared_ += x*x;
        numValues_++;
    }

    public double getMean()
    {
        return totalX_ / numValues_;
    }

    public double getVariance()
    {
        // mean(xsquared) - (mean(X))^2
        double meanXSquared = totalXSquared_/numValues_;
        double meanX = totalX_/numValues_;
        return (((double)numValues_)/(numValues_-1)) * (meanXSquared - meanX*meanX);
    }
}
