package org.utd.cs.mln.alchemy.util;

import java.util.Arrays;

/**
 * This class is an abstract class, which defines loss associated with a state of an MLN.
 * @author Happy
 * @since 04/17/18
 */
public abstract class Loss {
    // priorLambda is lambda for weight regularization. For each weight, there can be different lambda. If an entry of prior
    // is zero, it means there is no regularization on that corresponding weight.
    // priorMeans and priorStdDevs are the mean and standard deviation vectors for gaussian prior. Note that standard
    // deviation is a vector here, meaning its a diagonal matrix.
    // By default, priorMeans = 0, and priorStdDevs = 1, in that case, gaussian prior is just l2 regularization.
    public double [] priorLambda, priorMeans, priorStdDevs;
    public abstract double getLossValue();

    public abstract double[] getGradient();

    public void init(int numWts) {
        priorMeans = new double[numWts];
        priorStdDevs = new double[numWts];
        Arrays.fill(priorStdDevs,2.0);
        priorLambda = new double[numWts];
        Arrays.fill(priorLambda,0.0);
    }
}
