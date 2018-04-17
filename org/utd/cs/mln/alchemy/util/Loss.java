package org.utd.cs.mln.alchemy.util;

import java.util.Arrays;

/**
 * This class is an abstract class, which defines loss associated with a state of an MLN.
 * @author Happy
 * @since 04/17/18
 */
public abstract class Loss {
    double loss_value;
    public abstract double getLossValue();
    public abstract double[] getGradient();
}
