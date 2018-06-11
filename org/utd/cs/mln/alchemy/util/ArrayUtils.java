package org.utd.cs.mln.alchemy.util;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by Happy on 2/28/18.
 */
public class ArrayUtils {

    public static void bubbleSort(ArrayList<Integer> values, int[] indices) {
        for(int k = 0; k < values.size(); k++) {
            for(int l = k + 1; l < values.size(); l++) {
                if(values.get(k) > values.get(l)) {
                    int temp_value = values.get(k);
                    values.set(k, values.get(l));
                    values.set(l, temp_value);
                    int temp_index = indices[k];
                    indices[k] = indices[l];
                    indices[l] = temp_index;
                }
            }
        }
    }

    public static double dotprod(double[] v1, double[] v2) {
        double total = 0.0;
        MyAssert.assume(v1.length == v2.length);
        for (int i = 0; i < v1.length ; i++) {
            total += v1[i]*v2[i];
        }
        return total;
    }

    public static double getMean(double [] v)
    {
        double sum = 0;
        for (int i = 0; i < v.length; i++) {
            sum += v[i];
        }
        return sum/v.length;
    }

    public static double getVariance(double [] v)
    {
        double [] vSq = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            vSq[i] = v[i]*v[i];
        }
        double E_v = getMean(v);
        return getMean(vSq) - E_v*E_v;
    }

}
