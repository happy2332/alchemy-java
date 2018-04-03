package org.utd.cs.mln.alchemy.util;
import java.util.ArrayList;

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
}
