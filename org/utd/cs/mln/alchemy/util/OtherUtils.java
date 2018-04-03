package org.utd.cs.mln.alchemy.util;

import java.util.*;

/**
 * Created by Happy on 3/27/18.
 */
public class OtherUtils {
    public static List<List<Integer>> cartesianProd(List<Set<Integer>> sets){
        List<List<Integer>> result = new ArrayList<List<Integer>>();
        int numTuples = 1;
        for(int i = 0 ; i < sets.size() ; i++){
            numTuples *= sets.get(i).size();
        }

        // add numTuples lists into result
        for(int i = 0 ; i < numTuples ; i++){
            result.add(new ArrayList<Integer>());
        }

        int numRepeats = numTuples;
        // Now fill each element of input sets
        for(int i = 0 ; i < sets.size() ; i++){
            numRepeats = numRepeats/sets.get(i).size();
            int numFilledEntries = 0;
            while(numFilledEntries != numTuples){
                for(Integer elem : sets.get(i)){
                    for(int j = 0 ; j < numRepeats ; j++){
                        result.get(numFilledEntries).add(elem);
                        numFilledEntries++;
                    }
                }
            }
        }
        return result;
    }

    // Given a list of double values, normalize them into probability distribution

    public static List<Double> getProbabilityDistribution(List<Double> vals){
        List<Double> probabilities = new ArrayList<Double>();
        double maxVal = Collections.max(vals);
        //First calculate sum of all exponentiated satWeights
        double sum = 0.0;
        for(Double val : vals) {
            sum += Math.exp(val-maxVal);
        }
        // Now calculate probabilities
        for(Double val : vals) {
            probabilities.add(Math.exp(val-maxVal)/sum);
        }
        return probabilities;
    }
    public static void main(String args[])
    {
        List<Set<Integer>> myList = new ArrayList<>();
        myList.add(new HashSet<>(Arrays.asList(1,2)));
        myList.add(new HashSet<>(Arrays.asList(3,4)));
        myList.add(new HashSet<>(Arrays.asList(5,6)));
        System.out.println(OtherUtils.cartesianProd(myList));
    }
}

