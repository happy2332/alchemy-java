package org.utd.cs.mln.alchemy.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class comb {
	
	public static double findComb(int n, int k){
		return Math.exp(fact(n) - fact(n-k) - fact(k));
	}
	public static double fact(int n){
		double result = 0.0;
		for(int i = 1 ; i <= n ; i++){
			result += Math.log(i);
		}
		return result;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.println(findComb(10,5));
		List<Integer> a = new ArrayList<>();
		a.add(5);
		Set<Integer> b = new HashSet<>(a);
		b.clear();
		System.out.println("b = " + b);

	}

}
