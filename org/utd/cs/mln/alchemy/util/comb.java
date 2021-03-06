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
		Set<Integer> a = new HashSet<>();
		a.add(25);
		a.add(50);
		Set<Integer>b = (a);
		b.remove(25);
		System.out.println("a = " + a);
		MyAssert.assume(a.size() == 0);
	}

}
