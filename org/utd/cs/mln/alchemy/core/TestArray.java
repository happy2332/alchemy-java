package org.utd.cs.mln.alchemy.core;

import java.util.Arrays;

/**
 * Created by piyush on 22/1/18.
 */
public class TestArray {
    public static void main(String[] args) {
        int[][][] b = new int[4][][];
        b[0] = test();
        b[1] = test();
        b[2] = test();
        b[3] = test();
        System.out.println(Arrays.toString(b[0][0]));
        System.out.println(Arrays.toString(b[0][1]));
        System.out.println(Arrays.toString(b[0][2]));
        System.out.println(Arrays.toString(b[1][0]));
        System.out.println(Arrays.toString(b[1][1]));
        System.out.println(Arrays.toString(b[1][2]));
        System.out.println(Arrays.toString(b[2][0]));
        System.out.println(Arrays.toString(b[2][1]));
        System.out.println(Arrays.toString(b[2][2]));
        System.out.println(Arrays.toString(b[3][0]));
        System.out.println(Arrays.toString(b[3][1]));
        System.out.println(Arrays.toString(b[3][2]));
    }
    private static int[][] test(){
        int[][] a = new int[][]{{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
        return a;
    }
}
