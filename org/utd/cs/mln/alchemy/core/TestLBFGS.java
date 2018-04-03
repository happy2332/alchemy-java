package org.utd.cs.mln.alchemy.core;

import util.*;
import util.LBFGS.Function;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * Created by piyush on 19/1/18.
 */
public class TestLBFGS {
    static void linreg_test(final double[][] X, final double[] Y) {
        final int Nfeat = X[0].length;
        Function f = new Function() {
            @Override
            public double evaluate(double[] beta, double[] g, int n, double step) {
                double sqloss = 0;
                Arr.fill(g, 0); // is this necessary?

                for (int i=0; i<X.length; i++) {
                    double pred = Arr.innerProduct(beta, X[i]);
                    double resid = pred - Y[i];
                    sqloss += Math.pow(resid, 2);
                    for (int j=0; j<Nfeat; j++) {
                        g[j] += X[i][j] * resid;
                    }
                }
//                System.out.println("n: "+n);
//                Arrays.fill(g, 1);
                System.out.println("g: "+ Arrays.toString(g));

                return sqloss;
            }
        };
        LBFGS.Params p = new LBFGS.Params();
        LBFGS.ProgressCallback cb = new LBFGS.ProgressCallback() {
            @Override
            public int apply(double[] x, double[] g, double fx, double xnorm,
                             double gnorm, double step, int n, int k, LBFGS.Status ls) {
                U.pf("ITER %d obj=%g sol=%.6g\n", k, fx, x[0]);
                return 0;
            }
        };
        double[] coef = new double[Nfeat];
        LBFGS.Result r = LBFGS.lbfgs(coef, f, cb, p);
        U.p(r);
        U.p(coef);
    }


    static void mean_test() {
        final int target = 3;
        Function f = new Function() {
            @Override
            public double evaluate(double[] x, double[] g, int n, double step) {
                double resid = x[0] - target;
                double sqloss = Math.pow(resid, 2);
                g[0] = resid;
                return sqloss;
            }
        };
        LBFGS.Params p = new LBFGS.Params();
        LBFGS.ProgressCallback cb = new LBFGS.ProgressCallback() {
            @Override
            public int apply(double[] x, double[] g, double fx, double xnorm,
                             double gnorm, double step, int n, int k, LBFGS.Status ls) {
                U.pf("ITER %d sol=%s obj=%g\n", k, Arr.toString(x), fx);
                return 0;
            }
        };
        double[] sol = new double[1];
        U.p(sol);
        LBFGS.Result r = LBFGS.lbfgs(sol, f, cb, p);
        U.p(r);
        U.p(sol);
    }

    public static void main(String[] args) {
//		mean_test();
        double[][] X = Arr.readDoubleMatrix(args[0]);
        double[] Y = Arr.readDoubleVector(args[1]);
        System.out.println("X.size(): "+X.length);
        System.out.println("X[0].size(): "+X[0].length);
        System.out.print("X: ["+X[0][1]);

        for (int i = 1; i < X.length; ++i){
            System.out.print(", "+X[i][1]);
        }
        System.out.println("]");
        System.out.println("Y: "+Arrays.toString(Y));
        linreg_test(X,Y);
    }
}
