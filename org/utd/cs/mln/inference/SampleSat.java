package org.utd.cs.mln.inference;

/**
 * This class implements SampleSat algorithm, which samples a satisfying assignment uninformly (given a list of
 * clauses)
 * The algo works as follows :
 * Start with a random state
 * while(not all clauses satisfied):
 *  with prob p, make walksat move, otherwise make Simulated Annealing move (SA Move) (we will keep p=0.8)
 * return state
 * Walksat move :
 *  Pick an unsatisfied clause randomly, and flip one of its literal (which makes max number of clauses satisfiable)
 * SA move :
 *  Pick a literal (ground pred) randomly. Find its delta, where delta is difference of satisfied clauses if we flip
 *  this literal and satisfied clauses if don't flip this literal.
 *  If delta > 0:
 *      Flip the literal with prob 1.0
 *  else:
 *      Flip the literal with prob exp(-delta/temp) (we set temp=14)
 *
 * Created by Happy on 5/28/18.
 */
public class SampleSat {
}
