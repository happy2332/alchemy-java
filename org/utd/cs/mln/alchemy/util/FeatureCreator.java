package org.utd.cs.mln.alchemy.util;

import com.beust.jcommander.ParameterException;
import org.utd.cs.mln.alchemy.core.*;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * This class creates feature vector for all the constants of a particular type, which can then be used to
 * cluster the constants.
 * For a constant c, the size of a feature vector fvec is the number of first order formulas in which
 * c can occur.
 * fvec[fo] is the number of satisfied groundings of formula fo, when c occurs in fo.
 * To calculate features, this class requires an MLN file, a database, and variable type name (only this type's
 * constants' features will be calculated)
 * @author Happy
 * @since 04/01/18
 */
public class FeatureCreator {
    public static void main(String []args) throws ParameterException, FileNotFoundException, PredicateNotFound {
        // args must be of length 3
        if(args.length != 4)
            throw new ParameterException("Wrong arguments !!!");
        String mlnFile = args[0];
        String dbFile = args[1];
        String outFile = args[2];
        String typeName = args[3];
        MLN mln = new MLN();
        Parser parser = new Parser(mln);
        parser.parseInputMLNFile(mlnFile);
        List<String> files = new ArrayList<>();
        files.add(dbFile);
        parser.collectDomain(files);
        mln.overWriteDomain();
        boolean closedWorld = true;
        GroundMLN groundMLN = createGroundMLN(mln, dbFile, mln.varTypeToDomainMap);
        Evidence truth = parser.parseEvidence(groundMLN,dbFile);
        Map<Integer, List<Integer>> features = createFeatures(mln, groundMLN, truth, typeName, mln.varTypeToDomainMap.get(typeName), closedWorld);
        writeFeatures(features, outFile);
    }

    private static void writeFeatures(Map<Integer, List<Integer>> features, String outFile) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(outFile);
        for(int constant : features.keySet())
        {
            for(int val : features.get(constant))
            {
                pw.write(val+",");
            }
            pw.write(constant+"\n");
        }
        pw.close();
    }

    private static GroundMLN createGroundMLN(MLN mln, String dbFile, Map<String, Set<Integer>> varTypeToDomain) {
        // This groundMLN will contain only list of groundPredicates, not groundformulas since we are not grounding MLN.
        GroundMLN groundMLN = new GroundMLN();
        List<GroundPredicate> groundPredicates = FullyGrindingMill.createGroundPredicates(mln);
        groundMLN.groundPredicates = groundPredicates;
        groundMLN.setGroundPredToIntegerMap();
        return groundMLN;
    }

    private static Map<Integer,List<Integer>> createFeatures(MLN mln, GroundMLN groundMLN, Evidence truth, String typeName, Set<Integer> constants, boolean closedWorld) {
        Map<Integer, List<Integer>> featureVectors = new HashMap<>();
        int numFormulas = mln.formulas.size();
        for(int constant : constants)
        {
            featureVectors.put(constant, new ArrayList<>(Collections.nCopies(numFormulas,0)));
        }
        // for each formula, compute sat counts
        for (int formulaId = 0; formulaId < mln.formulas.size(); formulaId++) {
            Formula formula = mln.formulas.get(formulaId);
            Set<Term> formulaWiseTermToGround = new HashSet<Term>();
            for (WClause clause : formula.clauses) {
                for (Atom atom : clause.atoms) {
                    for (int j = 0; j < atom.terms.size(); j++) {
                        Term term = atom.terms.get(j);
                        formulaWiseTermToGround.add(term);
                    }
                }
            }
            createFeaturesPerFormula(mln, groundMLN, formulaId, new ArrayList<Term>(formulaWiseTermToGround), truth, typeName, closedWorld, featureVectors);
        }
        return featureVectors;
    }

    private static void createFeaturesPerFormula(MLN mln, GroundMLN groundMLN, int formulaId, ArrayList<Term> terms,
                                                 Evidence truth, String typeName, boolean closedWorld, Map<Integer, List<Integer>> featureVectors) {
        int[][] permutations = FullyGrindingMill.permute(terms);
        Formula formula = mln.formulas.get(formulaId);
        for (int i = 0; i < permutations.length; i++) {
            GroundFormula newFormula = new GroundFormula();
            List<GroundClause> newGroundClauseList = new ArrayList<GroundClause>();
            Set<Integer> formulaGpIndices = new HashSet<>();
            for (int c = 0; c < formula.clauses.size(); c++) {
                WClause clause = formula.clauses.get(c);
                GroundClause newGroundClause = new GroundClause();
                Map<Integer, BitSet> gpIndexToSatVals = new HashMap<>();
                // Iterate over each first order atom, and create ground atom for it.
                for (int j = 0; j < clause.atoms.size(); j++) {
                    boolean sign = clause.sign.get(j);
                    Atom oldAtom = clause.atoms.get(j); // first order atom
                    int valTrue = clause.valTrue.get(j);
                    GroundPredicate gp = new GroundPredicate(); // GroundPredicate to create
                    gp.symbol = new GroundPredicateSymbol(oldAtom.symbol.id, oldAtom.symbol.symbol, oldAtom.symbol.values, oldAtom.symbol.variable_types);
                    // Fill in the terms with constants
                    for (Term term : oldAtom.terms) {
                        int termIndex = terms.indexOf(term);
                        gp.terms.add(permutations[i][termIndex]);
                    }

                    int numPossibleValues = oldAtom.symbol.values.values.size();
                    gp.numPossibleValues = numPossibleValues;
                    int gpIndex = groundMLN.groundPredToIntegerMap.get(gp);
                    formulaGpIndices.add(gpIndex);

                    // Check if this groundPredicate occurs first time in this ground clause. then update
                    // groundClause's data structures about this groundPredicate.
                    int gpIndexInClause = newGroundClause.groundPredIndices.indexOf(gpIndex);
                    //GroundAtom newGroundAtom = new GroundAtom(gpIndex, gpIndexInClause, valTrue, sign);
                    if (gpIndexInClause == -1) {
                        newGroundClause.groundPredIndices.add(gpIndex);
                        gpIndexInClause = newGroundClause.groundPredIndices.size() - 1;
                        newGroundClause.globalToLocalPredIndex.put(gpIndex, gpIndexInClause);
                        gpIndexToSatVals.put(gpIndexInClause, new BitSet(gp.numPossibleValues));
                    }

                    // Now once we have added new ground Atom, we need to check if ground clause gets satisfied or not.
                    BitSet gpBitSet = new BitSet(gp.numPossibleValues);
                    gpBitSet.set(valTrue);
                    if (sign == true)
                        gpBitSet.flip(0, gp.numPossibleValues);
                    gpBitSet.or(gpIndexToSatVals.get(gpIndexInClause));

                    // If all bits are set for this groundPred, then this clause will always be satisfied and hence,
                    // shouldn't be added into groundformula. Note that, although at this point, we know that
                    // this clause shouldn't be added, but still we shouldn't just break out of this loop, as we
                    // need to add groundPredicates, but we shouldn't add any clauseInfo into groundPredicates appearing
                    // in this clause.

                    gpIndexToSatVals.put(gpIndexInClause, gpBitSet);
                }
                for (int gpId = 0; gpId < newGroundClause.groundPredIndices.size(); gpId++) {
                    BitSet b = gpIndexToSatVals.get(gpId);
                    newGroundClause.grounPredBitSet.add(b);
                }
                newGroundClauseList.add(newGroundClause);
            }
            newFormula.groundClauses.addAll(newGroundClauseList);

            boolean isFormulaSatisfied = true;
            for(GroundClause gc : newFormula.groundClauses)
            {
                boolean isClauseSatisfied = false;
                for(Integer gpIndex : gc.groundPredIndices)
                {
                    int truthVal = 0;
                    if(truth.predIdVal.containsKey(gpIndex))
                    {
                        truthVal = truth.predIdVal.get(gpIndex);
                    }
                    else
                    {
                        if(!closedWorld)
                            continue;
                    }
                    int localIndex = gc.globalToLocalPredIndex.get(gpIndex);
                    isClauseSatisfied = gc.grounPredBitSet.get(localIndex).get(truthVal);
                    if(isClauseSatisfied)
                        break;
                }
                if(!isClauseSatisfied)
                {
                    isFormulaSatisfied = false;
                    break;
                }
            }
            // If formula is satisfied then update feature vector of each constant appearing in this formula
            if(isFormulaSatisfied)
            {
                for (GroundClause gc : newFormula.groundClauses)
                {
                    for (Integer gpIndex : gc.groundPredIndices)
                    {
                        GroundPredicate gp = groundMLN.groundPredicates.get(gpIndex);
                        for (int j = 0; j < gp.terms.size(); j++) {
                            if(gp.symbol.variable_types.get(j).equals(typeName))
                            {
                                int constant = gp.terms.get(j);
                                int numSatisfied = featureVectors.get(constant).get(formulaId);
                                featureVectors.get(constant).set(formulaId, numSatisfied+1);
                            }
                        }
                    }

                }
            }
        }
    }

}
