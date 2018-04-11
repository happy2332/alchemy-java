package org.utd.cs.mln.learning;

import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.inference.GibbsParams;
import org.utd.cs.mln.inference.GibbsSampler_v3;
import org.utd.cs.mln.inference.Inference;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 * Created by Happy on 4/11/18.
 */
public class DiscLearn extends WeightLearner {
    private List<Inference> inferences;
    private double[] gradient;
    private long time;

    // for each domain d, for each first order formula f, stores number of true groundings according to domain d.
    // formulaTrainCnts[d][f]
    private double[][] formulaTrainCnts;

    // for each domain d, for each first order formula f, stores number of total groundings according to domain d.
    // formulaNumGndings[d][f]
    private double[][] formulaNumGroundings;

    /**
     * Constructor : Initializes the fields
     *
     * @param mlnsParam List of MLNs of size number of databases. All MLNs are same except for domains for predicates.
     * @param lArgs
     */
    public DiscLearn(List<MLN> mlnsParam, List<GroundMLN> groundMlnsParam, List<GroundMLN> groundMlnsEMParam,
                     List<Evidence> truthsParam, List<Evidence> truthsEMParam, LearnArgs lArgs) {
        super(mlnsParam, lArgs);
        init(groundMlnsParam, truthsParam, groundMlnsEMParam, truthsEMParam, lArgs);
    }

    private void init(List<GroundMLN> groundMlnsParam, List<Evidence> truthsParam, List<GroundMLN> groundMlnsEMParam, List<Evidence> truthsEMParam, LearnArgs lArgs)
    {
        inferences = new ArrayList<>();
        for (int i = 0; i < domain_cnt; i++) {
            State state = new State(mlns.get(i), groundMlnsParam.get(i), truthsParam.get(i));
            GibbsParams gibbsparams = new GibbsParams();
            Inference inference = new GibbsSampler_v3(state,-1,true, priorSoftEvidence, gibbsparams);
            inferences.add(inference);
        }
        gradient = new double[numWts];
    }

    @Override
    public void learnWeights() {
        time = System.currentTimeMillis();
        if(!withEM)
            findFormulaTrainCnts();
    }

    private void findFormulaTrainCnts() {
        for (int i = 0; i < domain_cnt; i++) {
            GroundMLN gm = inferences.get(i).state.groundMLN;
            for(GroundFormula gf : gm.groundFormulas)
            {
                List<Integer> parentFormulaId = gf.parentFormulaId;
                List<Integer> numCopies = gf.numCopies;
                if(parentFormulaId.get(0) == -1)
                    continue;
                boolean isFormulaSatisfied = true;
                for(GroundClause gc : gf.groundClauses)
                {
                    boolean isClauseSatisfied = false;
                    for(int gpId : gc.groundPredIndices)
                    {
                        int trueVal = inferences.get(i).state.truthVals.get(gpId);
                        BitSet b = gc.grounPredBitSet.get(gc.globalToLocalPredIndex.get(gpId));
                        isClauseSatisfied |= b.get(trueVal);
                        if(isClauseSatisfied) {
                            break;
                        }
                    }
                    isFormulaSatisfied &= isClauseSatisfied;
                    if(!isFormulaSatisfied) {
                        break;
                    }
                }

                for (int j = 0; j < parentFormulaId.size(); j++) {
                    if(isFormulaSatisfied) {
                        formulaTrainCnts[i][parentFormulaId.get(j)] += numCopies.get(j);
                    }
                    formulaNumGroundings[i][parentFormulaId.get(j)] += numCopies.get(j);
                }

            }
        }
    }
}
