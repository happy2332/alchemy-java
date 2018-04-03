package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.mln.lmap.HyperCubeRoot;

import java.util.ArrayList;
import java.util.List;

public class Formula {

	public List<WClause> clauses = new ArrayList<>();
	public HyperCubeRoot root = new HyperCubeRoot();
	public List<Term> terms = new ArrayList<Term>();
	public boolean isEvidence;
	public LogDouble weight;
	public int formulaId;
	public boolean tautology;

	public Formula(List<WClause> clauses_, LogDouble weight_, boolean tautology) {
		this(clauses_, weight_, false, tautology);
	}

	public Formula(List<WClause> clauses_, LogDouble weight_, boolean isEvidence_, boolean tautology) {
		clauses = clauses_;
		weight = (weight_);
		isEvidence = (isEvidence_);
		this.tautology = tautology;
	}

}
