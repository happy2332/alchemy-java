package org.utd.cs.mln.alchemy.core;

import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.mln.lmap.HyperCubeRoot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Formula {

	public List<WClause> clauses = new ArrayList<>();
	public HyperCubeRoot root = new HyperCubeRoot();
	public List<Term> terms = new ArrayList<Term>();
	public boolean isEvidence;
	public LogDouble weight;
	public int formulaId;
	public int parentFormulaId = -1; // In case this is a formula of hyperCubeMLN, then store formulaId of first order formula from which it is generated.
	public boolean tautology;

	public Formula(List<WClause> clauses_, LogDouble weight_, boolean tautology) {
		this(clauses_, weight_, false, tautology);
	}

	@Override
	public String toString() {
		Map<Term, String> termReferenceMap = new HashMap<>();
		String result = "";
		for (int i = 0; i < clauses.size(); i++) {
			WClause clause = clauses.get(i);
			for (int j = 0; j < clause.atoms.size(); j++) {
				Atom atom = clause.atoms.get(j);
				if(clause.sign.get(j))
					result += "!";
				result += atom.symbol.symbol + "(";
				for (int k = 0; k < atom.terms.size(); k++) {
					Term term = atom.terms.get(k);
					if(!termReferenceMap.containsKey(term))
						termReferenceMap.put(term, Character.toString((char)(termReferenceMap.size()+97)));
					result += termReferenceMap.get(term);
					if(k < atom.terms.size()-1)
						result += ",";
				}
				result += ")";
				if(j < clause.atoms.size()-1)
					result += " | ";
			}
			if(i < clauses.size()-1)
				result += " ^ ";
		}
		result += (" :: " + weight.getValue());
		return result;
	}

	public Formula(List<WClause> clauses_, LogDouble weight_, boolean isEvidence_, boolean tautology) {
		clauses = clauses_;
		weight = (weight_);
		isEvidence = (isEvidence_);
		this.tautology = tautology;
	}

}
